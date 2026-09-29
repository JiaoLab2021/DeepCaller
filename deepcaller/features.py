"""
Candidate discovery and pileup feature encoding.

A region of a chromosome is streamed once through a pysam pileup. A sliding
window of WINDOW_LENGTH columns is kept in memory; whenever the centre column
passes the candidate filters, the reads covering the window are summarised into
a (ploidy, WINDOW_LENGTH, FEATURE_DIM) tensor that the network consumes.

Two properties drive the implementation:

* Reads are grouped per candidate allele. Each alternate allele owns one slot
  ("head") of the tensor, so the network sees one independent view per allele
  and can predict its copy number.
* The window slides by one column at a time, which means every column is part of
  up to WINDOW_LENGTH windows. Anything that depends on the column alone - the
  token class of each read, its mapping quality, its identity - is therefore
  computed once when the column enters the window and reused afterwards.
"""

import bisect
import gzip
import re
from collections import defaultdict, deque
from typing import List, NamedTuple

import numpy as np
import pysam

from .config import (ENCODING_BLOCK_SIZE, ENCODING_DTYPE, FEATURE_DIM,
                     MANDATORY_GROUPS, MAPQ_SCALE, PILEUP_FLAG_FILTER,
                     REFERENCE_PADDING, TOKEN_CLASSES, WINDOW_LENGTH, WINDOW_SIZE)

# Token classes: A, C, G, T, insertion, deletion, gap or unusable base.
_BASE_CATEGORY = {"A": 0, "C": 1, "G": 2, "T": 3}
_INSERTION_CATEGORY = 4
_DELETION_CATEGORY = 5
_GAP_CATEGORY = 6

_INDEL_TOKEN_RE = re.compile(r"^(\d+)([ACGTNacgtn]+)")

# Fast path for the overwhelmingly common single-character pileup token; any
# other single character (gap, N, ...) carries no allele and maps to None.
_SINGLE_BASE = {"A": "A", "C": "C", "G": "G", "T": "T",
                "a": "A", "c": "C", "g": "G", "t": "T"}

# Reference bases are normalised to upper case; anything that is not a plain
# nucleotide becomes '*' so that it can never be called as a candidate.
_REFERENCE_TABLE = {}
for _code in range(256):
    _char = chr(_code)
    _REFERENCE_TABLE[_code] = _char.upper() if _char in "ACGTacgt" else "*"


# --------------------------------------------------------------------------- #
# Target regions
# --------------------------------------------------------------------------- #
class BedIndex:
    """
    Interval index over a BED file, grouped by chromosome.

    Overlapping intervals are merged at load time so that a membership test is a
    single binary search, which matters because it runs once per pileup column.
    """

    def __init__(self, bed_path):
        raw = defaultdict(list)
        opener = gzip.open if bed_path.endswith(".gz") else open
        with opener(bed_path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("track", "browser", "#")):
                    continue
                chrom, start, end, *_ = line.split("\t")
                # BED is 0-based half-open; the pipeline works in 1-based coordinates.
                start, end = int(start) + 1, int(end) + 1
                if end > start:
                    raw[chrom].append((start, end))

        self._index = {}
        for chrom, intervals in raw.items():
            intervals.sort()
            merged = []
            for start, end in intervals:
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            self._index[chrom] = (
                [start for start, _ in merged],
                [end for _, end in merged],
            )

    def contains(self, chrom, position):
        """Return True when a 1-based position falls inside a target interval."""
        entry = self._index.get(chrom)
        if entry is None:
            return False
        starts, ends = entry
        i = bisect.bisect_right(starts, position) - 1
        return i >= 0 and position < ends[i]


class ReferenceSlice:
    """
    A contiguous slice of the reference sequence held as a plain string.

    Storing the slice as a string rather than a {position: base} dictionary keeps
    the memory cost at one byte per base instead of the ~100 bytes per entry that
    a boxed dictionary would need over a region of several megabases.
    """

    __slots__ = ("_sequence", "_first", "_last")

    def __init__(self, fasta, chrom, first, last):
        sequence = fasta.fetch(chrom, first - 1, last)
        self._sequence = sequence.translate(_REFERENCE_TABLE)
        self._first = first
        self._last = first + len(self._sequence) - 1

    def base(self, position):
        """Reference base at a 1-based position, or '*' outside the slice."""
        if self._first <= position <= self._last:
            return self._sequence[position - self._first]
        return "*"

    def segment(self, first, last):
        """Inclusive 1-based substring, or None when it leaves the slice."""
        if first < self._first or last > self._last:
            return None
        return self._sequence[first - self._first:last - self._first + 1]


# --------------------------------------------------------------------------- #
# Pileup tokens
# --------------------------------------------------------------------------- #
def group_labels(ploidy):
    """Names of the allele slots of the feature tensor: alt1 .. altP."""
    return ["alt{}".format(i) for i in range(1, ploidy + 1)]


def parse_indel_token(text):
    """
    Split a pysam indel suffix such as '3ACG' or '2AT' into (length, bases).

    pysam appends the indel to the aligned base, e.g. 'A+3ACG'; this helper is
    fed the part that follows the sign.
    """
    match = _INDEL_TOKEN_RE.match(text)
    length = int(match.group(1))
    return length, match.group(2).upper()[:length]


def token_category(token):
    """
    Map a pileup token to its class index, or None when it carries no signal.

    Insertions and deletions are recognised before the leading base is examined
    because pysam reports them as a base followed by the event.
    """
    if not token or token == "*":
        return _GAP_CATEGORY
    if "+" in token:
        return _INSERTION_CATEGORY
    if "-" in token:
        return _DELETION_CATEGORY
    return _BASE_CATEGORY.get(token[0].upper())


def _normalise_token(token, reference, position):
    """
    Convert a raw pileup token into the internal allele representation:
    a single base for substitutions, '+BASES' for insertions and '-BASES' for
    deletions, where BASES are the deleted reference bases. Returns None when
    the token cannot be represented (e.g. a deletion running past the slice).
    """
    if not token or token == "*":
        return None
    if "+" in token:
        _, inserted = parse_indel_token(token.split("+")[1])
        return "+" + inserted
    if "-" in token:
        length, _ = parse_indel_token(token.split("-")[1])
        deleted = reference.segment(position + 1, position + length)
        return None if deleted is None else "-" + deleted
    base = token[0].upper()
    return base if base in _BASE_CATEGORY else None


def select_candidate_alleles(tokens, ref_base, depth, position, reference,
                             min_af, max_indel_length, max_alleles):
    """
    Rank the non-reference alleles observed at the centre column.

    Alleles are counted, filtered on indel length and on the minimum allele
    fraction, then truncated to the number of slots the tensor provides.
    Returns a list of (allele token, supporting read count) sorted by support.
    """
    counter = {}
    for token in tokens:
        allele = (_SINGLE_BASE.get(token) if len(token) == 1
                  else _normalise_token(token, reference, position))
        if allele is None or allele == ref_base:
            continue
        counter[allele] = counter.get(allele, 0) + 1

    candidates = [
        (allele, support) for allele, support in counter.items()
        if allele[0] not in "+-" or len(allele) - 1 <= max_indel_length
    ]
    passing = [item for item in candidates if item[1] / depth >= min_af]
    passing.sort(key=lambda item: item[1], reverse=True)
    return passing[:max_alleles]


def assign_read_groups(tokens, keys, ref_base, position, reference, alleles):
    """
    Label every read covering the centre column as 'ref', 'altN' or 'other'.

    The labels drive which slot of the feature tensor a read contributes to for
    the whole window, not only for the centre column.
    """
    slot_of_allele = {allele: i + 1 for i, (allele, _) in enumerate(alleles)}
    groups = {}
    for token, key in zip(tokens, keys):
        allele = (_SINGLE_BASE.get(token) if len(token) == 1
                  else _normalise_token(token, reference, position))
        if allele == ref_base:
            groups[key] = "ref"
        elif allele in slot_of_allele:
            groups[key] = "alt{}".format(slot_of_allele[allele])
        else:
            groups[key] = "other"
    return groups


# --------------------------------------------------------------------------- #
# Feature tensor assembly
# --------------------------------------------------------------------------- #
class FeatureBuffer:
    """
    Append-only store for per-locus feature tensors.

    Loci are written into pre-allocated blocks and the blocks are consolidated
    into one contiguous array at the end, each block being released as soon as it
    has been copied. Peak memory is therefore the size of the final array plus a
    single block, instead of the twofold peak of `np.array(list_of_tensors)`.
    """

    __slots__ = ("_shape", "_blocks", "_fill", "_block_size")

    def __init__(self, locus_shape, block_size=ENCODING_BLOCK_SIZE):
        self._shape = tuple(locus_shape)
        self._block_size = block_size
        self._blocks = []
        self._fill = 0

    def __len__(self):
        return (len(self._blocks) - 1) * self._block_size + self._fill if self._blocks else 0

    def next_slot(self):
        """Return a zeroed view of the next free locus, ready to be filled."""
        if not self._blocks or self._fill == self._block_size:
            self._blocks.append(
                np.zeros((self._block_size,) + self._shape, dtype=ENCODING_DTYPE))
            self._fill = 0
        slot = self._blocks[-1][self._fill]
        slot[...] = 0
        self._fill += 1
        return slot

    def consolidate(self):
        """Return one contiguous array and drop every intermediate block."""
        total = len(self)
        result = np.empty((total,) + self._shape, dtype=ENCODING_DTYPE)
        offset = 0
        for index, block in enumerate(self._blocks):
            rows = self._fill if index == len(self._blocks) - 1 else self._block_size
            result[offset:offset + rows] = block[:rows]
            offset += rows
            self._blocks[index] = None  # release as we go, not at the end
        self._blocks = []
        self._fill = 0
        return result


class PileupWindow:
    """
    Sliding window of pileup columns.

    For every column the window stores the raw pileup tokens, the identity and
    the mapping quality of each read. The per-read view consumed by the encoder -
    (token class, mapping quality, read key) triples - is built on first use and
    cached: a column is shared by up to WINDOW_LENGTH windows, but most columns
    never belong to a window that yields a candidate, and those must not pay for
    the conversion.
    """

    __slots__ = ("positions", "ref_bases", "tokens", "keys", "mapqs", "depths", "_views")

    def __init__(self, length):
        self.positions = deque(maxlen=length)
        self.ref_bases = deque(maxlen=length)
        self.tokens = deque(maxlen=length)
        self.keys = deque(maxlen=length)
        self.mapqs = deque(maxlen=length)
        self.depths = deque(maxlen=length)
        self._views = deque(maxlen=length)

    def __len__(self):
        return len(self.positions)

    def push(self, position, ref_base, tokens, keys, mapqs, depth):
        """Append a column; the oldest one leaves the window automatically."""
        self.positions.append(position)
        self.ref_bases.append(ref_base)
        self.tokens.append(tokens)
        self.keys.append(keys)
        self.mapqs.append(mapqs)
        self.depths.append(depth)
        self._views.append(None)

    def view(self, index):
        """(token class, mapping quality, read key) triples of one column."""
        cached = self._views[index]
        if cached is None:
            cached = list(zip([token_category(token) for token in self.tokens[index]],
                              self.mapqs[index], self.keys[index]))
            self._views[index] = cached
        return cached


class WindowScratch:
    """
    Buffers reused by every locus of a region.

    Read counts and mapping-quality sums are accumulated in a flat Python list
    of integers and converted to float32 once per locus, rather than per window
    position: a single conversion of the whole window costs far less than the
    per-column numpy calls it replaces, and integer accumulation is exact.
    """

    __slots__ = ("flat", "blank", "stats", "flat_view", "counts_divisor",
                 "mapq_divisor", "stride", "ploidy")

    def __init__(self, window_length, ploidy):
        self.ploidy = ploidy
        self.stride = 2 * TOKEN_CLASSES
        size = window_length * ploidy * self.stride
        self.flat = [0] * size
        self.blank = [0] * size
        self.stats = np.zeros((window_length, ploidy, self.stride), dtype=np.float32)
        self.flat_view = self.stats.reshape(-1)
        self.counts_divisor = np.zeros(window_length, dtype=np.float32)
        self.mapq_divisor = np.zeros(window_length, dtype=np.float32)


def encode_window(target, window, groups, chrom_depth, ploidy, scratch):
    """
    Fill one (ploidy, WINDOW_LENGTH, FEATURE_DIM) tensor in place.

    Every read of every column is visited once and dispatched to the slot of the
    allele it supports, so the inner loop does not repeat itself per allele; the
    whole window is then normalised with a handful of vectorised operations.

    A slot stays entirely zero when no read supports it, unless it is listed in
    MANDATORY_GROUPS; the network reads that all-zero pattern as a masked slot.
    """
    labels = group_labels(ploidy)
    slot_of_label = {label: index for index, label in enumerate(labels)}

    slot_of_key = {}
    occupied = [False] * ploidy
    for key, label in groups.items():
        slot = slot_of_label.get(label)
        if slot is not None:
            slot_of_key[key] = slot
            occupied[slot] = True

    active = [i for i in range(ploidy) if occupied[i] or labels[i] in MANDATORY_GROUPS]
    if not active:
        return

    # Channel 0: the normalised depth profile, shared by every active slot.
    depths = window.depths
    depth_profile = np.asarray(depths, dtype=np.float64) / chrom_depth
    for slot in active:
        target[slot, :, 0] = depth_profile

    if not slot_of_key:
        return

    counts_divisor = scratch.counts_divisor
    mapq_divisor = scratch.mapq_divisor
    counts_divisor[:] = depths
    np.maximum(counts_divisor, 1, out=counts_divisor)
    mapq_divisor[:] = counts_divisor
    mapq_divisor *= MAPQ_SCALE

    flat = scratch.flat
    flat[:] = scratch.blank
    stride = scratch.stride
    column_view = window.view
    slot_lookup = slot_of_key.get

    for index in range(len(depths)):
        if depths[index] <= 0:
            continue
        base = index * ploidy * stride
        for category, mapq, key in column_view(index):
            if category is None:
                continue
            slot = slot_lookup(key)
            if slot is None:
                continue
            offset = base + slot * stride + category
            flat[offset] += 1
            flat[offset + TOKEN_CLASSES] += mapq

    statistics = scratch.stats
    scratch.flat_view[:] = flat
    statistics[:, :, :TOKEN_CLASSES] /= counts_divisor[:, None, None]
    statistics[:, :, TOKEN_CLASSES:] /= mapq_divisor[:, None, None]
    target[:, :, 1:FEATURE_DIM] = statistics.transpose(1, 0, 2)


# --------------------------------------------------------------------------- #
# Region encoding
# --------------------------------------------------------------------------- #
class RegionEncoding(NamedTuple):
    """Everything one worker returns for one region of a chromosome."""

    features: np.ndarray       # (n, ploidy, WINDOW_LENGTH, FEATURE_DIM) float16
    chrom: str
    positions: List[int]       # 1-based candidate positions
    ref_bases: List[str]
    alt_tokens: List[List[str]]
    depths: List[int]
    alt_counts: List[List[int]]
    ref_counts: List[int]

    def __len__(self):
        return len(self.positions)


def encode_region(chrom, start, end, chrom_depth, params):
    """
    Discover and encode every candidate locus in [start, end] (1-based inclusive).

    The function is the unit of work handed to the process pool, so it opens its
    own file handles and returns a self-contained RegionEncoding.
    """
    ploidy = params["ploidy"]
    min_af = params["min_af"]
    min_mapq = params["min_mapq"]
    max_indel_length = params["max_indel_length"]
    min_depth = max(round(chrom_depth / ploidy), params["depth_floor"])

    bed_index = BedIndex(params["bed_path"]) if params.get("bed_path") else None

    fasta = pysam.FastaFile(params["fasta_path"])
    chrom_length = fasta.get_reference_length(chrom)
    alignments = pysam.AlignmentFile(
        params["bam_path"], "rb", reference_filename=params["fasta_path"])

    reference = ReferenceSlice(
        fasta, chrom,
        max(1, start - REFERENCE_PADDING),
        min(chrom_length, end + REFERENCE_PADDING),
    )

    centre = WINDOW_SIZE
    window = PileupWindow(WINDOW_LENGTH)

    buffer = FeatureBuffer((ploidy, WINDOW_LENGTH, FEATURE_DIM))
    scratch = WindowScratch(WINDOW_LENGTH, ploidy)
    out_positions, out_ref, out_depth = [], [], []
    out_alt, out_alt_counts, out_ref_counts = [], [], []

    for column in alignments.pileup(
            chrom,
            max(0, start - WINDOW_SIZE - 1),
            min(end + WINDOW_SIZE, chrom_length),
            min_base_quality=0,
            min_mapping_quality=min_mapq,
            flag_filter=PILEUP_FLAG_FILTER,
            truncate=True,
            multiple_iterators=False):

        pileup_reads = column.pileups
        tokens = column.get_query_sequences(
            mark_matches=False, mark_ends=False, add_indels=True)

        # Read identity and mapping quality must be extracted while the column
        # is current: pysam invalidates the pileup records once it advances.
        keys, mapqs = [], []
        for pileup_read in pileup_reads:
            record = pileup_read.alignment
            keys.append((record.query_name, record.is_read1))
            mapqs.append(record.mapping_quality)

        position = column.pos + 1
        window.push(position, reference.base(position), tokens, keys, mapqs,
                    column.get_num_aligned())

        if len(window) < WINDOW_LENGTH:
            continue

        ref_base = window.ref_bases[centre]
        depth = window.depths[centre]
        centre_position = window.positions[centre]

        if ref_base not in _BASE_CATEGORY or depth < min_depth:
            continue
        if bed_index is not None and not bed_index.contains(chrom, centre_position):
            continue

        centre_tokens = window.tokens[centre]
        centre_keys = window.keys[centre]

        alleles = select_candidate_alleles(
            centre_tokens, ref_base, depth, centre_position, reference,
            min_af, max_indel_length, ploidy)
        if not alleles:
            continue

        groups = assign_read_groups(
            centre_tokens, centre_keys, ref_base, centre_position, reference, alleles)

        encode_window(buffer.next_slot(), window, groups, chrom_depth, ploidy, scratch)

        out_positions.append(centre_position)
        out_ref.append(ref_base)
        out_depth.append(depth)
        out_alt.append([allele for allele, _ in alleles])
        out_alt_counts.append([support for _, support in alleles])
        out_ref_counts.append(sum(1 for label in groups.values() if label == "ref"))

    alignments.close()
    fasta.close()

    return RegionEncoding(
        features=buffer.consolidate(),
        chrom=chrom,
        positions=out_positions,
        ref_bases=out_ref,
        alt_tokens=out_alt,
        depths=out_depth,
        alt_counts=out_alt_counts,
        ref_counts=out_ref_counts,
    )
