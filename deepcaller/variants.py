"""
Variant records: quality scores, intermediate storage and VCF assembly.

Genotyped loci are streamed into a per-chromosome Parquet file as soon as they
are scored, and streamed back out row group by row group when the VCF is built.
Neither stage ever holds a whole chromosome in memory.
"""

import os

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

# Explicit schema: it keeps the files written by different chromosomes
# compatible and avoids type inference surprises on empty or ragged chunks.
VARIANT_SCHEMA = pa.schema([
    ("chrom", pa.string()),
    ("pos", pa.int64()),
    ("ref", pa.string()),
    ("alt", pa.list_(pa.string())),       # candidate allele tokens, ranked
    ("rd", pa.int32()),                   # read depth at the locus
    ("alt_num", pa.list_(pa.int32())),    # support of each candidate allele
    ("ref_num", pa.int32()),              # reads supporting the reference
    ("dosage", pa.list_(pa.int8())),      # predicted copy number per slot
    ("gq", pa.float32()),
    ("qual", pa.float32()),
])

VCF_READ_BATCH = 65536


# --------------------------------------------------------------------------- #
# Quality scores
# --------------------------------------------------------------------------- #
def compute_gq_qual(probs, mask=None, num_alts=None):
    """
    Derive the genotype call and its quality scores from the network output.

    Args:
        probs: (n, ploidy, ploidy + 1) softmax distributions.
        mask, num_alts: accepted for call compatibility and deliberately unused;
            the joint probability is taken over all slots, masked or not.

    Returns:
        dosages: (n, ploidy) int, the argmax copy number of each slot.
        gq: joint probability of the called combination, on the Phred scale.
            The slots are not strictly independent, so this is an approximation.
        qual: -10 log10 P(all slots are reference), i.e. the confidence that the
            locus is not homozygous reference.
    """
    n, ploidy, _ = probs.shape
    dosages = np.argmax(probs, axis=-1)

    row_index = np.arange(n)[:, None]
    slot_index = np.arange(ploidy)[None, :]
    picked = probs[row_index, slot_index, dosages]
    joint = np.clip(np.prod(picked, axis=1), 1e-8, 1 - 1e-8)
    gq = (-10.0 * np.log10(1.0 - joint)).astype(np.float32)

    all_reference = np.clip(np.prod(probs[:, :, 0], axis=1), 1e-8, 1 - 1e-8)
    qual = (-10.0 * np.log10(all_reference)).astype(np.float32)

    return dosages, gq, qual


# --------------------------------------------------------------------------- #
# Intermediate storage
# --------------------------------------------------------------------------- #
class VariantWriter:
    """
    Incremental Parquet writer for one chromosome.

    The underlying file is only created on the first non-empty batch, so a
    chromosome without candidates leaves no artefact behind.
    """

    def __init__(self, path):
        self.path = path
        self.rows = 0
        self._writer = None

    def write(self, chrom, positions, ref_bases, alt_tokens, depths,
              alt_counts, ref_counts, dosages, gq, qual):
        """Append one batch of genotyped loci; all sequences share one order."""
        rows = len(positions)
        if rows == 0:
            return

        table = pa.Table.from_arrays(
            [
                pa.array([chrom] * rows, type=pa.string()),
                pa.array(positions, type=pa.int64()),
                pa.array(ref_bases, type=pa.string()),
                pa.array(alt_tokens, type=pa.list_(pa.string())),
                pa.array(depths, type=pa.int32()),
                pa.array(alt_counts, type=pa.list_(pa.int32())),
                pa.array(ref_counts, type=pa.int32()),
                pa.array(np.asarray(dosages, dtype=np.int8).tolist(),
                         type=pa.list_(pa.int8())),
                pa.array(np.asarray(gq, dtype=np.float32), type=pa.float32()),
                pa.array(np.asarray(qual, dtype=np.float32), type=pa.float32()),
            ],
            schema=VARIANT_SCHEMA,
        )

        if self._writer is None:
            os.makedirs(os.path.dirname(self.path), exist_ok=True)
            self._writer = pq.ParquetWriter(self.path, VARIANT_SCHEMA,
                                            compression="snappy")
        self._writer.write_table(table)
        self.rows += rows

    def close(self):
        if self._writer is not None:
            self._writer.close()
            self._writer = None


# --------------------------------------------------------------------------- #
# VCF records
# --------------------------------------------------------------------------- #
def build_ref_alt(ref_base, called_tokens):
    """
    Turn internal allele tokens into a VCF REF/ALT pair sharing one REF span.

    Internal tokens are a single base for substitutions, '+BASES' for insertions
    and '-BASES' for deletions. The longest deletion of the locus defines a
    common suffix that is appended to REF and to every ALT, so that a locus
    carrying both a substitution and a deletion still yields consistent lengths.

    This assumes that when several deletions are called at one locus, the shorter
    ones are prefixes of the longest. That holds for the overlapping deletions
    seen in real data, but it is not guaranteed mathematically.
    """
    deletions = [token for token in called_tokens if token.startswith("-")]
    suffix = max(deletions, key=len)[1:] if deletions else ""

    vcf_ref = ref_base + suffix
    vcf_alts = []
    for token in called_tokens:
        if token.startswith("+"):
            vcf_alts.append(ref_base + token[1:] + suffix)
        elif token.startswith("-"):
            vcf_alts.append(ref_base + suffix[len(token) - 1:])
        else:
            vcf_alts.append(token + suffix)
    return vcf_ref, vcf_alts


def format_vcf_record(chrom, pos, ref_base, alt_tokens, dosage, ploidy,
                      rd, alt_num, ref_num, gq, qual):
    """
    Render one genotyped locus as a VCF data line.

    Only alleles that received at least one copy reach ALT and GT. A locus where
    every slot was called reference is still reported, with ALT '.' and the
    RefCall filter, which keeps the output usable for downstream concordance
    analysis.
    """
    num_alts = len(alt_tokens)
    called = [(j, alt_tokens[j], int(dosage[j]))
              for j in range(num_alts) if int(dosage[j]) > 0]
    alt_copies = sum(int(dosage[j]) for j in range(num_alts))
    ref_dosage = max(ploidy - alt_copies, 0)

    if not called:
        ref_out = ref_base
        alt_field = "."
        gt_field = "/".join(["0"] * ploidy)
        filter_field = "RefCall"
        ad_field = str(int(ref_num))
        af_field = "0.0"
    else:
        ref_out, alt_out = build_ref_alt(ref_base, [token for _, token, _ in called])
        alt_field = ",".join(alt_out)

        # GT lists ref_dosage zeros followed by each called allele repeated as
        # often as its copy number; the genotype is unphased.
        gt_alleles = [0] * ref_dosage
        for allele_index, (_, _, copies) in enumerate(called, start=1):
            gt_alleles += [allele_index] * copies
        gt_field = "/".join(str(a) for a in sorted(gt_alleles))
        filter_field = "PASS"

        ad_field = ",".join([str(int(ref_num))] +
                            [str(int(alt_num[j])) for j, _, _ in called])
        af_field = ",".join(str(round(int(alt_num[j]) / max(int(rd), 1), 4))
                            for j, _, _ in called)

    return "\t".join([
        chrom, str(pos), ".", ref_out, alt_field,
        "{:.1f}".format(qual), filter_field,
        "DP={}".format(int(rd)),
        "GT:GQ:DP:AD:AF",
        "{}:{:.1f}:{}:{}:{}".format(gt_field, gq, int(rd), ad_field, af_field),
    ])


def _vcf_header(chrom_lengths, sample_name):
    lines = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth at position">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt1,alt2,...)">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">',
    ]
    for chrom, length in chrom_lengths.items():
        lines.append("##contig=<ID={},length={}>".format(chrom, length))
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name)
    return "\n".join(lines)


def generate_vcf(ref_path, chrom_list, out_vcf_path, ploidy, work_dir, sample_name):
    """
    Concatenate the per-chromosome Parquet files into a bgzip-compressed VCF.

    Row groups are streamed one at a time, so the memory footprint is bounded by
    VCF_READ_BATCH records regardless of how many variants a chromosome holds.

    Returns:
        (path of the indexed VCF, number of records written)
    """
    with pysam.FastaFile(ref_path) as fasta:
        chrom_lengths = {c: fasta.get_reference_length(c) for c in chrom_list}

    parquet_dir = os.path.join(work_dir, "parquet")
    plain_vcf = out_vcf_path[:-3] if out_vcf_path.endswith(".gz") else out_vcf_path
    written = 0

    with open(plain_vcf, "w") as out:
        out.write(_vcf_header(chrom_lengths, sample_name) + "\n")

        for chrom in chrom_list:
            path = os.path.join(parquet_dir, "{}.parquet".format(chrom))
            if not os.path.exists(path):
                continue  # no candidate locus survived the filters here

            parquet_file = pq.ParquetFile(path)
            for batch in parquet_file.iter_batches(batch_size=VCF_READ_BATCH):
                columns = {name: batch.column(name).to_pylist()
                           for name in VARIANT_SCHEMA.names}
                out.writelines(
                    format_vcf_record(
                        chrom=columns["chrom"][i], pos=columns["pos"][i],
                        ref_base=columns["ref"][i], alt_tokens=columns["alt"][i],
                        dosage=columns["dosage"][i], ploidy=ploidy,
                        rd=columns["rd"][i], alt_num=columns["alt_num"][i],
                        ref_num=columns["ref_num"][i],
                        gq=columns["gq"][i], qual=columns["qual"][i],
                    ) + "\n"
                    for i in range(batch.num_rows)
                )
                written += batch.num_rows

    return pysam.tabix_index(plain_vcf, preset="vcf", force=True), written
