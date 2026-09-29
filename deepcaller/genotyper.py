"""
Per-chromosome encoding and genotyping stage.

Run as a subprocess by the command line driver, one chromosome at a time, so
that the operating system reclaims everything the chromosome used as soon as it
is done.

Inside a chromosome the work is a pipeline rather than two separate phases:

    region encoding (process pool)  ->  network inference  ->  Parquet append

Regions are submitted with a bounded number of them in flight, and each result
is genotyped and written out before the next one is collected. Only a couple of
region encodings are ever resident, instead of the whole chromosome plus the
copy that concatenating them would require, and inference overlaps with the
encoding of the following regions instead of waiting for all of them.
"""

import argparse
import os
import time

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import multiprocessing
from collections import deque

import pysam

from . import console
from .config import WINDOW_LENGTH, WINDOW_SIZE
from .features import encode_region
from .variants import VariantWriter, compute_gq_qual

# Regions allowed in flight on top of one per worker. Two is enough to keep the
# pool busy while the driver genotypes a result, and small enough that completed
# encodings cannot pile up in the driver.
QUEUE_SLACK = 2


def parse_args():
    parser = argparse.ArgumentParser(
        description="Encode and genotype the candidate loci of a single chromosome")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--bam", required=True, help="Coordinate-sorted BAM")
    parser.add_argument("--bed", default=None, help="Optional target region file")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--chrom_dp", required=True, type=float,
                        help="Mean depth of this chromosome")
    parser.add_argument("--model_path", required=True, help="Absolute path of the weights")
    parser.add_argument("--ploidy", required=True, type=int)
    parser.add_argument("--min_af", type=float, default=0.10)
    parser.add_argument("--rd_floor", type=int, default=8)
    parser.add_argument("--min_mq", type=int, default=5)
    parser.add_argument("--max_id_len", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=8192)
    parser.add_argument("--cpus", type=int, default=4,
                        help="Worker processes used to encode this chromosome")
    parser.add_argument("--work_dir", required=True,
                        help="Directory receiving <work_dir>/parquet/<chrom>.parquet")
    return parser.parse_args()


def build_regions(chrom_length, num_workers):
    """
    Split a chromosome into regions of work.

    Twice as many regions as workers gives the pool something to rebalance with:
    a worker that finishes early picks up the next region instead of being bound
    to a fixed share of the chromosome. A trailing region that would be very
    short is merged into its predecessor.
    """
    region_count = max(num_workers * 2, 1)
    region_size = max(chrom_length // region_count, 1)

    regions = []
    start = 1
    while start <= chrom_length:
        end = min(start + region_size - 1, chrom_length)
        regions.append((start, end))
        start = end + 1

    if len(regions) >= 2:
        regions[-2] = (regions[-2][0], regions[-1][1])
        regions.pop()

    return regions


def genotype_chromosome(args):
    """Encode, genotype and store every candidate locus of one chromosome."""
    with pysam.FastaFile(args.ref) as fasta:
        chrom_length = fasta.get_reference_length(args.chrom)

    regions = build_regions(chrom_length, args.cpus)
    params = {
        "bam_path": args.bam,
        "fasta_path": args.ref,
        "bed_path": args.bed,
        "min_af": args.min_af,
        "depth_floor": args.rd_floor,
        "ploidy": args.ploidy,
        "max_indel_length": args.max_id_len,
        "min_mapq": args.min_mq,
    }

    parquet_path = os.path.join(args.work_dir, "parquet",
                                "{}.parquet".format(args.chrom))
    writer = VariantWriter(parquet_path)

    console.step(args.chrom, "Scanning alignments and encoding candidate loci "
                             "across {} ...".format(console.bases(chrom_length)))

    started_at = time.perf_counter()
    model = None
    reported = 0

    def report_progress(done):
        """Emit at most one progress line per quarter of the chromosome."""
        nonlocal reported
        if done >= len(regions):
            return  # the completion line below already reports the final state
        quarter = int(4 * done / len(regions))
        if quarter > reported:
            reported = quarter
            console.substep(args.chrom, "progress {:3d}% | {} loci genotyped".format(
                quarter * 25, console.count(writer.rows)))

    def consume(encoding):
        """Genotype one encoded region and append it to the Parquet file."""
        nonlocal model
        if len(encoding) == 0:
            return
        if model is None:
            from .network import configure_threading, load_model
            configure_threading(args.cpus)
            load_timer = console.Timer()
            with load_timer:
                model = load_model(args.model_path, WINDOW_LENGTH, args.ploidy)
            console.substep(args.chrom, "Genotyping model ready ({})".format(
                console.duration(load_timer.seconds)))

        probs, mask = model.predict_in_batches(encoding.features, args.batch_size)
        dosages, gq, qual = compute_gq_qual(
            probs, mask, [len(tokens) for tokens in encoding.alt_tokens])
        writer.write(
            chrom=encoding.chrom, positions=encoding.positions,
            ref_bases=encoding.ref_bases, alt_tokens=encoding.alt_tokens,
            depths=encoding.depths, alt_counts=encoding.alt_counts,
            ref_counts=encoding.ref_counts,
            dosages=dosages, gq=gq, qual=qual,
        )

    try:
        if args.cpus > 1 and len(regions) > 1:
            # The pool is forked before TensorFlow is imported anywhere in this
            # process: forking a runtime that already owns thread pools is not
            # safe. `consume` imports the network lazily, after the fork.
            context = multiprocessing.get_context("fork")
            with context.Pool(processes=min(args.cpus, len(regions))) as pool:
                in_flight = deque()
                submitted = 0
                capacity = min(args.cpus, len(regions)) + QUEUE_SLACK

                while submitted < min(capacity, len(regions)):
                    start, end = regions[submitted]
                    in_flight.append(pool.apply_async(
                        encode_region,
                        (args.chrom, start, end, args.chrom_dp, params)))
                    submitted += 1

                completed = 0
                while in_flight:
                    encoding = in_flight.popleft().get()
                    if submitted < len(regions):
                        start, end = regions[submitted]
                        in_flight.append(pool.apply_async(
                            encode_region,
                            (args.chrom, start, end, args.chrom_dp, params)))
                        submitted += 1
                    consume(encoding)
                    del encoding  # release the region before collecting the next
                    completed += 1
                    report_progress(completed)
        else:
            for completed, (start, end) in enumerate(regions, start=1):
                encoding = encode_region(args.chrom, start, end, args.chrom_dp, params)
                consume(encoding)
                del encoding
                report_progress(completed)
    finally:
        writer.close()

    return writer.rows, time.perf_counter() - started_at


def main():
    args = parse_args()

    if args.ploidy < 1:
        raise SystemExit("[ERROR] --ploidy must be a positive integer")
    if WINDOW_SIZE < 1:
        raise SystemExit("[ERROR] WINDOW_SIZE must be a positive integer")

    rows, elapsed = genotype_chromosome(args)

    if rows == 0:
        console.step(args.chrom, "No candidate locus passed the filters "
                                 "({})".format(console.duration(elapsed)))
    else:
        console.step(args.chrom, "Genotyped {} candidate loci in {}".format(
            console.count(rows), console.duration(elapsed)))


if __name__ == "__main__":
    main()
