"""
Preparation of the reference and of the alignments.

Every function here is a thin, well defined wrapper around samtools, mosdepth or
pysam: index validation, chromosome selection, depth estimation and optional
genome-wide downsampling. They all fail fast with an actionable message rather
than letting a malformed input reach the encoder.
"""

import gzip
import os
import subprocess
import sys

import pandas as pd
import pysam

from . import console


def check_reference(ref_path):
    """Ensure the reference FASTA carries a .fai index, creating it if needed."""
    fai_path = "{}.fai".format(ref_path)
    if os.path.isfile(fai_path):
        console.ok("Found FASTA index: {}".format(fai_path))
        return

    console.info("Generating FASTA index for {} ...".format(ref_path))
    try:
        subprocess.run(["samtools", "faidx", ref_path], check=True,
                       stderr=subprocess.PIPE)
    except FileNotFoundError:
        sys.exit("[ERROR] samtools was not found in PATH")
    except subprocess.CalledProcessError as exc:
        sys.exit("[ERROR] FASTA index generation failed:\n{}".format(
            exc.stderr.decode().strip()))

    if not os.path.isfile(fai_path):
        sys.exit("[ERROR] FASTA index was not created: {}".format(fai_path))
    console.ok("Generated FASTA index: {}".format(fai_path))


def get_chrom_list(ref_path, bed_path=None, user_chroms=None):
    """
    Resolve the chromosomes to process.

    Precedence is BED file, then an explicit --chroms list, then every contig of
    the reference. Names are always validated against the reference.
    """
    try:
        with pysam.FastaFile(ref_path) as fasta:
            ref_chroms = list(fasta.references)
    except Exception as exc:
        sys.exit("[ERROR] Failed to read {}\nReason: {}".format(ref_path, exc))

    if not ref_chroms:
        sys.exit("[ERROR] Empty reference: {}".format(ref_path))

    if bed_path:
        chrom_list, seen = [], set()
        try:
            opener = gzip.open if bed_path.endswith(".gz") else open
            with opener(bed_path, "rt") as handle:
                for line in handle:
                    if not line.strip() or line.startswith(("track", "browser", "#")):
                        continue
                    chrom = line.split("\t")[0]
                    if chrom not in seen:
                        chrom_list.append(chrom)
                        seen.add(chrom)
        except Exception as exc:
            sys.exit("[ERROR] Failed to read {}\nReason: {}".format(bed_path, exc))

        if not chrom_list:
            sys.exit("[ERROR] No valid chromosome found in BED file: {}".format(bed_path))

        unknown = set(chrom_list) - set(ref_chroms)
        if unknown:
            sys.exit("[ERROR] Chromosomes absent from the reference: {}\n"
                     "  Valid options: {}".format(sorted(unknown), " ".join(ref_chroms)))

        console.info("Chromosomes taken from the BED file: {}".format(", ".join(chrom_list)))
        return chrom_list

    if user_chroms:
        unknown = set(user_chroms) - set(ref_chroms)
        if unknown:
            sys.exit("[ERROR] Chromosomes absent from the reference: {}\n"
                     "  Valid options: {}".format(sorted(unknown), " ".join(ref_chroms)))
        console.info("User-specified chromosomes: {}".format(", ".join(user_chroms)))
        return list(user_chroms)

    console.info("No BED file or chromosome list given, using every reference contig")
    return ref_chroms


def check_alignments(bam_path, num_threads):
    """
    Verify that the BAM is coordinate-sorted and indexed.

    Sorting is not performed here: for whole-genome inputs it is an expensive
    operation that the user should schedule explicitly.
    """
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            sort_order = bam.header.get("HD", {}).get("SO", "unspecified")
    except Exception as exc:
        raise RuntimeError("Failed to read BAM file: {}\nReason: {}".format(bam_path, exc))

    if sort_order != "coordinate":
        sys.exit(
            "[ERROR] BAM file is not coordinate-sorted: {bam}\n"
            "Please sort it first:\n"
            "  samtools sort -@ {threads} -o sorted.bam {bam}\n"
            "then re-run with the sorted file.".format(bam=bam_path, threads=num_threads)
        )

    bai_path = "{}.bai".format(bam_path)
    if os.path.exists(bai_path):
        console.ok("Found BAM index: {}".format(bai_path))
        return bam_path

    console.info("Generating index for {} ...".format(bam_path))
    try:
        subprocess.run(["samtools", "index", "-@", str(num_threads), bam_path],
                       check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as exc:
        sys.exit("[ERROR] Index generation failed:\n{}\n"
                 "You can create it manually with:\n"
                 "  samtools index -@ {} {}".format(
                     exc.stderr.decode().strip(), num_threads, bam_path))
    console.ok("Created BAM index: {}".format(bai_path))
    return bam_path


def compute_depth(bam_path, work_dir, num_threads, prefix="genome"):
    """
    Estimate mean coverage with mosdepth.

    Returns:
        chrom_depth_map: {chromosome: (length, mean depth)}
        genome_mean_depth: mean depth of the 'total' row of the summary
    """
    output_prefix = os.path.join(work_dir, prefix)
    command = ["mosdepth", "-t", str(min(num_threads, 4)), "-x", "-n",
               output_prefix, bam_path]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        raise RuntimeError("mosdepth was not found in PATH")
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("mosdepth failed on {}\nCMD: {}\nSTDERR:\n{}".format(
            bam_path, " ".join(command), exc.stderr))

    summary_path = "{}.mosdepth.summary.txt".format(output_prefix)
    if not os.path.exists(summary_path):
        raise RuntimeError("mosdepth summary was not produced: {}".format(summary_path))

    summary = pd.read_csv(summary_path, sep="\t")

    total_row = summary[summary["chrom"] == "total"]
    if total_row.empty:
        raise RuntimeError("No 'total' row in mosdepth summary: {}".format(summary_path))
    genome_mean_depth = float(total_row["mean"].iloc[0])

    chrom_depth_map = {
        row.chrom: (int(row.length), float(row.mean))
        for row in summary.itertuples()
        if row.chrom != "total"
    }
    if not chrom_depth_map:
        raise RuntimeError("No per-chromosome depth in mosdepth summary: {}".format(
            summary_path))

    return chrom_depth_map, genome_mean_depth


def downsample_genome(bam_path, work_dir, fraction, seed, num_threads):
    """
    Subsample the whole BAM in a single samtools pass.

    Sampling genome-wide rather than per chromosome keeps read pairs and the
    relative coverage between chromosomes consistent.

    Returns:
        (path of the downsampled BAM, chromosome depth map recomputed on it)
    """
    if not 0 < fraction < 1:
        raise RuntimeError("Downsampling fraction must lie in (0, 1), got {}".format(fraction))

    output_bam = os.path.join(work_dir, "downsampled.bam")
    # samtools encodes seed and fraction in a single SEED.FRACTION argument.
    sampling_argument = "{}.{}".format(seed, "{:.3f}".format(fraction).split(".")[1])

    command = ["samtools", "view", "-b", "-@", str(num_threads),
               "-s", sampling_argument, "-o", output_bam, bam_path]
    try:
        subprocess.run(command, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("Genome-wide downsampling failed\nReason: {}".format(
            exc.stderr.decode().strip()))

    try:
        pysam.index(output_bam)
    except Exception as exc:
        raise RuntimeError("Indexing failed for {}\nReason: {}".format(output_bam, exc))

    chrom_depth_map, _ = compute_depth(output_bam, work_dir, num_threads,
                                       prefix="downsampled")
    return output_bam, chrom_depth_map
