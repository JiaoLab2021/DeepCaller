#!/usr/bin/env python3
"""
DeepCaller command line driver.

The driver validates the inputs, prepares the alignments, then hands each
chromosome to a dedicated subprocess and finally merges the per-chromosome
results into a single VCF. Running the chromosomes out of process is deliberate:
peak memory is then set by the largest chromosome alone, because the operating
system reclaims everything when the subprocess exits.
"""

import os

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import argparse
import multiprocessing
import shutil
import subprocess
import sys
import time
import uuid
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

from . import __version__, console
from .alignment import (check_alignments, check_reference, compute_depth,
                        downsample_genome, get_chrom_list)
from .config import (DOWNSAMPLE_TARGET_DEPTH, MODEL_REGISTRY, MODELS_ROOT,
                     SUPPORTED_PLOIDY, default_species, describe_species_options,
                     resolve_weights, species_choices)
from .variants import generate_vcf

GENOTYPER_MODULE = "deepcaller.genotyper"

try:  # optional, purely cosmetic in `ps` output
    from setproctitle import setproctitle

    setproctitle("DeepCaller " + " ".join(sys.argv[1:]))
except ImportError:
    pass


# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
def cpu_count_type(value):
    """argparse type accepting a positive core count or -1 for 'all cores'."""
    count = int(value)
    if count == -1:
        return None
    if count <= 0:
        raise argparse.ArgumentTypeError(
            "must be -1 or a positive integer, got {}".format(count))
    return count


def build_parser():
    parser = argparse.ArgumentParser(
        prog="deepcaller",
        description="DeepCaller {} - small variant discovery and genotyping "
                    "for polyploid genomes".format(__version__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version",
                        version="DeepCaller {}".format(__version__))

    required = parser.add_argument_group("Required arguments")
    required.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    required.add_argument("-b", "--bam", required=True, help="Input BAM file")
    required.add_argument("-p", "--ploidy", required=True, type=int,
                          choices=list(SUPPORTED_PLOIDY), help="Ploidy level")

    io_group = parser.add_argument_group("Input/output configuration")
    io_group.add_argument("-c", "--chroms", nargs="+", default=None,
                          help="Chromosomes to include")
    io_group.add_argument("-o", "--out", default="output.vcf", help="Output VCF file")
    io_group.add_argument("-l", "--bed", default=None,
                          help="Optional BED file. If provided, --chroms is ignored")
    io_group.add_argument("--work_dir", default=None,
                          help="Temp directory; defaults to a unique directory "
                               "under the working directory, removed on success")

    processing = parser.add_argument_group("Processing options")
    processing.add_argument("-t", "--cpus", type=cpu_count_type, default=24,
                            help="Number of CPU cores. Use -1 for all available")
    processing.add_argument("--downsample", action="store_true",
                            help="Downsample the BAM to a target depth "
                                 "(ploidy 4: 50X, ploidy 6: 80X) when the "
                                 "genome-wide depth exceeds it")
    processing.add_argument("--species", default=None, help=describe_species_options())
    processing.add_argument("-v", "--min_af", type=float, default=0.10,
                            help="Minimum allele fraction of a candidate allele")
    processing.add_argument("-d", "--rd_floor", type=int, default=8,
                            help="Minimum read depth of a candidate locus")
    processing.add_argument("--sample", default="SAMPLE", help="Sample name written to the VCF")
    processing.add_argument("--min_mq", type=int, default=5,
                            help="Minimum mapping quality kept during pileup")
    processing.add_argument("--max_id_len", type=int, default=50,
                            help="Maximum indel length considered as a candidate allele")
    processing.add_argument("--batch_size", type=int, default=8192,
                            help="Model inference batch size")
    processing.add_argument("--seed", type=int, default=42,
                            help="Random seed used by samtools view -s")
    processing.add_argument("--keep_tmp", action="store_true",
                            help="Keep the temp directory for debugging")
    return parser


# --------------------------------------------------------------------------- #
# Stages
# --------------------------------------------------------------------------- #
def validate_inputs(args):
    """Reject anything that would only fail much later in the run."""
    console.banner("Input validation")

    if not os.path.isfile(args.ref):
        sys.exit("[ERROR] Reference FASTA not found: {}".format(args.ref))
    if not os.path.isfile(args.bam):
        sys.exit("[ERROR] BAM file not found: {}".format(args.bam))
    if args.bed is not None and not os.path.isfile(args.bed):
        sys.exit("[ERROR] BED file not found: {}".format(args.bed))

    output_dir = os.path.dirname(os.path.abspath(args.out)) or os.getcwd()
    if not os.path.isdir(output_dir):
        sys.exit("[ERROR] Output directory does not exist: {}".format(output_dir))

    if not 0 < args.min_af < 1:
        sys.exit("[ERROR] --min_af must lie strictly between 0 and 1, got {}".format(args.min_af))
    if args.rd_floor < 1:
        sys.exit("[ERROR] --rd_floor must be a positive integer, got {}".format(args.rd_floor))
    if args.min_mq < 0:
        sys.exit("[ERROR] --min_mq must be >= 0, got {}".format(args.min_mq))
    if args.max_id_len < 1:
        sys.exit("[ERROR] --max_id_len must be a positive integer, got {}".format(args.max_id_len))
    if args.batch_size < 1:
        sys.exit("[ERROR] --batch_size must be a positive integer, got {}".format(args.batch_size))
    if not args.sample or any(char.isspace() for char in args.sample):
        sys.exit("[ERROR] --sample must be non-empty and free of whitespace, got {!r}".format(
            args.sample))

    if args.bed is not None and args.chroms is not None:
        console.warn("--chroms is ignored because --bed was provided")

    registry = MODEL_REGISTRY[args.ploidy]
    if args.species is None:
        args.species = default_species(args.ploidy)
    elif args.species not in registry["choices"]:
        sys.exit(
            "[ERROR] --species '{species}' is not available for ploidy {ploidy}.\n"
            "  Available choices: {choices}\n"
            "  Default: {default}".format(
                species=args.species, ploidy=args.ploidy,
                choices=", ".join(species_choices(args.ploidy)),
                default=default_species(args.ploidy))
        )

    weights_path = resolve_weights(args.species, args.ploidy)
    if not os.path.isfile(weights_path):
        sys.exit(
            "[ERROR] Model weights not found: {weights}\n"
            "  Current models root: {root}\n"
            "  Set DEEPCALLER_MODELS_ROOT if the models directory lives elsewhere:\n"
            "    export DEEPCALLER_MODELS_ROOT=/path/to/models".format(
                weights=weights_path, root=MODELS_ROOT)
        )

    console.info("Ploidy: {}  |  Species model: {}".format(args.ploidy, args.species))
    console.ok("Input validation complete")
    return weights_path


def configure_hardware(requested_cpus):
    """Decide how many cores the run may use and report the choice."""
    console.banner("Hardware configuration")

    available = multiprocessing.cpu_count()
    if requested_cpus is None:
        used = available
        console.info("Auto-selected CPU threads: {}/{}".format(used, available))
    elif requested_cpus > available:
        used = available
        console.warn("Requested {} CPU threads but only {} are available; "
                     "using {}".format(requested_cpus, available, available))
    else:
        used = requested_cpus
        console.info("CPU threads: {}/{}".format(used, available))

    console.info("Execution device: CPU")
    console.ok("Hardware configuration complete")
    return used


def run_genotyper(args, chrom, chrom_depth, weights_path, num_threads, work_dir):
    """
    Genotype one chromosome in a dedicated subprocess.

    Isolation is what bounds the peak memory of a whole-genome run: whatever the
    chromosome allocated is returned to the system the moment it exits.
    """
    command = [
        sys.executable, "-m", GENOTYPER_MODULE,
        "--ref", args.ref,
        "--bam", args.bam_in_use,
        "--chrom", chrom,
        "--chrom_dp", str(chrom_depth),
        "--model_path", weights_path,
        "--ploidy", str(args.ploidy),
        "--min_af", str(args.min_af),
        "--rd_floor", str(args.rd_floor),
        "--min_mq", str(args.min_mq),
        "--max_id_len", str(args.max_id_len),
        "--batch_size", str(args.batch_size),
        "--cpus", str(num_threads),
        "--work_dir", work_dir,
    ]
    if args.bed is not None:
        command += ["--bed", args.bed]

    # Keep the package importable when DeepCaller is run straight from a source
    # tree rather than from an installed distribution.
    environment = dict(os.environ)
    package_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    existing = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        package_root if not existing else package_root + os.pathsep + existing)

    subprocess.run(command, check=True, env=environment)


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def main():
    args = build_parser().parse_args()
    started_at = time.perf_counter()

    work_dir = args.work_dir or os.path.join(
        os.getcwd(), "DeepCaller_tmp_{}".format(uuid.uuid4().hex[:8]))
    os.makedirs(work_dir, exist_ok=True)

    vcf_path = args.out if os.path.isabs(args.out) else os.path.join(os.getcwd(), args.out)

    # ---- Input validation -------------------------------------------------- #
    weights_path = validate_inputs(args)

    # ---- Hardware ---------------------------------------------------------- #
    try:
        num_threads = configure_hardware(args.cpus)
    except Exception as exc:
        sys.exit("\n[ERROR] Hardware initialisation failed:\n{}".format(exc))

    # ---- Reference --------------------------------------------------------- #
    try:
        console.banner("Reference preparation")
        console.info("Reference: {}".format(args.ref))
        check_reference(args.ref)
        console.ok("Reference preparation complete")
    except SystemExit:
        raise
    except Exception as exc:
        sys.exit("\n[ERROR] Reference preparation failed:\n{}".format(exc))

    # ---- Chromosome selection ---------------------------------------------- #
    try:
        console.banner("Chromosome selection")
        chrom_list = get_chrom_list(args.ref, args.bed, args.chroms)
        console.info("Selected {} chromosome(s)".format(len(chrom_list)))
        console.ok("Chromosome selection complete")
    except SystemExit:
        raise
    except Exception as exc:
        sys.exit("\n[ERROR] Chromosome selection failed:\n{}".format(exc))

    # ---- Alignments -------------------------------------------------------- #
    try:
        console.banner("Alignment preparation")
        console.info("Alignments: {}".format(args.bam))
        check_alignments(args.bam, num_threads)
        console.ok("Alignment preparation complete")
    except SystemExit:
        raise
    except Exception as exc:
        sys.exit("\n[ERROR] Alignment preparation failed:\n{}".format(exc))

    # ---- Coverage ---------------------------------------------------------- #
    try:
        console.banner("Coverage estimation")
        chrom_depth_map, genome_mean_depth = compute_depth(
            args.bam, work_dir, num_threads, prefix="genome")
        console.info("Genome-wide mean depth: {:.1f}X".format(genome_mean_depth))
        args.bam_in_use = args.bam
        console.ok("Coverage estimation complete")
    except Exception as exc:
        sys.exit("\n[ERROR] Coverage estimation failed:\n{}".format(exc))

    # ---- Optional downsampling --------------------------------------------- #
    if args.downsample:
        try:
            console.banner("Downsampling")
            target_depth = DOWNSAMPLE_TARGET_DEPTH[args.ploidy]
            if genome_mean_depth > target_depth:
                fraction = target_depth / genome_mean_depth
                console.info("Genome-wide depth {:.1f}X exceeds the {}X target, "
                             "sampling {:.2%} of the reads".format(
                                 genome_mean_depth, target_depth, fraction))
                args.bam_in_use, chrom_depth_map = downsample_genome(
                    args.bam, work_dir, fraction, args.seed, num_threads)
                console.info("Downsampled alignments: {}".format(args.bam_in_use))
            else:
                console.info("Genome-wide depth {:.1f}X is at or below the {}X target, "
                             "downsampling skipped".format(genome_mean_depth, target_depth))
            console.ok("Downsampling complete")
        except Exception as exc:
            sys.exit("\n[ERROR] Downsampling failed:\n{}".format(exc))

    chrom_depth_map = {c: d for c, d in chrom_depth_map.items() if c in chrom_list}
    missing = set(chrom_list) - set(chrom_depth_map)
    if missing:
        sys.exit("[ERROR] No coverage information for chromosome(s): {}".format(
            sorted(missing)))

    # ---- Encoding and genotyping ------------------------------------------- #
    console.banner("Variant encoding and genotyping")

    for index, chrom in enumerate(chrom_list, start=1):
        chrom_length, chrom_depth = chrom_depth_map[chrom]
        console.rule("{chrom} [{i}/{n}] | {size} | mean depth {depth:.1f}X".format(
            chrom=chrom, i=index, n=len(chrom_list),
            size=console.bases(chrom_length), depth=chrom_depth))
        try:
            run_genotyper(args, chrom, chrom_depth, weights_path, num_threads, work_dir)
        except subprocess.CalledProcessError as exc:
            sys.exit("\n[ERROR] Encoding and genotyping failed on {}:\n{}".format(chrom, exc))

    console.rule()
    console.ok("All chromosomes processed")

    # ---- VCF --------------------------------------------------------------- #
    try:
        console.banner("VCF generation")
        with console.Timer() as timer:
            indexed_vcf, record_count = generate_vcf(
                ref_path=args.ref,
                chrom_list=chrom_list,
                out_vcf_path=vcf_path,
                ploidy=args.ploidy,
                work_dir=work_dir,
                sample_name=args.sample,
            )
        console.info("Records written: {}".format(console.count(record_count)))
        console.info("Elapsed: {}".format(console.duration(timer.seconds)))
        console.ok("VCF generated: {}".format(indexed_vcf))
    except Exception as exc:
        sys.exit("\n[ERROR] VCF generation failed:\n{}".format(exc))

    # ---- Cleanup ----------------------------------------------------------- #
    if args.keep_tmp:
        console.info("Temp directory kept for debugging: {}".format(work_dir))
    else:
        try:
            shutil.rmtree(work_dir)
        except Exception as exc:
            console.warn("Cleanup incomplete: {}".format(exc))

    console.banner("Variant calling successful")
    console.info("Total runtime: {}".format(console.duration(time.perf_counter() - started_at)))


if __name__ == "__main__":
    main()
