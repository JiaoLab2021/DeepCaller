#!/usr/bin/env python3
import os
import sys
import uuid

os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

from setproctitle import setproctitle
setproctitle("DeepCaller " + " ".join(sys.argv[1:]))

import warnings
warnings.filterwarnings('ignore', category=UserWarning)

import pysam
import shutil
import argparse
import multiprocessing
import tensorflow as tf

from multiprocessing import Pool
from .process_chrom import main as process_chrom_main
from .utils_bam import (
    check_reference,
    get_chrom_list,
    sort_bam,
    filter_bam,
    summarize_genome_depth,
)
from .utils_vcf import generate_vcf
from DeepCaller import __version__


SPECIES_MAP = {
    'potato':       ('01_C88_Potato_default',            [4]),
    'alfalfa':      ('02_Bolivia_Alfalfa',               [4]),
    'rose':         ('03_Samantha_Rose',                 [4]),
    'sweetpotato':  ('04_Tanzania_Sweetpotato_default',  [6]),
    'syn_potato':   ('05_SyntheticPotato_Potato',        [6]),
}

PLOIDY_DEFAULT_SPECIES = {
    4: 'potato',
    6: 'sweetpotato',
}


def check_hardware(requested_cpus=None):
    available_cpus = multiprocessing.cpu_count()

    if requested_cpus is None:
        used_cpus = min(24, available_cpus)
        print(f"* Auto-selected CPU threads: {used_cpus}/{available_cpus}")
    elif requested_cpus > available_cpus:
        used_cpus = available_cpus
        print(
            f"[WARNING] Requested CPUs ({requested_cpus}) exceed available "
            f"({available_cpus}), using {available_cpus}"
        )
    else:
        used_cpus = requested_cpus
        print(f"* Use CPU threads: {used_cpus}/{available_cpus}")

    tf.config.set_visible_devices([], 'GPU')
    tf.config.threading.set_intra_op_parallelism_threads(used_cpus)
    tf.config.threading.set_inter_op_parallelism_threads(used_cpus)

    print(f"* Running on CPU only")

    return used_cpus


def cpu_count_type(value):
    v = int(value)
    if v == -1:
        return None
    if v <= 0:
        raise argparse.ArgumentTypeError(f"must be -1 or a positive integer, got {v}")
    return v


def main():

    parser = argparse.ArgumentParser(
        description="Tool for identifying polyploidy small variation sites",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"DeepCaller {__version__}",
    )

    required = parser.add_argument_group('Required arguments')
    required.add_argument("-r", "--ref",    required=True, help="Reference FASTA file")
    required.add_argument("-b", "--bam",    required=True, help="Input BAM file")
    required.add_argument("-p", "--ploidy", required=True, type=int, choices=[4, 6], help="Ploidy level")

    io_group = parser.add_argument_group('Input/output configuration')
    io_group.add_argument("-c", "--chroms", nargs="+",        help="Chromosomes to include")
    io_group.add_argument("-o", "--output", default="output.vcf", help="Output VCF file")
    io_group.add_argument("-l", "--bed",    default=None,     help="Optional BED file. If provided, --chroms will be ignored")

    processing = parser.add_argument_group('Processing options')
    processing.add_argument(
        "-t", "--cpus",
        type=cpu_count_type,
        default=24,
        help="Number of CPU cores. Use -1 to use all available",
    )
    processing.add_argument(
        "-s", "--species",
        default=None,
        choices=list(SPECIES_MAP.keys()),
        help=(
            "Species model to use. "
            "Tetraploid (ploidy=4): potato (default), alfalfa, rose. "
            "Hexaploid (ploidy=6): sweetpotato (default), syn_potato. "
            "If not specified, defaults to the primary model for the given ploidy."
        ),
    )
    processing.add_argument(
        "-m", "--mode",
        default="speed",
        choices=["speed", "performance"],
        help="Inference mode: speed (faster) or performance (more accurate)",
    )
    processing.add_argument("--min_af",   type=float, default=0.10, help="Threshold of allele frequency at candidate variants")
    processing.add_argument("--rd_floor", type=int,   default=10,   help="Threshold of read depth at candidate variants")

    args = parser.parse_args()

    if args.species is None:
        args.species = PLOIDY_DEFAULT_SPECIES[args.ploidy]

    work_dir = os.path.join(os.getcwd(), f"DeepCaller_tem_dir_{uuid.uuid4().hex[:8]}")
    os.makedirs(work_dir, exist_ok=True)

    vcf_file = os.path.join(os.getcwd(), args.output)

    # ============================================================
    # STEP 0: INPUT VALIDATION
    # ============================================================
    print("=" * 80)
    print("INPUT VALIDATION".center(80))
    print("=" * 80)

    if not os.path.isfile(args.ref):
        sys.exit(f"[ERROR] Reference FASTA not found: {args.ref}")

    if not os.path.isfile(args.bam):
        sys.exit(f"[ERROR] BAM file not found: {args.bam}")

    if args.bed is not None and not os.path.isfile(args.bed):
        sys.exit(f"[ERROR] BED file not found: {args.bed}")

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if not os.path.isdir(output_dir):
        sys.exit(f"[ERROR] Output directory does not exist: {output_dir}")

    if not (0 < args.min_af < 1):
        sys.exit(f"[ERROR] --min_af must be between 0 and 1 (exclusive), got {args.min_af}")

    if args.rd_floor < 1:
        sys.exit(f"[ERROR] --rd_floor must be a positive integer, got {args.rd_floor}")

    if args.bed is not None and args.chroms is not None:
        print("[WARNING] --chroms will be ignored because --bed is provided")

    _, supported_ploidies = SPECIES_MAP[args.species]
    if args.ploidy not in supported_ploidies:
        tetra = [k for k, v in SPECIES_MAP.items() if 4 in v[1]]
        hexa  = [k for k, v in SPECIES_MAP.items() if 6 in v[1]]
        sys.exit(
            f"[ERROR] Species '{args.species}' does not support ploidy {args.ploidy}.\n"
            f"  Tetraploid (ploidy=4): {', '.join(tetra)}\n"
            f"  Hexaploid  (ploidy=6): {', '.join(hexa)}"
        )

    print(f"* Species: {args.species}  |  Mode: {args.mode}  |  Ploidy: {args.ploidy}")
    print("[OK] Input validation complete!")

    # ============================================================
    # STEP 1: HARDWARE CONFIGURATION
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("HARDWARE CONFIGURATION".center(80))
        print("=" * 80)

        num_threads = check_hardware(requested_cpus=args.cpus)

        print("[OK] Hardware configuration complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] Hardware initialization failed:\n{str(e)}")

    # ============================================================
    # STEP 2: REFERENCE VALIDATION
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("REFERENCE VALIDATION".center(80))
        print("=" * 80)

        print(f"* Reference input: {args.ref}")
        check_reference(args.ref)

        print("[OK] Reference validation complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] Reference validation failed:\n{str(e)}")

    # ============================================================
    # STEP 3: CHROMOSOME SELECTION
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("CHROMOSOME SELECTION".center(80))
        print("=" * 80)

        chrom_list = get_chrom_list(args.ref, args.bed, args.chroms)

        print(f"* Selected: {len(chrom_list)} chromosome(s)")
        print("[OK] Chromosome selection complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] Chromosome selection failed:\n{str(e)}")

    # ============================================================
    # STEP 4: BAM SORTING
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("BAM SORTING".center(80))
        print("=" * 80)

        print(f"* BAM file input: {args.bam}")
        sorted_bam = sort_bam(args.bam, num_threads, work_dir)

        print("[OK] BAM sorting complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] BAM processing failed:\n{str(e)}")

    # ============================================================
    # STEP 5: BAM FILTERING
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("BAM FILTERING".center(80))
        print("=" * 80)

        with pysam.FastaFile(args.ref) as fasta:
            all_chroms = fasta.references

        print("[INFO] Filtering chromosomes")

        task_args = [(sorted_bam, chrom, work_dir) for chrom in all_chroms]

        with Pool(processes=num_threads) as pool:
            pool.map(filter_bam, task_args)

        print("[OK] BAM filtering complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] BAM filtering failed:\n{str(e)}")

    # ============================================================
    # STEP 6: SEQUENCING DEPTH CALCULATION
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("SEQUENCING DEPTH CALCULATION".center(80))
        print("=" * 80)

        bam_depth = summarize_genome_depth(work_dir, all_chroms)
        print(f"* Use {bam_depth}X as the sequencing depth")
        print("[OK] Sequencing depth calculation complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] Depth calculation failed:\n{str(e)}")

    # ============================================================
    # STEP 7: ENCODING AND PREDICTION
    # ============================================================
    print("\n" + "=" * 80)
    print("ENCODING AND PREDICTION".center(80))
    print("=" * 80)

    for chrom in chrom_list:

        process_args = [
            "--ref",       args.ref,
            "--bam",       os.path.join(work_dir, f"temp_{chrom}.bam"),
            "--chrom",     chrom,
            "--ploidy",    str(args.ploidy),
            "--cpus",      str(num_threads),
            "--bam_depth", str(bam_depth),
            "--species",   args.species,
            "--mode",      args.mode,
            "--min_af",    str(args.min_af),
            "--rd_floor",  str(args.rd_floor),
            "--work_dir",  work_dir,
        ]

        if args.bed is not None:
            process_args.extend(["--bed", args.bed])

        original_argv = sys.argv
        try:
            sys.argv = ["process_chrom"] + process_args
            process_chrom_main()
        except Exception as e:
            sys.exit(f"\n[ERROR] Encoding and prediction failed on {chrom}:\n{str(e)}")
        finally:
            sys.argv = original_argv

    print("[OK] Chromosomes processing complete!")

    # ============================================================
    # STEP 8: VCF FILE GENERATION
    # ============================================================
    try:
        print("\n" + "=" * 80)
        print("VCF FILE GENERATION".center(80))
        print("=" * 80)

        generate_vcf(
            args.ref,
            chrom_list,
            num_threads,
            vcf_file,
            args.ploidy,
            work_dir,
        )

        print(f"[OK] VCF generated: {vcf_file}.gz")
        print("[OK] VCF file generation complete!")

    except Exception as e:
        sys.exit(f"\n[ERROR] VCF file generation failed:\n{str(e)}")

    # ============================================================
    # FINISH
    # ============================================================
    print("\n" + "=" * 80)

    try:
        shutil.rmtree(work_dir)
    except Exception as e:
        print(f"[WARNING] Cleanup incomplete: {str(e)}")

    print("VARIANTS CALL SUCCESSFUL".center(80))
    print("=" * 80)


if __name__ == "__main__":
    main()
