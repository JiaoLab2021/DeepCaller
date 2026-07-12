import os
import time
import pysam
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Pool
from .utils import tetra_gt_list, hexa_gt_list


def format_vcf_record(args):
    """
    Convert a variant record dictionary to a VCF format string.

    Args:
        args: Tuple of (row_dict, ploidy) where row_dict contains:
            - chrom (str): Chromosome name
            - pos (int): Genomic position (1-based)
            - ref (str): Reference allele
            - alt (str): Alternate allele
            - pred_label (int): Genotype prediction label
            - pred_prob (float): Prediction probability [0, 1]
            - ref_num (int): Reference allele read count
            - alt_num (int): Alternate allele read count
            - rd (int): Read depth

    Returns:
        str: Formatted VCF line ending with newline
    """
    row_dict, ploidy = args

    try:
        chrom  = row_dict['chrom']
        pos    = int(row_dict['pos'])
        ref    = row_dict['ref']
        alt    = row_dict['alt']
        label  = int(row_dict['pred_label'])
        alts    = alt if label != 0 else "."
        filters = 'PASS' if label != 0 else "RefCall"

        probs    = np.clip(row_dict['pred_prob'], 1e-8, 1.0)
        qual     = np.round(-10 * np.log10(1 - probs), 2)
        gq       = int(qual)
        genotype = (tetra_gt_list if ploidy == 4 else hexa_gt_list)[label]

        ref_num = row_dict['ref_num']
        alt_num = row_dict['alt_num']
        rd      = max(row_dict['rd'], 1)
        ad      = f"{ref_num},{alt_num}"
        af      = round(alt_num / rd, 4)

        return (
            f"{chrom}\t{pos}\t.\t{ref}\t{alts}\t{qual}\t{filters}\t.\t"
            f"GT:GQ:DP:AD:AF\t{genotype}:{gq}:{rd}:{ad}:{af}\n"
        )

    except KeyError as e:
        raise KeyError(f"Missing required field in row_dict: {str(e)}")

    except (ValueError, TypeError) as e:
        raise ValueError(f"Invalid data type in row_dict: {str(e)}")


def generate_vcf(ref_path, chrom_list, num_threads, vcf_file, ploidy, work_dir):
    """
    Generate a compressed and indexed VCF file from per-chromosome Parquet data.

    Args:
        ref_path:    Path to reference genome FASTA file
        chrom_list:  List of chromosome names to process
        num_threads: Number of threads for parallel processing
        vcf_file:    Output VCF path (will be bgzip-compressed to .gz)
        ploidy:      Ploidy level (4 or 6)
        work_dir:    Directory containing per-chromosome Parquet files

    Returns:
        None: Writes output to {vcf_file}.gz and creates tabix index

    Raises:
        FileNotFoundError: If input Parquet files or reference FASTA are missing
        RuntimeError: If bgzip/tabix commands fail
    """
    try:
        print('[INFO] Generating VCF file...')

        # Get chromosome lengths from reference
        chrom_lengths = {}
        with pysam.FastaFile(ref_path) as fa:
            for chrom in chrom_list:
                chrom_lengths[chrom] = fa.get_reference_length(chrom)

        # Write VCF header and records
        with open(vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=DeepCaller\n")
	        f.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
            f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
            f.write("##FILTER=<ID=RefCall,Description=\"Genotyping model thinks this site is reference.\">\n")
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            f.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
            f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n")
            f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
            f.write("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n")

            for chrom, length in chrom_lengths.items():
                f.write(f"##contig=<ID={chrom},length={length}>\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

            for chrom in chrom_list:
                parquet_path = os.path.join(work_dir, f"{chrom}.parquet")

                if not os.path.exists(parquet_path):
                    raise FileNotFoundError(f"Missing Parquet file for {chrom}: {parquet_path}")

                chrom_df     = pd.read_parquet(parquet_path)
                records_args = [(row._asdict(), ploidy) for row in chrom_df.itertuples(index=False)]

                with Pool(processes=num_threads) as pool:
                    records = pool.map(format_vcf_record, records_args)

                f.writelines(records)

        # Compress and index
        subprocess.run(f"bgzip -f {vcf_file}", shell=True, check=True)
        subprocess.run(f"tabix -p vcf {vcf_file}.gz", shell=True, check=True)

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"bgzip/tabix failed: {e.stderr.strip()}")

    except pysam.SamtoolsError as e:
        raise RuntimeError(f"Reference FASTA error: {str(e)}")

    except Exception as e:
        raise RuntimeError(f"VCF generation failed: {str(e)}")
