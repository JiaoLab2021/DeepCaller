[README.md](https://github.com/user-attachments/files/27900478/README.md)
# Demo

This directory contains a small dataset for verifying your DeepCaller installation.

## Dataset

| File | Description |
|------|-------------|
| `DM8.1_chr10_100000_1100000.fa` | Reference FASTA — chromosome 10, positions 100,000–1,100,000 |
| `DM8.1_chr10_100000_1100000.fa.fai` | FASTA index |
| `C88_20x_chr10_100000_1100000.bam` | Tetraploid potato (C88) short-read alignments at ~20× coverage |
| `C88_20x_chr10_100000_1100000.bam.bai` | BAM index |

## Run

From the `Demo/` directory:

```bash
DeepCaller \
    -r DM8.1_chr10_100000_1100000.fa \
    -b C88_20x_chr10_100000_1100000.bam \
    -p 4 \
    --mode speed \
    -o demo_output.vcf
```

`--species` is omitted here — DeepCaller defaults to `potato` when `--ploidy 4` is set.

## Expected output

- `demo_output.vcf.gz` — bgzip-compressed VCF
- `demo_output.vcf.gz.tbi` — tabix index
