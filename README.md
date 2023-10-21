# GenomonFisher
Genomon fisher exact test mutation caller, modified for use in ogawa-lab

## Dependency
Python (>= 2.7), pysam, scipy, builtins,
samtools (>= 1.10)

## Run
Disease sample vs. Control sample Comparison
```
usage: fisher_0.2.1 comparison [-h] -1 BAM1 -2 BAM2 -o OUTPUT -r REF_FA -s SAMTOOLS_PATH
                         [-S SAMTOOLS_PARAMS] [-Q BASE_QUALITY]
                         [-m MIN_ALLELE_FREQ] [-M MAX_ALLELE_FREQ]
                         [-f FISHER_VALUE] [-d MIN_DEPTH]
                         [-v MIN_VARIANT_READ] [-R REGION]
                         [-e] [-g LOG_FILE] [-l LOG_LEVEL]

```
Single sample mutation calling
```
usage: fisher_0.2.1 single [-h] -1 BAM1 -o OUTPUT -r REF_FA -s SAMTOOLS_PATH
                     [-S SAMTOOLS_PARAMS] [-Q BASE_QUALITY]
                     [-m MIN_ALLELE_FREQ] [-p POST_10_Q] [-d MIN_DEPTH]
                     [-v MIN_VARIANT_READ] [-R REGION]
                     [-e] [-g LOG_FILE] [-l LOG_LEVEL]
```

