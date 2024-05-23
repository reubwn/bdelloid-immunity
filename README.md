# bdelloid-immunity
_Scripts for analysis of DE experiment_

A repository for all custom analysis and plotting scripts for the manuscript 'Bdelloid rotifers deploy horizontally acquired biosynthetic genes against a fungal pathogen' (Nowell et al. 2024).

## 0. Install various software

Via conda:
```
>> conda install -c bioconda bbmap ncbi-datasets-cli sra-tools gnu-parallel
```
+ SRA Tools documentation [here](https://github.com/ncbi/sra-tools)
+ NCBI Datasets documentation [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (such a great program!)
+ BBTools documentation [here](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/).

## 1. Download raw reads and run sequence QC

### 1.1. Download raw reads from SRA

Link to raw data (SRA Run Selector permanent URL) [here](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=3&WebEnv=MCID_664cb51e8626ff46afab21f6&o=acc_s%3Aa&s=ERR4469891,ERR4469902,ERR4469903,ERR4469904,ERR4469905,ERR4469906,ERR4469907,ERR4469908,ERR4471099,ERR4471100,ERR4471101,ERR4471102,ERR4471104,ERR4471105,ERR4471106,ERR4471107,ERR4471108,ERR4471109,ERR4471110,ERR4471111,ERR4471113,ERR4471114,ERR4471115,ERR4471116#).

Or use list of accessions below:
```
>> while read acc; do fasterq-dump $acc; done < SRR_Acc_List.txt
```
[SRR_Acc_List.txt](misc/SRR_Acc_List.txt)

### 1.2. File compression (optional):
```
>> ls *fastq | parallel gzip {}
```

### 1.3. Run QC pipeline:

#### Quality and adapter trimming:
```
## run bbduk
>> bbduk.sh -Xmx60g t=$THREADS \
    in1=$READS1 in2=$READS2 \
    out1=$OUT1 out2=$OUT2 \
    ref=path/to/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo stats=bbduk.contaminants

## run fastqc
>> fastqc -d ./ --threads 6 --nogroup $READS1 $READS2 $OUT1 $OUT2
```
Note user defined parameters e.g. `$THREADS`, `$READS1` etc! Edit as required for your system.

#### Filter reads mapping to rRNA databases or _E. coli_ OP50 genome (rotifer food):

Link to SILVA rRNA database [here](https://www.arb-silva.de/).

```
## download SILVA LSU and SSU rrna databases
>> wget https://ftp.arb-silva.de/release_132/Exports/SILVA_132_LSUParc_tax_silva.fasta.gz
>> wget https://ftp.arb-silva.de/release_132/Exports/SILVA_132_LSUParc_tax_silva.fasta.gz.md5 && md5sum -c SILVA_132_LSUParc_tax_silva.fasta.gz.md5
>> wget https://ftp.arb-silva.de/release_132/Exports/SILVA_132_SSUParc_tax_silva.fasta.gz
>> wget https://ftp.arb-silva.de/release_132/Exports/SILVA_132_SSUParc_tax_silva.fasta.gz.md5 && md5sum -c SILVA_132_SSUParc_tax_silva.fasta.gz.md5

## download e. coli OP50 genome
>> datasets download genome accession GCF_009496595.1

## cat them together
>> cat ncbi_dataset/data/GCF_009496595.1/GCF_009496595.1_ASM949659v1_genomic.fna \
    <(zcat SILVA_132_LSUParc_tax_silva.fasta.gz SILVA_132_SSUParc_tax_silva.fasta.gz | perl -lane 'if(/>/){print $F[0]}else{s/U/T/g;print}') \
    | gzip > rrna_op50_db.fa.gz

## keep reads which DON'T map using outu
>> bbmap.sh -Xmx120g \
    threads=12 \
    nodisk=t \
    ref=rrna_op50_db.fa.gz \
    in1=$READS1 \
    in2=$READS2 \
    local=t \
    outu=${PREFIX}_#.filtered.fq.gz
```

## 2. Run differential expression analysis

### 2.1. Generate target files

'Gentrome' files are concatenated fasta files containing the reference transcriptome + the genomic scaffolds for each species, _A. vaga_ (Av13, GCA_000513175.1; [Flot et al. 2013](http://dx.doi.org/10.1038/nature12326)) and _A. ricciae_ (Ar18, GCA_900240375.1; [Nowell et al. 2018](http://dx.doi.org/10.1371/journal.pbio.2004830)). 'Decoy' files determine genomic scaffolds from target transcriptome sequences. For more information see the Salmon docs [here](https://salmon.readthedocs.io/en/latest/).

Transcriptome files (filtered for length >= 150 bases):
+ _A. vaga_: [cds_Av.fa.gz](data/cds_Av.fa.gz)
+ _A. ricciae_: [cds_Ar.fa.gz](data/cds_Ar.fa.gz)

Gentrome files (concatenated transcriptome + genome):
+ _A. vaga_: [gentrome_Av.fa.gz](data/gentrome_Av.fa.gz)
+ _A. ricciae_: [gentrome_Ar.fa.gz](data/gentrome_Ar.fa.gz)

Decoy files (plain text list of genome sequence IDs):
+ _A. vaga_: [decoys_Av.txt](data/decoys_Av.txt)
+ _A. ricciae_: [decoys_Ar.txt](data/decoys_Ar.txt)

#### Salmon indexing

Make Salmon index files:
```
>> parallel salmon index -t gentrome_{}.fa.gz -i cds_{}.salmon.idx -d decoys_{}.txt -p 4 ::: Ar Av
```

### 2.2. Quantification using Trinity pipeline

The Trinity RNA-seq assembler wiki has a great documentation on DE analysis, see [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki). Many of the steps below follow the Trinity tutorial.

#### Quantification using Salmon

Using the script `align_and_estimate_abundance.pl` from the [Trinity package](https://github.com/trinityrnaseq/trinityrnaseq/tree/master).

First need to generate a space-delim samples file specifying the grouping hierarchy and absolute paths to the fq read files, one for each species.

_A. vaga_ samples file:
```
AD8X24	AD8X24b	/path/to/AD8X24b/AD8X24b_1.filtered.fq.gz	/path/to/AD8X24b/AD8X24b_2.filtered.fq.gz
AD8X24	AD8X24c	/path/to/AD8X24c/AD8X24c_1.filtered.fq.gz	/path/to/AD8X24c/AD8X24c_2.filtered.fq.gz
AD8X24	AD8X24d	/path/to/AD8X24d/AD8X24d_1.filtered.fq.gz	/path/to/AD8X24d/AD8X24d_2.filtered.fq.gz
AD8X7	AD8X7a	/path/to/AD8X7a/AD8X7a_1.filtered.fq.gz	/path/to/AD8X7a/AD8X7a_2.filtered.fq.gz
AD8X7	AD8X7b	/path/to/AD8X7b/AD8X7b_1.filtered.fq.gz	/path/to/AD8X7b/AD8X7b_2.filtered.fq.gz
AD8X7	AD8X7c	/path/to/AD8X7c/AD8X7c_1.filtered.fq.gz	/path/to/AD8X7c/AD8X7c_2.filtered.fq.gz
AD8Y24	AD8Y24a	/path/to/AD8Y24a/AD8Y24a_1.filtered.fq.gz	/path/to/AD8Y24a/AD8Y24a_2.filtered.fq.gz
AD8Y24	AD8Y24b	/path/to/AD8Y24b/AD8Y24b_1.filtered.fq.gz	/path/to/AD8Y24b/AD8Y24b_2.filtered.fq.gz
AD8Y24	AD8Y24c	/path/to/AD8Y24c/AD8Y24c_1.filtered.fq.gz	/path/to/AD8Y24c/AD8Y24c_2.filtered.fq.gz
AD8Y7	AD8Y7a	/path/to/AD8Y7a/AD8Y7a_1.filtered.fq.gz	/path/to/AD8Y7a/AD8Y7a_2.filtered.fq.gz
AD8Y7	AD8Y7c	/path/to/AD8Y7c/AD8Y7c_1.filtered.fq.gz	/path/to/AD8Y7c/AD8Y7c_2.filtered.fq.gz
AD8Y7	AD8Y7d	/path/to/AD8Y7d/AD8Y7d_1.filtered.fq.gz	/path/to/AD8Y7d/AD8Y7d_2.filtered.fq.gz
```

_A. ricciae_ samples file:
```
AD1X24	AD1X24b	/path/to/AD1X24b/AD1X24b_1.filtered.fq.gz	/path/to/AD1X24b/AD1X24b_2.filtered.fq.gz
AD1X24	AD1X24c	/path/to/AD1X24c/AD1X24c_1.filtered.fq.gz	/path/to/AD1X24c/AD1X24c_2.filtered.fq.gz
AD1X24	AD1X24d	/path/to/AD1X24d/AD1X24d_1.filtered.fq.gz	/path/to/AD1X24d/AD1X24d_2.filtered.fq.gz
AD1X7	AD1X7a	/path/to/AD1X7a/AD1X7a_1.filtered.fq.gz	/path/to/AD1X7a/AD1X7a_2.filtered.fq.gz
AD1X7	AD1X7b	/path/to/AD1X7b/AD1X7b_1.filtered.fq.gz	/path/to/AD1X7b/AD1X7b_2.filtered.fq.gz
AD1X7	AD1X7c	/path/to/AD1X7c/AD1X7c_1.filtered.fq.gz	/path/to/AD1X7c/AD1X7c_2.filtered.fq.gz
AD1Y24	AD1Y24b	/path/to/AD1Y24b/AD1Y24b_1.filtered.fq.gz	/path/to/AD1Y24b/AD1Y24b_2.filtered.fq.gz
AD1Y24	AD1Y24c	/path/to/AD1Y24c/AD1Y24c_1.filtered.fq.gz	/path/to/AD1Y24c/AD1Y24c_2.filtered.fq.gz
AD1Y24	AD1Y24d	/path/to/AD1Y24d/AD1Y24d_1.filtered.fq.gz	/path/to/AD1Y24d/AD1Y24d_2.filtered.fq.gz
AD1Y7	AD1Y7a	/path/to/AD1Y7a/AD1Y7a_1.filtered.fq.gz	/path/to/AD1Y7a/AD1Y7a_2.filtered.fq.gz
AD1Y7	AD1Y7c	/path/to/AD1Y7c/AD1Y7c_1.filtered.fq.gz	/path/to/AD1Y7c/AD1Y7c_2.filtered.fq.gz
AD1Y7	AD1Y7d	/path/to/AD1Y7d/AD1Y7d_1.filtered.fq.gz	/path/to/AD1Y7d/AD1Y7d_2.filtered.fq.gz
```

Filename codes are as follows:
+ AD8 = _A. vaga_
+ AD1 = _A. ricciae_
+ X = control sample
+ Y = treatment sample
+ 7 = timepoint 7h
+ 24 = timepoint 24h
+ a,b,c,d in second column = replicate samples within each treatment group (3 per group)

### 2.3. Differential expression analyses

#### Main results - DESeq2

#### Exploration of parameter space (alternative programs and significance thresholds)

## 3. Collate HGTc information into final results files

Use custom script `collate_DE_results.pl`
