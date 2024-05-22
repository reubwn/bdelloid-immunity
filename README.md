# bdelloid-immunity
_Scripts for analysis of DE experiment_

A repository for all custom analysis and plotting scripts for the manuscript 'Bdelloid rotifers deploy horizontally acquired biosynthetic genes against a fungal pathogen' (Nowell _et al._ 2024).

## 1. Download raw reads and run sequence QC

### 1.1. Download raw reads from SRA 

Link to raw data (SRA Run Selector permanent URL) [here](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=3&WebEnv=MCID_664cb51e8626ff46afab21f6&o=acc_s%3Aa&s=ERR4469891,ERR4469902,ERR4469903,ERR4469904,ERR4469905,ERR4469906,ERR4469907,ERR4469908,ERR4471099,ERR4471100,ERR4471101,ERR4471102,ERR4471104,ERR4471105,ERR4471106,ERR4471107,ERR4471108,ERR4471109,ERR4471110,ERR4471111,ERR4471113,ERR4471114,ERR4471115,ERR4471116#).

```
>> while read acc; do fasterq-dump $acc; done < SRR_Acc_List.txt
```
[SRR_Acc_List.txt](SRR_Acc_List.txt)

### 1.2. File compression (optional):
```
>> ls *fastq | parallel gzip {}
```

### 1.3. Run QC pipeline:

#### Install BBTools via conda:
```
>> conda install -c bioconda bbmap
```
BBTools documentation [here](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/).

#### Run basic QC pipeline for each read pair:
```
## run bbduk
>> bbduk.sh -Xmx60g t=$THREADS \
    in1=$READS1 in2=$READS2 \
    out1=$OUT1 out2=$OUT2 \
    ref=path/to/resources/adapters.fa \
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
    outu=${PREFIX}_#.bbduk.ecc.filtered.fq.gz
```

## 2. Run differential expression analysis

**2.1.** Generate target 'gentrome.fa'

Following Salmon documentation available here: 

For _A. vaga_:
```
commands
```

For _A. ricciae_
```
commands
```

**2.2.** Run Salmon quantification

```
commands
```

## 3. Collate HGTc information into Salmon results files

Use custom script `collate_DE_results.pl` 
