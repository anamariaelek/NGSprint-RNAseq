# Requirements

Create conda enviroment with all dependencies (you can choose to install only some of these).

```bash
conda create -n rna
conda activate rna
conda install -c bioconda sra-tools fastqc trimmomatic star hisat2 subread bioconductor-rsubread samtools bedtools
conda env export --name rna > rna_env.yml
```

THe entire enviroment comes up to 380 Mb.

Recommended is to have at least 6 cores available to run the steps described below.
How to find out how many cores do you have? Run `lscpu` command and look for these lines in the output:

```
CPU(s):                          8
Thread(s) per core:              2
Core(s) per socket:              4
Socket(s):                       1
```

This system has 1 physical CPU (socket) which has 4 cores (cores per socket) and each core has 2 threads - this all sums up to 1 x 4 x 2 = 8 CPUs. You can get this total number with `nproc`, too.  
Never use all of CPUs for the analysis you are running - your system needs some to run, too!

# Analysis

## 1. QC

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```bash
nth=6
name=SRR6078292
fastq_r1=raw_reads/${name}_1.fastq.gz
fastq_r2=raw_reads/${name}_2.fastq.gz
mkdir logs
fastqc -t ${nth} ${fastq_r1} ${fastq_r2} |& tee -a logs/fastqc.log
```

This should take a few minutes per file.

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```bash
mkdir trimmed_reads

# create mock FASTA with common adapters
echo -e ">PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" > trimmed_reads/adapters_NexteraPE-PE.fasta

# trim reads
trimmomatic PE \
	-threads ${nth} \
	-phred33 \
 	${fastq_r1} \
	${fastq_r2} \
	trimmed_reads/${name}_1.tp.fastq.gz \
	trimmed_reads/${name}_1.tu.fastq.gz \
	trimmed_reads/${name}_2.tp.fastq.gz \
	trimmed_reads/${name}_2.tu.fastq.gz \
	SLIDINGWINDOW:4:20 \
	ILLUMINACLIP:trimmed_reads/adapters_NexteraPE-PE.fasta:2:30:10 |& tee -a logs/trimmomatic.log
	
# re-check with fastqc
fastqc -t ${nth} trimmed_reads/${name}_1.tp.fastq.gz trimmed_reads/${name}_2.tp.fastq.gz |& tee -a ${log}
```

Trimmomatic can take up 15-20 minutes for a pair of files.

## 2. Mapping

Get reference genome and annotation from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/). You can download the files prepared for use in pipelines and analyses, those include some of the pregenerated genome indices.

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/*no_alt_analysis_set.fna.gz -O reference/GRCh38.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/*no_alt_analysis_set.fna.fai -O reference/GRCh38.fna.fai
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/*no_alt_analysis_set.fna.hisat2_index.tar.gz -O reference/GRCh38.fna.hisat2_index.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/*gff* -O reference/GRCh38.refseq_annotation.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/*gtf* -O reference/GRCh38.refseq_annotation.gtf.gz

```

### [STAR](https://github.com/alexdobin/STAR)

```
# genome index
mkdir -p reference/STAR_index
STAR --runThreadN 3 \
	--runMode genomeGenerate \
	--genomeDir reference/STAR_index \
	--genomeFastaFiles reference/GRCh38.fna \
	--sjdbGTFfile reference/GRCh38.refseq_annotation.gtf \
	--sjdbOverhang 99 \\
	|& tee -a logs/genomeindex.HISAT2.log

# mapping
mkdir -p mapping/STAR
STAR --genomeDir reference/STAR_index \
	--runThreadN ${nth} \
	--readFilesIn ${fastq_r1} ${fastq_r2} \
	--outFileNamePrefix mapping/STAR/${name} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard |& tee -a logs/mapping.STAR.log

```

For human genome you need you need approximately 30Gb of RAM (check how many megabytes of memory you have with `free -mega`)
Tweaking of memory usage parameters may also be required: `--limitGenomeGenerateRAM` determines total RAM available for index building and `--genomeChrBinNbits` to reduce RAM consumption per scaffold (to accomodate processing large number of scaffolds).  

### [HISAT2](http://daehwankimlab.github.io/hisat2/)

HISAT2 needs 8Gb RAM for mapping to human genome, making it potentially usable on a personal laptop. With 8Gb and 6 cores, it takes about an hour to build index, mapping about 30-40 min.  

```bash
mkdir -p mapping/HISAT2

hisat2-build reference/GRCh38.fa reference/GRCh38.fa |& tee -a logs/genomeindex.HISAT2.log

hisat2 -x reference/GRCh38.fa \
	-1 ${fastq_r1} \
	-2 ${fastq_r2} \
	-S mapping/HISAT2/${name}.sam \
	--phred33 \
	--rna-strandness RF \
	-t -p $nth |& tee -a logs/mapping.HISAT.log

samtools view -S -b -@ ${nth} mapping/HISAT2/${name}.sam > mapping/HISAT2/${name}.bam |& tee -a logs/mapping.HISAT.log
```

It is possible to run HISAT2 on SRA samples "directly", by passing comma-separated list of accession numbers as `--sra-acc` argument.

### [Subread](http://subread.sourceforge.net/)

Building full index needs 15Gb of memory for human/mouse genome, and with this index, alignment requires 17.8Gb of memory for read mapping. On personal computer, it is recommended to build gapped index instead, which needs 5.7GB of memory for human/mouse (and takes 15 minutes, as compared with 40 for building full index), and 8.2GB of memory for mapping.  
Memory use for index building and read mapping can be further reduced by building a split index using the `-B` and `-M` options.

```bash
mkdir mapping/Subread

subread-buildindex -o GRCh38.fa.Subread_index -M 8000 reference/GRCh38.fa |& tee -a logs/genomeindex.Subread.log

subread-align -t 0 \
    -T ${nth} \
    -i reference/GRCh38.fa.Subread_index \
    -r ${fastq_r1} \
    -R ${fastq_r2} \
    -o mapping/Subread/${name}.bam |& tee -a logs/Subread.log

```

### [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)

R implementaton of `Subread`

```r
library(Rsubread)
dir.create("mapping/Rsubread")

ref <- "reference/GRCh38.fa"	
name <- "SRR6078292"
fastq_r1 <- sprintf("raw_reads/%s_1.fastq.gz",name)
fastq_r2 <- sprintf("raw_reads/%s_2.fastq.gz",name)

buildindex(basename="GRCh38.fa.Rsubread_index",reference=ref)

align.stat <- align(
  index="GRCh38.fa.Rsubread_index", 
  readfile1=fastq_r1, 
  readfile2=fastq_r2,
  output_file=sprintf("mapping/Rsubread/%s.bam",name),
  phredOffset=64
)

```

## 3. Read counting 

`featureCounts`

```bash
bam=mapping/Subread/${name}.bam

featureCounts -a reference/GRCh38.refseq_annotation.gff \
    -o <output_file> \
    -T ${nth} \
    $bam |& tee -a logs/featureCounts.log
```

If using `Rsubread`, then `Rsubread::featureCounts()` can be run from R. You can use externally downloaded annotation, but for mm9, mm10, hg19 and hg38 genomes there is an in-built NCBI RefSeq gene annotation available, too.

```r
require(Rsubread)

name <- "SRR6078292"
bam_files <- list.files("mapping/Rsubread", pattern="*bam", full.names=TRUE)
fcounts <- featureCounts(bam_files, annot.inbuilt = "hg38", isPairedEnd=TRUE)

```

# Other resources

* [mapping wiith STAR on HPC](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)
* [mapping with HISAT2](https://wikis.utexas.edu/display/bioiteam/Mapping+with+HISAT2)
* [Rsubread paper](https://academic.oup.com/nar/article/47/8/e47/5345150)
* [benchmarking RNA-seq algners paper](https://www.nature.com/articles/nmeth.4106)

