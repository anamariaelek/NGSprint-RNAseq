# Requirements

Create conda enviroment with all dependencies (you can choose toinstall only some of these).

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

1. QC

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```bash
nth=6
nname="SRR6078292"
fastq_r1="raw_reads/${name}_1.fastq.gz"
fastq_r2="raw_reads/${name}_2.fastq.gz"
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
	ILLUMINACLIP:trimmed_reads/adapters_NexteraPE-PE.fasta:2:30:10 |& tee -a ${log}
	
# re-check with fastqc
fastqc -t ${nth} trimmed_reads/${name}_1.tp.fastq.gz trimmed_reads/${name}_2.tp.fastq.gz |& tee -a ${log}
```

Trimmomatic can take up 15-20 minutes for a pair of files.

2. Mapping

Get reference genome and annotation from [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39).

[STAR](https://github.com/alexdobin/STAR)

```
mkdir -p reference/hg38_index

# genome index
STAR --runThreadN 6 \
	--runMode genomeGenerate \
	--genomeDir reference/hg38_index \
	--genomeFastaFiles reference/Homo_sapiens.GRCh38.dna.fa \
	--sjdbGTFfile reference/Homo_sapiens.GRCh38.92.gtf \
	--sjdbOverhang 99

# mapping

```

[HISAT2](http://daehwankimlab.github.io/hisat2/)
[Subread](http://subread.sourceforge.net/)


# Other resources

* [mapping wiith STAR on HPC](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)
* [mapping with HISAT2](https://wikis.utexas.edu/display/bioiteam/Mapping+with+HISAT2)
* [Rsubread paper](https://academic.oup.com/nar/article/47/8/e47/5345150)
* [benchmarking RNA-seq algners paper](https://www.nature.com/articles/nmeth.4106)

