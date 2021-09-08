# NGSprint2021 RNA-seq hackathon

When: 8 - 15.09.2021

Where: Online

More info about the event can be found on the [website](https://ngschool.eu/ngsprint). 

# Requirements

For this hackathon, it is recommended to work on a linux operating system with minimum 6 cores available to run the basic analyses.  

Below is a suggestion for how to install some of the tools you might find useful during the hackathon, however, you may use any other tool to solve the problems as well.

You can create a conda enviroment with dependencies and potentially useful tools:

```bash
conda create -y -n rna python=3.6
conda activate rna
conda install -y -c conda-forge r=4.1.0
conda install -y -c bioconda sra-tools parallel-fastq-dump fastqc multiqc trimmomatic trim-galore star hisat2 subread samtools bedtools
conda env export --name rna > rna_env.yml
```

After activating the enviroment, install R packages directly from R (not using conda here because of the version conflicts).

```r
install.packages("BiocManager")
BiocManager::install(pkgs = c("Rsubread","DESeq2","fgsea","clusterprofiler","topGO"))
```

# Resources

* [paper](https://doi.org/10.1242/dmm.030536)
* [raw data](https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE104288)
* [RNA-seq analysis review](https://doi.org/10.1186/s13059-016-0881-8)
