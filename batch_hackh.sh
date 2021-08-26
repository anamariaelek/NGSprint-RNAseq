#!/bin/bash

ri=$1
nth=$2 

function hackh {

        local name=$1
        local nth=$2
	mkdir logs
	mkdir trimmed_reads

	# download fastq
        echo "# Downloading $srr"
	parallel-fastq-dump --sra-id ${name} \
		--threads ${nth} \
		--outdir raw_reads \
		--split-files --gzip |& tee -a logs/fastq-dump.${name}.tp.log

	# fastqc
	fastq_r1=raw_reads/${name}_1.fastq.gz
	fastq_r2=raw_reads/${name}_2.fastq.gz
	fastqc -t ${nth} ${fastq_r1} ${fastq_r2} |& tee -a logs/fastqc.${name}.log
	
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
		ILLUMINACLIP:trimmed_reads/adapters_NexteraPE-PE.fasta:2:30:10 |& tee -a logs/trimmomatic.${name}.log
		
	# re-check with fastqc
	fastqc -t ${nth} trimmed_reads/${name}_1.tp.fastq.gz trimmed_reads/${name}_2.tp.fastq.gz |& tee -a logs/fastqc.${name}.tp.log
	
	# map
	subread-align -t 0 \
	    -T ${nth} \
	    -i reference/GRCh38.fa.Subread_index \
	    -r trimmed_reads/${name}_1.tp.fastq.gz \
	    -R trimmed_reads/${name}_2.tp.fastq.gz \
	    -o mapping/Subread/${name}.bam |& tee -a logs/Subread.${name}.log
	

}



#ri="SRR6078288,SRR6078289,SRR6078290,SRR6078291,SRR6078292,SRR6078293"
echo "# Will download $ri to $out"
mkdir -p ${out}
ris=$( echo $ri | tr "," " " )
for rii in $ris
do
	hackh ${rii} ${nth}
done

