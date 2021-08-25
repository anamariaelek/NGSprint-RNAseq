#!/bin/bash

function download_sra {

        local srr=$1
        local out=$2
        local nth=$3

        if [ -s "${out}/${srr}_2.fastq.gz" ] || [ -s "${out}/${srr}.fastq.gz" ] ; then 
                echo "# Already found: ${out}/${srr} reads $(ls ${out}/${srr}*.fastq.gz)"
        else
                # check if file already exists and size is non-zero
                fasterq-dump ${srr} -O ${out} --temp ${out} -e ${nth}
                echo "Zipping files"
                if [ -s "${out}/${srr}_1.fastq" ] ; then
                        gzip ${out}/${srr}_1.fastq
                        gzip ${out}/${srr}_2.fastq 
                else
                        gzip ${out}/${srr}.fastq 
                fi
        fi

}

ri=$1
out=$2
nth=$3
echo "# I will download $ri to $out"
mkdir -p ${out}
ris=$( echo $ri | tr "," " " )
for srr in $ris
do
 echo "# Download $rii"
 download_sra ${srr} ${out} ${nth}
done
wait
echo "DONE"
