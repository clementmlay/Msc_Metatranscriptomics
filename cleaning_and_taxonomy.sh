!#/bin/bash

##creating a conda environment
conda create -n metatranscriptomics

## Activate working environment
conda activate metatranscriptomics

######################### FOR GOOGLE CLOUD USERS ###########
## Install gcsfuse tool to access data on google bucket
export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
sudo apt-get update
sudo apt-get install gcsfuse

## mount google cloud bucket to acces data
gcsfuse --implicit-dirs coral_metatranscriptomics_data rawdata/
############################       END GOOGLE        #############

## filtering low quality reads and trimming

samplenames={ls | cut -d _ -f1,2,3}

for sampleid in $samplenames

do

printf "\nPIPELINE ALERT: Trimming %s\n" "${sampleid}"

fastp -q 5 -l 25 -3 -M 5 -i ${sampleid}_R1_001.fastq.gz -I ${sampleid}_R2_001.fastq.gz -o ./clean_data/${sampleid}_R1_001.clean_fastq.gz -O ./clean_data/${sampleid}_R2_001.clean_fastq.gz --html ./clean_data/fastqc_results/${sampleid}_fastq_trim_report.html

done


## determine the quality of trimmed data
fastqc ./clean_data/*.trimmed_fastq.gz -o ../fastqc_results/after_trimming

## Removing host contaminants

#samplenames={ls ../results/fastp_results | cut -d _ -f1,2,3}

for sampleid in $samplenames

do

bbsplit.sh in1=${sampleid}_R1_001.clean_fastq.gz in2=${sampleid}_R2_001.clean_fastq.gz ref=../refferences/coral.fna.gz basename=./contaminats_free/out_%.fq.gz outu1=./contaminats_free/${sampleid}_R1_001.contfree_fastq.gz outu2=./contaminats_free/${sampleid}_R2_001.contfree_fastq.gz

done

## taxinomoic classification

#samplenames={ls ./contaminants_free | cut -d _ -f1,2,3}

for sampleid in $samplenames

do

kraken2 --use-names --threads 20 --db ../databases/minikraken_8GB_20200312 --fastq-input --report ../kraken_results/mysample.kreport2 --gzip-compressed --paired ${sampleid}_R1_001.contfree_fastq.gz ${sampleid}_R2_001.contfree_fastq.gz > ../kraken_results/${sampleid}.krk

done

#Visualizing kraken results

samplenames={ls ../kraken_results}

for sampleid in $samplenames

do

rcf -n /home/cley/recentrifuge/taxdump -k ${sampleid} -c 1 -o ../recentrifuge_results/${sampleid}.html

done

echo "ANALYSIS COMPLETEDS"
