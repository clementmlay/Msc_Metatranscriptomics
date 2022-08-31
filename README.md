# Msc_Metatranscriptomics
Metatranscriptomics data analysis


 This repository describes step by step the analysis of the Metatranscriptome dataset 

  A continuous effort to create a reliable pipeline for;


   1. Quality control


   2. Removing host RNAs


   3. Taxonomic assignment


   4. Denovo assembly


   5. Differential gene expression


   6. Functional annotation


   7. Pathway analysis


# codes to be automated for the pipeline

Quality control

`fastp -q 5 -l 25 -3 -M 5 -i ${sampleid}_R1_001.fastq.gz -I ${sampleid}_R2_001.fastq.gz -o ./clean_data/${sampleid}_R1_001.clean_fastq.gz -O ./clean_data/${sampleid}_R2_001.clean_fastq.gz --html ./clean_data/fastqc_results/${sampleid}_fastq_trim_report.html`

Removing host RNAs

`bbsplit.sh in1=${sampleid}_R1_001.clean_fastq.gz in2=${sampleid}_R2_001.clean_fastq.gz ref=../refferences/coral.fna.gz basename=./contaminats_free/out_%.fq.gz outu1=./contaminats_free/${sampleid}_R1_001.contfree_fastq.gz outu2=./contaminats_free/${sampleid}_R2_001.contfree_fastq.gz`

Taxonomic assignment

`kraken2 --use-names --threads 20 --db ../databases/minikraken_8GB_20200312 --fastq-input --report ../kraken_results/mysample.kreport2 --gzip-compressed --paired ${sampleid}_R1_001.contfree_fastq.gz ${sampleid}_R2_001.contfree_fastq.gz > ../kraken_results/${sampleid}.krk`

Taxonomy visualization

`rcf -n /home/cley/recentrifuge/taxdump -k ${sampleid} -c 1 -o ../recentrifuge_results/${sampleid}.html`

Denovo assembly

`Trinity --seqType fq --max_memory 100G -–left all_reads_R1.fa.fq –-right all_reads_R2.fa.fq --SS_lib_type RF --CPU 48 --output trinity_results
`
