#!/bin/bash

echo "Trinity assembly in progress"

Trinity --seqType fq --max_memory 100G -–left all_reads_R1.fa.fq –-right all_reads_R2.fa.fq --SS_lib_type RF --CPU 48 --output trinity_results

echo "analysis completed"