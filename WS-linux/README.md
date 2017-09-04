# Introduction


## paste <(zcat file1.gz) <(zcat file2.gz)

## ln -s

## awk

## xargs

## sed

## make

Create a `Makefile` with the following content:
```
STAR_FLAGS := --genomeDir ref/GRCm38_85_index
STAR_FLAGS += --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1
#STAR_FLAGS += --clip3pAdapterSeq CTGTCTCTTATACACATCT

.PRECIOUS:%_umi.fastq.gz
%_umi.fastq.gz:%_R1.fastq %_R2.fastq; paste -d'\t' $(word 1,$^) $(word 2,$^) | awk 'NR%4==1{if($$1!=$$3){print "ERROR: header mismatch at line" NR,$$0 > "/dev/stderr";exit 1} else hdr = $$3;next}NR%4==2{bc = substr($$1,1,6);umi = substr($$1,1,5);print hdr "_" umi}{print $$2}' | gzip > $@

.PRECIOUS:%.star/Aligned.sortedByCoord.out.bam %.star.bam %.star.bam.bai
%.star/Aligned.sortedByCoord.out.bam:SHELL:=bsub -K -n1 -R "span[ptile=1]" -M 32000000 -R "rusage[mem=32000]" -e "./tmp/log/stderr_%J" -o "./tmp/log/stdout_%J" /bin/bash
%.star/Aligned.sortedByCoord.out.bam:%.fastq.gz; module add UHTS/Aligner/STAR/2.5.3a;mkdir -p $(@D);STAR $(STAR_FLAGS) --runThreadN 1 --readFilesCommand gunzip -c --readFilesIn $^ --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $(@D)/
%.star.bam:%.star/Aligned.sortedByCoord.out.bam; ln -s $(notdir $*).star/Aligned.sortedByCoord.out.bam $@
%.star.bam.bai:%.star.bam; module add UHTS/Analysis/samtools/0.1.19;samtools index $^
%.star.umitools_dedup.bam:%.star.bam %.star.bam.bai; umi_tools dedup -I $< -S $@
```


