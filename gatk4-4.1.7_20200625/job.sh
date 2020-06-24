#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l mem=40G
#PBS -l walltime=60:00:00
#PBS -q medium_ext
#PBS -N cfDNA

cd $PBS_O_WORKDIR

# do something  &&\
export PATH=/home/lizhixin/softwares/anaconda3/bin:$PATH
export PERL5LIB=/home/lizhixin/softwares/anaconda3/lib/perl5/site_perl/5.22.0
export R_LIBS_USER=/home/lizhixin/softwares/R_lib_361
export R_LIBS=/home/lizhixin/softwares/R_lib_361

# AdapterRemoval --file1 161230_I136_FCHGK5HBBXX_L6_WHRDHUMouwUAHAAAB-52_1.fq.gz 161230_I136_FCHGK5HBBXX_L6_WHRDHUMouwUAHABAB-54_1.fq.gz 170101_I137_FCHGK5GBBXX_L3_WHRDHUMouwUAHAAAB-52_1.fq.gz 170101_I137_FCHGK5GBBXX_L3_WHRDHUMouwUAHABAB-54_1.fq.gz 170110_I188_FCCA1F3ANXX_L6_WHRDHUMouwUAFAAAB-52_1.fq.gz 170110_I188_FCCA1F3ANXX_L7_WHRDHUMouwUAFABAB-54_1.fq.gz --file2 161230_I136_FCHGK5HBBXX_L6_WHRDHUMouwUAHAAAB-52_2.fq.gz 161230_I136_FCHGK5HBBXX_L6_WHRDHUMouwUAHABAB-54_2.fq.gz 170101_I137_FCHGK5GBBXX_L3_WHRDHUMouwUAHAAAB-52_2.fq.gz 170101_I137_FCHGK5GBBXX_L3_WHRDHUMouwUAHABAB-54_2.fq.gz 170110_I188_FCCA1F3ANXX_L6_WHRDHUMouwUAFAAAB-52_2.fq.gz 170110_I188_FCCA1F3ANXX_L7_WHRDHUMouwUAFABAB-54_2.fq.gz --threads 10 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA --basename zstar_cfDNA --gzip &&\
# fastqc -t 10 zstar_cfDNA.pair1.truncated.gz zstar_cfDNA.pair2.truncated.gz &&\
#bwa  mem -t 10 -M ~/references/index/bwa/hg19 zstar_cfDNA.pair1.truncated.gz zstar_cfDNA.pair2.truncated.gz | samtools sort -@ 10 - -o zstar.sorted.bam &&\
#samtools index -@ 10 zstar.sorted.bam &&\


#picard MarkDuplicates \
#java -Xmx20g -XX:ParallelGCThreads=8 -jar /home/lizhixin/softwares/gatk-4.1.7.0/picard.jar MarkDuplicates \
#        I=zstar.sorted.addhead.bam \
#        O=zstar.sorted.addhead.rmdup.bam \
#        VALIDATION_STRINGENCY=LENIENT \
#        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
#        REMOVE_DUPLICATES= true \
#        M=zstar.sorted.addhead.rmdup.metric &&\

#samtools indes -@ 10 zstar.sorted.addhead.rmdup.bam &&\

#time gatk BaseRecalibrator \
#-I zstar.sorted.addhead.rmdup.bam \
#-R /home/lizhixin/references/genome/hg19/hg19.fa \
#--known-sites /home/lizhixin/annotations/gatk/1000G_phase1.indels.hg19.sites.vcf \
#--known-sites /home/lizhixin/annotations/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
#--known-sites /home/lizhixin/annotations/gatk/dbsnp_138.hg19.vcf \
#-O zstar.sorted.markdup.recal_data.table && echo "** sorted.markdup.recal_data.table done **" &&\

#time gatk ApplyBQSR \
#--bqsr-recal-file zstar.sorted.markdup.recal_data.table \
#-R /home/lizhixin/references/genome/hg19/hg19.fa \
#-I zstar.sorted.addhead.rmdup.bam \
#-O zstar.sorted.markdup.BQSR.bam && echo "** ApplyBQSR done **" &&\

#time samtools index zstar.sorted.markdup.BQSR.bam && echo "** sorted.markdup.BQSR.bam index done **" &&\

#time gatk HaplotypeCaller \
#-R /home/lizhixin/references/genome/hg19/hg19.fa \
#-I zstar.sorted.markdup.BQSR.bam \
#-O zstar.HC.vcf.gz && echo echo "** HC.vcf.gz done ** " &&\

#time gatk VariantRecalibrator \
#--reference /home/lizhixin/references/genome/hg19/hg19.fa \
#--variant zstar.HC.vcf.gz \
#--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/lizhixin/annotations/gatk/hapmap_3.3.hg19.sites.vcf \
#--resource:omini,known=false,training=true,truth=false,prior=12.0 /home/lizhixin/annotations/gatk/1000G_omni2.5.hg19.sites.vcf \
#--resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/lizhixin/annotations/gatk/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
#--resource:dbsnp,known=true,training=false,truth=false,prior=6.0 /home/lizhixin/annotations/gatk/dbsnp_138.hg19.vcf \
#--use-annotation DP --use-annotation QD --use-annotation FS --use-annotation SOR --use-annotation ReadPosRankSum --use-annotation MQRankSum \
#--mode SNP \
#--truth-sensitivity-tranche 100.0 --truth-sensitivity-tranche 99.9 --truth-sensitivity-tranche 99.0 --truth-sensitivity-tranche 95.0 --truth-sensitivity-tranche 90.0 \
#--rscript-file zstar.HC.snps.plots.R \
#--tranches-file zstar.HC.snps.tranches \
#--output zstar.HC.snps.recal && \

time gatk ApplyVQSR \
--reference /home/lizhixin/references/genome/hg19/hg19.fa \
--variant zstar.HC.vcf.gz \
--truth-sensitivity-filter-level 99.0 \
--tranches-file zstar.HC.snps.tranches \
--recal-file zstar.HC.snps.recal \
--mode SNP \
--output zstar.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **" &&\

# have not test indel calling
#time $gatk VariantRecalibrator \
#    -R $reference/Homo_sapiens_assembly38.fasta \
#    -input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
#    -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
#    -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
#    -mode INDEL \
#    --max-gaussians 6 \
#    -rscriptFile $outdir/gatk/${sample}.HC.snps.indels.plots.R \
#    --tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
#    -O $outdir/gatk/${sample}.HC.snps.indels.recal && \

#time $gatk ApplyVQSR \
#    -R $reference/Homo_sapiens_assembly38.fasta \
#    -input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
#    --ts_filter_level 99.0 \
#    --tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
#    -recalFile $outdir/gatk/${sample}.HC.snps.indels.recal \
#    -mode INDEL \
#    -O $outdir/gatk/${sample}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

echo job-done
# multiqc .
