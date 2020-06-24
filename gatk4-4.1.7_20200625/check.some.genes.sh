cat /home/lizhixin/references/genome/hg19/key_gene.list | while read id;

do

# chr=$(echo $id |cut -d" " -f 1|sed 's/chr//' )
chr=$(echo $id |cut -d" " -f 1 )

start=$(echo $id |cut -d" " -f 2 )

end=$(echo $id |cut -d" " -f 3 )

gene=$(echo $id |cut -d" " -f 4 )

echo $chr:$start-$end  $gene

samtools mpileup -r  $chr:$start-$end   -ugf /home/lizhixin/references/genome/hg19/hg19.fa /home/lizhixin/project2/WGS/cfDNA/fastq/zstar.sorted.bam  | bcftools call -vmO z -o $gene.vcf.gz

done

