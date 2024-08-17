
#IFS=$'\n' 
#IFS=`echo -e "\n"`

awk -F'\t' '{print $2}' ../celltype.txt |sort |uniq >type.txt

for line in `cat type.txt`
do
echo $line

#less celltype.txt |grep $line |awk '{print $1}' >tmp.txt

#zcat pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz |grep -f tmp.txt > ${line}.bed

#awk -F '\t' '{sum[$4]++}END{for(i in sum) print i "\t" sum[i]}' ${line}.bed > ${line}_cellnumber.bed

macs2 callpeak -g 1.87e9 --name ${line} --treatment ${line}.bed --outdir ./ -B --format BED --shift -100 --extsize 200 --nomodel --call-summits --nolambda --keep-dup  all

sort -k1,1 -k2,2n ${line}_treat_pileup.bdg > ${line}_sorted.bdg

cat ${line}_sorted.bdg|grep -E 'chr' > ${line}_sorted_adjust.bdg
bedClip ${line}_sorted_adjust.bdg  /home/chengww/data/project/database/scrna-atac/pbmc_10k/bwfile/hg38.chrom.sizes ${line}_sorted_adjust_clip.bdg

bedGraphToBigWig ${line}_sorted_adjust_clip.bdg /home/chengww/data/project/database/scrna-atac/pbmc_10k/bwfile/hg38.chrom.sizes ${line}.bw

done
