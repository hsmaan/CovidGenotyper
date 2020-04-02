SNPEFF=/usr/local/bin/snpEff/snpEff.jar
  
cd ../data

DATE=`date +"%Y-%m-%d"`

snp-sites -cv -o 12_aligned_$DATE.vcf dec_aligned_fasta_filtered*;

vcftools --vcf 12_aligned_$DATE.vcf --freq --out 12_var_freq_$DATE;

vcftools --vcf 12_aligned_$DATE.vcf --hap-r2 --out 12_var_linkage_$DATE;


