cd ../data

DATE=`date +"%Y-%m-%d"`

rm *.vcf *.frq *.log *.hap.ld

snp-sites -v -o aligned_$DATE.vcf dec_aligned_plus_ref_filtered*;

java -jar /usr/local/bin/snpEff/snpEff.jar COVID aligned_$DATE.vcf > aligned_annotated_$DATE.vcf;

vcftools --vcf aligned_annotated_$DATE.vcf --freq --out aligned_annotated_freq_$DATE;

vcftools --vcf aligned_annotated_$DATE.vcf --hap-r2 --out aligned_annotated_linkage_$DATE;

rm aligned_$DATE.vcf
