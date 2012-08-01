# This utility converts .gtf genomic annotation to cps file containing exon boundaries
# $1 = gtf, stdout = cps
awk '$3=="exon"' $1 | cut -f 1,4,5,7 | awk '{print $1"\t" $2 "\t" $4 "\n" $1 "\t" $3 "\t" $4}' | sort -k 2,2n | uniq 
