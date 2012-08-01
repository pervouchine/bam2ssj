# This utility converts .gtf genomic annotation to cps file containing exon boundaries
# $1 = gtf, stdout = cps
awk '$3=="exon"' $1 | perl -e 'while(<>) { if(m/G(\d+)/) {@a = split /\t/; print "$a[0]\t$a[3]\t$a[6]\t$1\n$a[0]\t$a[4]\t$a[6]\t$1\n"}}' | sort -k 2,2n | uniq 
