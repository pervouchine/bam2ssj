SAMDIR=~/samtools/
GCC=gcc
OPTN=-lz

.PHONY: all

all: bam2ssj

EXPORT = bam2ssj-dp-1.6

export:
	mkdir $(EXPORT)/
	cp -f *.c $(EXPORT)/
	cp -f *.h $(EXPORT)/
	cp -f *.sh $(EXPORT)/
	cp -f README $(EXPORT)/
	cp -f VERSION $(EXPORT)/
	cp -f LICENCE $(EXPORT)/
	cp makefile $(EXPORT)/
	tar -cf $(EXPORT).tar $(EXPORT)/
	gzip $(EXPORT).tar
	rm -f -r $(EXPORT)/
	mv $(EXPORT).tar.gz ..

$(SAMDIR)libbam.a:
	# You need to install samtools
	# Get it by svn:
	# svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools
	# enter the dir and type 'make all'
	# don't forget to update the SAMDIR varibale in this makefile
	exit 1	

progressbar.o:	progressbar.c progressbar.h
	gcc -c progressbar.c 

list.o:	list.c list.h
	gcc -c list.c

bam2ssj:	bam2ssj.c progressbar.o list.o $(SAMDIR)libbam.a
	$(GCC) $(OPTN) -I $(SAMDIR) bam2ssj.c progressbar.o list.o $(SAMDIR)libbam.a -o bam2ssj


clean:
	rm -f -r list.o progressbar.o bam2ssj

