SAMDIR=~/samtools/
GCC=g++

.PHONY: all

all: bam2ssj sjcount

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
	$(GCC) -c progressbar.c 

list.o : list.c list.h
	$(GCC) -c list.c

bam2ssj:	bam2ssj.c progressbar.o list.o $(SAMDIR)libbam.a
	$(GCC) -I $(SAMDIR) bam2ssj.c progressbar.o list.o $(SAMDIR)libbam.a -lz -o bam2ssj

sjcount : sjcount.c progressbar.o $(SAMDIR)libbam.a
	$(GCC) -I $(SAMDIR) sjcount.c progressbar.o $(SAMDIR)libbam.a -lz -o sjcount

clean:
	rm -f -r progressbar.o list.o bam2ssj sjcount

