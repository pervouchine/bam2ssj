EXAMPLE=example/

all: ${EXAMPLE}output.ssj


${EXAMPLE}gencode.v10.annotation.gtf:
	mkdir -p ${EXAMPLE}
	wget ftp://ftp.sanger.ac.uk/pub/gencode/release_10/gencode.v10.annotation.gtf.gz -O ${EXAMPLE}gencode.v10.annotation.gtf.gz
	gunzip ${EXAMPLE}gencode.v10.annotation.gtf.gz

${EXAMPLE}BjCellPapAlnRep1.bam:
	mkdir -p ${EXAMPLE}
	wget http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/releaseLatest/wgEncodeCshlLongRnaSeqBjCellPapAlnRep1.bam -O ${EXAMPLE}BjCellPapAlnRep1.bam

${EXAMPLE}hg19.cps: ${EXAMPLE}gencode.v10.annotation.gtf
	./gtf2cps_with_gene_id.sh ${EXAMPLE}gencode.v10.annotation.gtf > ${EXAMPLE}hg19.cps

${EXAMPLE}hg19.cps.srt: ${EXAMPLE}hg19.cps
	sort -n -k 2 ${EXAMPLE}hg19.cps > ${EXAMPLE}hg19.cps.srt

${EXAMPLE}output.ssj: ${EXAMPLE}hg19.cps.srt ${EXAMPLE}BjCellPapAlnRep1.bam
	./bam2ssj -cps ${EXAMPLE}hg19.cps.srt -bam ${EXAMPLE}BjCellPapAlnRep1.bam -read1 0 -read2 0 -out ${EXAMPLE}output.ssj


clean:
	rm -f -r ${EXAMPLE}
