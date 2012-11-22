//	Copyright 2012 Dmitri Pervouchine (dp@crg.eu), Lab Roderic Guigo
//	Bioinformatics and Genomics Group @ Centre for Genomic Regulation
//	Parc de Recerca Biomedica: Dr. Aiguader, 88, 08003 Barcelona
//	
//	This file is a part of the 'bam2ssj' package.
//	'bam2ssj' package is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//	
//	'bam2ssj' package is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with 'bam2ssj' package.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/ioctl.h>
#include <bam.h>
#include "progressbar.h"
#include "list.h"

#define MAXFILEBUFFLENGTH 1000
#define ARRAY_MARGIN 2
#define INFTY 65535
#define WINDOW 1
#define MIN2(A,B) ((A<B) ? (A) : (B))

#define RT_OVRLAP 0
#define RT_GENOME 1
#define RT_KJUNCT 2
#define RT_UJUNCT 3
#define RT_OTHER  4

#define MIN2MAX(A,B) ((A)<(B) ? (double)(A)/(B) : (double)(B)/(A))

#define BAM_UNIQUE_MAP 0x800

#define N_READ_TYPES 5

char read_type_descr[N_READ_TYPES][MAXFILEBUFFLENGTH] = { "genomic overlapping sj", "genomic non-overlapping sj (ignored)", 
 "known junctions", "unknown junctions (ignored)", "ignored by other reasons"};


typedef struct {
    int pos;
    int label;
    int count00[2];
    int count5X[2];
    int countX3[2];
    list_element* junctions;
}   splice_site;

int verbose=1;

int main(int argc,char* argv[]) {
    time_t timestamp, current;
    int i,j,k;
    int a,n;
    char *pc;

    FILE *input_file;
    FILE *output_file;

    bamFile bam_input;
    bam_header_t *header;
    bam1_t* b;
    bam1_core_t *c;


    char cps_file_name[MAXFILEBUFFLENGTH]="";
    char bam_file_name[MAXFILEBUFFLENGTH]="";
    char out_file_name[MAXFILEBUFFLENGTH]="";

    char buff[MAXFILEBUFFLENGTH];
    char chr[MAXFILEBUFFLENGTH];
    int beg, beg_prev, end, pos, offset; 
    int ref_id, ref_id_prev, label;
    int s, side;
    int read_type, mapped_strand;
    char ch;

    int limit_counts = 0;

    int* contig_count[2];
    int* contig_index[2];
    splice_site** contig_sites[2];

    long int n_reads[N_READ_TYPES][2];

    long int n_total_reads = 0;
    long int n_skipped_reads = 0;

    int max_intron_length=0;
    int min_intron_length=0;
    int ignore_gene_labels = 0;
    int stranded = 1;
    int rev_compl[2] = {1,0};

    int other_end, donor_id, acceptor_id;

    int *cigar;
    int flagged = 0;


    /** reading input from the command line **/

    timestamp = time(NULL);

    if(argc==1) {
	fprintf(stderr, "BAM2SSJ is the utility for fast counting reads covering splice junctions\nCommand line use:\n");
        fprintf(stderr, "%s -cps <cps_file> -bam <bam_file> [-out <out_file>] [-maxlen <max_intron_length>] [-minlen <min_intron_length>] ",argv[0]);
	fprintf(stderr, "[-v suppress verbose output] [-read1 0/1] [-read2 0/1] [-g ignore gene labels] [-u unstranded] [-f count reads flagged 0x800 only]\ntype %s -h for more info\n",argv[0]);
        exit(1);
    }

    for(i=1;i<argc;i++) {
        pc = argv[i];
        if(*pc == '-') {
            if(strcmp(pc+1,"cps") == 0) sscanf(argv[++i], "%s", &cps_file_name[0]);
	    if(strcmp(pc+1,"bam") == 0) sscanf(argv[++i], "%s", &bam_file_name[0]);
	    if(strcmp(pc+1,"out") == 0) sscanf(argv[++i], "%s", &out_file_name[0]);

            if(strcmp(pc+1,"read1") == 0) sscanf(argv[++i], "%i", &rev_compl[0]);
            if(strcmp(pc+1,"read2") == 0) sscanf(argv[++i], "%i", &rev_compl[1]);

	    if(strcmp(pc+1,"lim") == 0) sscanf(argv[++i], "%i", &limit_counts);
	    if(strcmp(pc+1,"minlen") == 0) sscanf(argv[++i], "%i", &min_intron_length);
	    if(strcmp(pc+1,"maxlen") == 0) sscanf(argv[++i], "%i", &max_intron_length);

	    if(strcmp(pc+1,"v") == 0) verbose = 0;
	    if(strcmp(pc+1,"g") == 0) ignore_gene_labels = 1;
	    if(strcmp(pc+1,"u") == 0) stranded = 0;
	    if(strcmp(pc+1,"f") == 0) flagged = 1;

	    if(strcmp(pc+1,"h") ==0 ) {
		fprintf(stderr, "Input:  (1) sorted BAM file\n");
		fprintf(stderr, "\t(2) CPS (chromosome-position-strand) tab-delimited file sorted by position (chr1 100 + etc)\n\n");
        	fprintf(stderr, "\tIn order to get CPS file from gtf, use the utility gtf2cps.sh\n");
        	fprintf(stderr, "\tImportant: CPS must be sorted by position ONLY!\n\n");
        	fprintf(stderr, "\tIf the 4th column contains (a numeric) gene label then only splice junctions within the same gene will be considered (unless the '-g' option is active)\n");
		fprintf(stderr, "\tThe utility to generate CPS with gene labels is gtf2cps_with_gene_id.sh (or update the script accordingly if you are using genome other than human)\n\n");
		fprintf(stderr, "Options:\n");
        	fprintf(stderr, "\t-maxlen <upper limit on intron length>; 0 = no limit (default=%i)",max_intron_length);
		fprintf(stderr, "\t-minlen <lower limit on intron length>; 0 = no limit (default=%i)",min_intron_length);
        	fprintf(stderr, "\t-read1 0/1, reverse complement read1 no/yes (default=%i)\n",rev_compl[0]);
        	fprintf(stderr, "\t-read2 0/1, reverse complement read2 no/yes (default=%i)\n",rev_compl[1]);
        	fprintf(stderr, "\t-g ignore gene labels (column 4 of cps), default=%s\n", ignore_gene_labels ? "ON" : "OFF");
        	fprintf(stderr, "\t-u ignore strand (all reads map to the correct strand), default=%s\n", stranded ? "OFF" : "ON");
		fprintf(stderr, "\t-f count reads flagged 0x800 only (uniquely mapped reads), default=%s\n", flagged ? "ON" : "OFF");
		fprintf(stderr, "Output: tab-delimited  (default=stdout)\n");
        	fprintf(stderr, "\tColumn 1 is splice_junction_id\n");
        	fprintf(stderr, "\tColumns 2-6 are counts of 53, 5X, X3, 50, and 03 reads for the correct (annotated) strand\n");
        	fprintf(stderr, "\tColumns 7-11 are similar counts for the incorrect (opposite to annotated) strand\n");
		fprintf(stderr, "Descriptive read statistics are reported to stderr\n");
		exit(1);
	    }
	}
    }

    if(bam_file_name[0]==0) {
	fprintf(stderr,"Bam not specified, exiting\n");
	exit(1); 
    }

    if(cps_file_name[0]==0) {
        fprintf(stderr,"Input not specified, exiting\n");
        exit(1);
    }

    if(out_file_name[0]==0) {
	fprintf(stderr,"[Warning: output set to stdout]\n");
	output_file = stdout;
    }
    else {
	output_file = fopen(out_file_name,"w");
	if(output_file == NULL) {
	    fprintf(stderr,"[Warning: output set to stdout]\n");
            output_file = stdout;
	}
    }

    if(max_intron_length>0) {
	if(verbose) fprintf(stderr,"[Warning: set max intron length=%i]\n",max_intron_length);
    }

    if(ignore_gene_labels) {
	if(verbose) fprintf(stderr,"[Warning: ignoring gene labels (column 4)]\n");
    }

    if(flagged) {
	if(verbose) fprintf(stderr,"[Warning: only look at reads flagged 0x800\n");
    }

    if(verbose) {
	for(s = 0; s < 2; s++) if(rev_compl[s]) fprintf(stderr,"[Warning: take reverse complement of read %i]\n", s+1);
	fprintf(stderr,"[Warning: stranded = %s]\n", stranded ? "TRUE" : "FALSE (always correct strand)");
	if(ignore_gene_labels) fprintf(stderr,"[Warning: ignore gene labels (column 4)]\n");
    }


    for(i = 0; i < N_READ_TYPES; i++) for(s = 0; s < 2; s++) n_reads[i][s] = 0;

    /** initatializing BAM and header **/
   
    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    if(bam_input == NULL || header == NULL) {
        fprintf(stderr,"BAM can't be opened or contains no header, exiting\n");
        exit(1);
    }

    /** reading input from CPS **/

    input_file = fopen(cps_file_name, "r");
    if(input_file == NULL) {
	fprintf(stderr,"CPS can't be opened, exiting\n");
        exit(1);
    }

    /** populating gene structure arrays **/

    for(s = 0; s < 2; s++) {
    	contig_count[s] = (int*) malloc(sizeof(int) * (header->n_targets + ARRAY_MARGIN));
    	contig_index[s] = (int*) malloc(sizeof(int) * (header->n_targets + ARRAY_MARGIN));
    	contig_sites[s] = (splice_site**) malloc(sizeof(splice_site*) * (header->n_targets + ARRAY_MARGIN));

    	if(contig_count[s] == NULL || contig_sites[s] == NULL || contig_index[s] == NULL) {
	    fprintf(stderr, "Not enought memory, exiting\n");
            exit(1);
    	}
    }

    for(s = 0; s < 2; s++)
        for(i=0; i < header->n_targets; i++) 
	    contig_count[s][i] = contig_index[s][i] = 0;

    if(verbose) fprintf(stderr, "Reading %s pass1", cps_file_name);
    while(fgets(buff, MAXFILEBUFFLENGTH, input_file)) {
	sscanf(buff, "%s %*i %c", &chr[0], &ch);
	bam_parse_region(header, chr, &i, &beg, &end);
	s = (ch == '+' ? 0 : 1);
	if(i < header->n_targets) contig_count[s][i]++;
    }

    for(s = 0; s < 2; s++) {
    	for(i = 0;i < header->n_targets; i++) {
	    contig_sites[s][i] = (splice_site*) malloc(sizeof(splice_site) * (contig_count[s][i] + ARRAY_MARGIN));
	    if(contig_sites[s][i] == NULL) {
	    	fprintf(stderr, "Not enought memory, exiting\n");
            	exit(1);
	    }
	}
    }
    if(verbose) fprintf(stderr, "\n");

    if(verbose) fprintf(stderr, "Reading %s pass2",cps_file_name);
    fseek(input_file, 0, SEEK_SET);
    while(fgets(buff, MAXFILEBUFFLENGTH, input_file)) {
        sscanf(buff, "%s %i %c %i", &chr[0], &pos, &ch, &label);
	bam_parse_region(header, chr, &i, &beg, &end);
	s = (ch == '+' ? 0 : 1);
	if(i < header->n_targets) {
	    if(contig_index[s][i]>0) {
		if(pos < contig_sites[s][i][contig_index[s][i]-1].pos) {
		    fprintf(stderr, "Splice sites weren't sorted, exiting\n");
		    exit(1);
		}
	    }
	    contig_sites[s][i][contig_index[s][i]].pos = pos;
	    contig_sites[s][i][contig_index[s][i]].label = ignore_gene_labels ? 0 : label;
	    for(side = 0; side < 2; side++) {
                contig_sites[s][i][contig_index[s][i]].count00[side] = 0;
                contig_sites[s][i][contig_index[s][i]].count5X[side] = 0;
                contig_sites[s][i][contig_index[s][i]].countX3[side] = 0;
		contig_sites[s][i][contig_index[s][i]].junctions = NULL;
	    }
	    contig_index[s][i]++;
	}
    }
    if(verbose) fprintf(stderr, "\n");

    for(s = 0; s < 2; s++)
    	for(i = 0;i < header->n_targets; i++) 
	    contig_index[s][i] = 0;

    /** analysis starts here **/

    b = bam_init1();
    k = 0;
    ref_id_prev = -1;
    beg_prev = -1;
    while(bam_read1(bam_input, b)>=0) {
        c   = &b->core;
	ref_id = c->tid;
	if(ref_id<0) continue;

	if(flagged && (c->flag & 0x800 == 0)) {
	    n_skipped_reads++;
	    continue;
	}

        if(stranded && ((c->flag & BAM_FREAD1) && (c->flag & BAM_FREAD2) || !(c->flag & BAM_FREAD1) && !(c->flag & BAM_FREAD2))) {
            n_skipped_reads++;
            continue;
        }

        cigar = bam1_cigar(b);

	if(ref_id != ref_id_prev  && ref_id_prev >= 0) {
	    if(contig_index[0][ref_id_prev] + contig_index[1][ref_id_prev] < contig_count[0][ref_id_prev] + contig_count[1][ref_id_prev]) 
		progressbar(1,1,header->target_name[ref_id_prev], verbose);
	    beg_prev = -1;
	}
	if(ref_id < ref_id_prev) {
	    fprintf(stderr,"BAM file wasn't sorted, exiting\n");
            exit(1);
	}

	ref_id_prev = ref_id;

	beg = c->pos + 1;
	if(beg < beg_prev) {
	    fprintf(stderr,"BAM file wasn't sorted, exiting\n");
	    exit(1);
	}
	beg_prev = beg;

	s = ((c->flag & BAM_FREVERSE)>0);
	mapped_strand = (c->flag & BAM_FREAD1) ? (s + rev_compl[0]) & 1 : (s + rev_compl[1]) & 1;

	for(s = 0; s < 1 + stranded; s++) {
            end = beg;
	    side = (s == mapped_strand) ? 0 : 1;
	    side *= stranded;

	    // keep reading until the currect site is on the same chromosome downstream of the read 

	    while(contig_sites[s][ref_id][contig_index[s][ref_id]].pos < beg && contig_index[s][ref_id] < contig_count[s][ref_id]) {
		contig_index[s][ref_id]++;
	    	progressbar(contig_index[0][ref_id]+contig_index[1][ref_id], contig_count[0][ref_id]+contig_count[1][ref_id], header->target_name[ref_id], verbose);
	    }

	    read_type = RT_OTHER;

            if(contig_index[s][ref_id]<contig_count[s][ref_id]) {
	    	// check if the read is a split read and find its other end
	    	read_type = RT_GENOME;
            	for(i = 0; i < c->n_cigar; i++) {
	    	    offset = cigar[i] >> 4;
	    	    switch(cigar[i] & 0x0F) {
		    	case BAM_CMATCH: 	end += offset;  // match to the reference
					 	break;
		    	case BAM_CINS:		end += 0;	// insertion to the reference, pointer stays unchanged
						break;
		    	case BAM_CDEL:		end += offset;	// deletion from the reference (technically the same as 'N') pointer moves
						break; 
		    	case BAM_CREF_SKIP:	other_end = end + offset;
						donor_id = acceptor_id = -INFTY;
						for(j = contig_index[s][ref_id]; contig_sites[s][ref_id][j].pos <= other_end && j < contig_count[s][ref_id];j++) {
						    if(contig_sites[s][ref_id][j].pos - end < min_intron_length && min_intron_length > 0) continue;
						    if(contig_sites[s][ref_id][j].pos - end > max_intron_length && max_intron_length > 0) break;
					    	    if(contig_sites[s][ref_id][j].label == contig_sites[s][ref_id][contig_index[s][ref_id]].label) {
					    	    	if(contig_sites[s][ref_id][j].pos == end - 1)   donor_id = j;
					    	    	if(contig_sites[s][ref_id][j].pos == other_end) acceptor_id = j;
					    	    }
					    	}
						if(donor_id>0 && acceptor_id>0) {
					    	    update_count(&contig_sites[s][ref_id][donor_id].junctions, acceptor_id, side);
					    	    contig_sites[s][ref_id][donor_id].count5X[side]++;
                                            	    contig_sites[s][ref_id][acceptor_id].countX3[side]++;
					    	    read_type = RT_KJUNCT;
						}
						else {
					    	    read_type = RT_UJUNCT;
						}
						end = other_end;
				 		break;
		    	case BAM_CSOFT_CLIP:
		    	case BAM_CHARD_CLIP:
		    	case BAM_CPAD:		break;
		    	default:		read_type = RT_OTHER;
	    	    }
            	}

	    	if(read_type == RT_GENOME) {
	            for(j=contig_index[s][ref_id]; beg<=contig_sites[s][ref_id][j].pos  && contig_sites[s][ref_id][j].pos<end && j<contig_count[s][ref_id]; j++) {
		    	contig_sites[s][ref_id][j].count00[side]++;
		    	read_type = RT_OVRLAP;
		    	k++;
	    	    }
	    	}
	    }

	    n_reads[read_type][side]++;
	}
	n_total_reads++;

	if(k>limit_counts && limit_counts>0) break;

    }

    if(ref_id_prev > 0) {
        if(contig_index[0][ref_id_prev] + contig_index[1][ref_id_prev] < contig_count[0][ref_id_prev] + contig_count[1][ref_id_prev]) 
	    progressbar(1, 1, header->target_name[ref_id_prev], verbose);
    }

    /**** output here *****/

    for(s = 0; s < 2; s++) {
    	for(i = 0; i < header->n_targets; i++) {
	    for(j = 0; j < contig_count[s][i]; j++) {
	    	list_element* ptr = contig_sites[s][i][j].junctions;
	    	while(ptr!=NULL) {
		    fprintf(output_file, "%s_%i_%i_%i", header->target_name[i], contig_sites[s][i][j].pos + 1, contig_sites[s][i][ptr->id].pos - 1,(s == 0) ? 1  : -1);
		    for(side = 0; side < 1 + stranded; side++) {
		    	fprintf(output_file, "\t%i\t%i\t%i\t%i\t%i", ptr->count53[side], 
			   contig_sites[s][i][j].count5X[side] - ptr->count53[side], 
			   contig_sites[s][i][ptr->id].countX3[side] - ptr->count53[side],
			   contig_sites[s][i][j].count00[side],
			   contig_sites[s][i][ptr->id].count00[side]);
		    }
		    fprintf(output_file, "\n");
	    	    ptr= ptr->next;
	    	}
	    }
	}
    }

    for(s = 0; s < 2; s++) {
    	for(i = 0;i < header->n_targets; i++) {
	    for(j = 0; j < contig_count[s][i]; j++) {
	    	destroy_list(contig_sites[s][i][j].junctions);
	    }
            free(contig_sites[s][i]);
        }
        free(contig_sites[s]);
        free(contig_index[s]);
        free(contig_count[s]);
    }

    bam_header_destroy(header);
    bam_close(bam_input);
    bam_destroy1(b);

    fclose(input_file);
    fclose(output_file);

    if(verbose) fprintf(stderr,"\n");

    current = time(NULL);
    fprintf(stderr,"Read statistics for %s\n",bam_file_name);
    for(s = 0; s < 1 + stranded; s++) fprintf(stderr,"%16s", s == 0 ? "correct" : "incorrect");
    if(stranded) fprintf(stderr,"\tmin/max");
    fprintf(stderr,"\tstrand\n");

    for(i = 0; i < N_READ_TYPES; i++) {
	for(s = 0; s < 1 + stranded; s++) fprintf(stderr,"%16li", n_reads[i][s]);
	if(stranded) fprintf(stderr,"\t%1.2lf",MIN2MAX(n_reads[i][1],n_reads[i][0]));
	fprintf(stderr,"\t%s\n",read_type_descr[i]);
    }
    for(s = 0; s < 1 + stranded; s++) fprintf(stderr,"%16li", n_skipped_reads);
    if(stranded) fprintf(stderr,"\t");
    fprintf(stderr,"\tskipped\n");
    for(s = 0; s < 1 + stranded; s++) fprintf(stderr,"%16li", n_total_reads + n_skipped_reads);
    if(stranded) fprintf(stderr,"\t");
    fprintf(stderr,"\ttotal reads\n");
    fprintf(stderr,"Completed in %1.0lf seconds\n",difftime(current,timestamp));

    return 0;
}
