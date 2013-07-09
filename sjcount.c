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

#define MAXFILEBUFFLENGTH 1000
#define ARRAY_MARGIN 2
#define INFTY 65535
#define WINDOW 1
#define MIN2(A,B) ((A<B) ? (A) : (B))

#define MIN2MAX(A,B) ((A)<(B) ? (double)(A)/(B) : (double)(B)/(A))

class site {
  public:
    int pos;
    int count[2];
    site *next;

    site(int end) {
	pos = end;
	count[0] = count[1] = 0;
	next = NULL;
    }
};

class junction {
  public:
    int pos;
    site *partner;
    junction* next;

    junction(int beg, int end) {
	pos = beg;
	partner = new site(end);
	next = NULL;
    }
};

void update_site(site **ptr, int end, int strand, int count) {
    while(*ptr != NULL && (*ptr)->pos < end) {
        ptr = &((*ptr)->next);
    }
    if(*ptr != NULL && (*ptr)->pos == end) {
    }
    else {
        site *next = (*ptr);
        (*ptr) = new site(end);
        (*ptr)->next = next;
    }
    (*ptr)->count[strand]+=count;
}

void update_jnxn(junction **ptr, int beg, int end, int strand, int count) {
    while(*ptr != NULL && (*ptr)->pos < beg) {
    	ptr = &((*ptr)->next);
    }
    if(*ptr != NULL && (*ptr)->pos == beg) {
	update_site(&((*ptr)->partner), end, strand, count);
    }
    else {
        junction *next = (*ptr);
	(*ptr) = new junction(beg, end);
	(*ptr)->next = next;
	(*ptr)->partner->count[strand]+=count;
    }
}


int main(int argc,char* argv[]) {
    time_t timestamp, current_time;
    int i,j,k,s,a,n;

    FILE *input_file;
    FILE *ssj_file;
    FILE *ssc_file;
    FILE* log_file=stderr;

    bamFile bam_input;
    bam_header_t *header;
    bam1_t* b;
    bam1_core_t *c;
    uint32_t *cigar;

    char bam_file_name[MAXFILEBUFFLENGTH]="";
    char ssj_file_name[MAXFILEBUFFLENGTH]="";
    char ssc_file_name[MAXFILEBUFFLENGTH]="";
    char log_file_name[MAXFILEBUFFLENGTH]="";

    char buff[MAXFILEBUFFLENGTH];
    char chr[MAXFILEBUFFLENGTH];
    int beg, end, pos, prev_pos, offset; 
    int ref_id, ref_id_prev;
    int read_type, mapped_strand;
    int flag;

    int output_strand[2] = {1, -1};

    int max_intron_length = 0;
    int min_intron_length = 50;
    int rev_compl[2] = {1,0};
    int margin = 4;
    int limit_counts = 0;
    int verbose = 1;

    int n_reads = 0;

    junction **root_junction, ***curr_junction;
    site **root_site, ***curr_site;


    timestamp = time(NULL);

    if(argc==1) {
        fprintf(stderr, "This utility counts splice junctions and reads covering exon boundaries; annotation-agnostic version\n");
        fprintf(stderr, "%s -bam <bam_file> [-ssj <junctions_output>] [-ssc <boundaries_output>] [-log <log_file>] ",argv[0]);
	fprintf(stderr, "[-maxlen <max_intron_length>] [-minlen <min_intron_length>] [-margin <length>] [-v suppress verbose output] [-read1 0/1] [-read2 0/1]\n");
        fprintf(stderr, "Type %s -h for more info\n",argv[0]);
        exit(1);
    }

    for(i=1;i<argc;i++) {
	if(strcmp(argv[i], "-bam") == 0) sscanf(argv[++i], "%s", &bam_file_name[0]);
	if(strcmp(argv[i], "-ssj") == 0) sscanf(argv[++i], "%s", &ssj_file_name[0]);
        if(strcmp(argv[i], "-ssc") == 0) sscanf(argv[++i], "%s", &ssc_file_name[0]);
        if(strcmp(argv[i], "-log") == 0) sscanf(argv[++i], "%s", &log_file_name[0]);

        if(strcmp(argv[i], "-read1") == 0) sscanf(argv[++i], "%i", &rev_compl[0]);
        if(strcmp(argv[i], "-read2") == 0) sscanf(argv[++i], "%i", &rev_compl[1]);

	if(strcmp(argv[i], "-lim") == 0)    sscanf(argv[++i], "%i", &limit_counts);
	if(strcmp(argv[i], "-minlen") == 0) sscanf(argv[++i], "%i", &min_intron_length);
	if(strcmp(argv[i], "-maxlen") == 0) sscanf(argv[++i], "%i", &max_intron_length);
	if(strcmp(argv[i], "-margin") == 0) sscanf(argv[++i], "%i", &margin);

	if(strcmp(argv[i], "-v") == 0) verbose = 0;

        if(strcmp(argv[i], "-h") ==0 ) {
            fprintf(stderr, "Input:  (1) sorted BAM file\n");
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "\t-maxlen <upper limit on intron length>; 0 = no limit (default=%i)\n",max_intron_length);
            fprintf(stderr, "\t-minlen <lower limit on intron length>; 0 = no limit (default=%i)\n",min_intron_length);
            fprintf(stderr, "\t-margin <length> minimum number of flanking nucleotides in the read in order to support SJ or EB, (default=%i)\n",margin);
            fprintf(stderr, "\t-read1 0/1, reverse complement read1 no/yes (default=%i)\n",rev_compl[0]);
            fprintf(stderr, "\t-read2 0/1, reverse complement read2 no/yes (default=%i)\n",rev_compl[1]);
            fprintf(stderr, "Output: (1) Junction counts, tab-delimited  (default=stdout)\n");
            fprintf(stderr, "\tColumns are: chr, begin, end, counts (+strand), counts(-strand)\n");
            fprintf(stderr, "Output: (2) Boundary counts, tab-delimited  (default=stdout)\n");
            fprintf(stderr, "\tColumns are: chr, position, counts (+strand), counts(-strand)\n");
            exit(1);
        }


    }

    if(log_file_name[0]==0) {
	log_file = stderr;
    }
    else {
	log_file = fopen(log_file_name,"w");
	if(log_file == NULL) log_file = stderr;
    }

    if(bam_file_name[0]==0) {
	fprintf(log_file,"Bam not specified, exiting\n");
	exit(1); 
    }

    if(ssj_file_name[0]==0) {
	fprintf(log_file,"[Warning: junction output set to stdout]\n");
	ssj_file = stdout;
    }
    else {
	ssj_file = fopen(ssj_file_name,"w");
	if(ssj_file == NULL) {
	    fprintf(log_file,"[Warning: junction output set to stdout]\n");
            ssj_file = stdout;
	} else {
	    fprintf(log_file,"[Junction counts: >%s]\n",ssj_file_name);
	}
    }

    if(ssc_file_name[0]==0) {
	fprintf(log_file,"[Warning: boundary output skipped]\n");
    }
    else {
        ssc_file = fopen(ssc_file_name,"w");
        if(ssc_file == NULL) {
            fprintf(log_file,"[Warning: boundary output set to stdout]\n");
            ssc_file = stdout;
        } else {
            fprintf(log_file,"[Boundary counts: >%s]\n",ssc_file_name);
        }
    }

    if(max_intron_length>0) {
	fprintf(log_file,"[Warning: set max intron length=%i]\n", max_intron_length);
    }

    if(min_intron_length>0) {
        fprintf(log_file,"[Warning: set min intron length=%i]\n", min_intron_length);
    }

    if(margin>0) {
	fprintf(log_file,"[Warning: read margin set to %i]\n", margin);
    }

    for(s = 0; s < 2; s++) {
	if(rev_compl[s]) fprintf(log_file,"[Warning: will take reverse complement of read %i]\n", s+1);
    }


    //*****************************************************************************************************************************//
   
    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    if(bam_input == NULL || header == NULL) {
        fprintf(log_file,"BAM can't be opened or contains no header, exiting\n");
        exit(1);
    }

    root_junction = (junction**)  malloc( sizeof(junction*)  * (header->n_targets + ARRAY_MARGIN) );
    curr_junction = (junction***) malloc( sizeof(junction**) * (header->n_targets + ARRAY_MARGIN) );

    root_site = (site**)  malloc( sizeof(site*)  * (header->n_targets + ARRAY_MARGIN) );
    curr_site = (site***) malloc( sizeof(site**) * (header->n_targets + ARRAY_MARGIN) );

    for(i=0; i < header->n_targets; i++) {
	root_junction[i] = NULL;
	curr_junction[i] = &root_junction[i];
	root_site[i] = NULL;
	curr_site[i] = &root_site[i];
    }

    b = bam_init1();
    k = 0;
    ref_id = -1;
    n_reads = 0;
    while(bam_read1(bam_input, b)>=0) {
        c   = &b->core;
        if(c->tid < 0 || c->tid >= header->n_targets) continue;

        ref_id_prev = ref_id;
        ref_id = c->tid;

        cigar = bam1_cigar(b);

	s = ((c->flag & BAM_FREVERSE)>0);

	mapped_strand = (c->flag & BAM_FREAD1) ? (s + rev_compl[0]) & 1 : (s + rev_compl[1]) & 1;

	pos = beg = c->pos + 1;
	end = bam_calend(c, cigar);

	while((*curr_junction[ref_id])!=NULL && (*curr_junction[ref_id])->pos < beg) curr_junction[ref_id] = &((*curr_junction[ref_id])->next);
        while((*curr_site[ref_id])!=NULL && (*curr_site[ref_id])->pos < beg) curr_site[ref_id] = &((*curr_site[ref_id])->next);

	if(ref_id != ref_id_prev  && ref_id_prev >= 0) {
	    progressbar(1, 1, header->target_name[ref_id_prev], verbose);
	    k=0;
	}

	for(;k<beg;k++) progressbar(k, header->target_len[ref_id], header->target_name[ref_id], verbose);

        for(i = 0; i < c->n_cigar; i++) {
	    offset = cigar[i] >> 4;
	    prev_pos = pos;
	    switch(cigar[i] & 0x0F) {
	    	case BAM_CMATCH: 	pos += offset;  // match to the reference
				 	break;
		case BAM_CINS:		pos += 0;	// insertion to the reference, pointer stays unchanged
					break;
		case BAM_CDEL:		pos += offset;	// deletion from the reference (technically the same as 'N') pointer moves
                                        break;
                case BAM_CREF_SKIP:     pos += offset;
					if(prev_pos - beg < margin) break;
					if(end - pos < margin) break;
					if(offset < min_intron_length && min_intron_length > 0) continue;
					if(offset > max_intron_length && max_intron_length > 0) break;
					update_jnxn(curr_junction[ref_id], prev_pos - 1, pos, mapped_strand, 1);
					update_site(curr_site[ref_id], prev_pos - 1, mapped_strand, 0);
					update_site(curr_site[ref_id], pos, mapped_strand, 0);
				 	break;
		case BAM_CSOFT_CLIP:
		case BAM_CHARD_CLIP:
		case BAM_CPAD:		break;
		default:		break;
	    }
       	}
	n_reads++;
	if(n_reads>limit_counts && limit_counts>0) break;
    }
    if(verbose) progressbar(1, 1, header->target_name[ref_id_prev], verbose); 

    for(i = 0; i < header->n_targets; i++) {
	junction* ptr = root_junction[i];
	while(ptr != NULL) {
	    site* qtr = ptr->partner;
	    while(qtr != NULL) {
		fprintf(ssj_file, "%s\t%i\t%i\t%i\t%i\n", header->target_name[i], ptr->pos, qtr->pos, qtr->count[0], qtr->count[1]);
		qtr = qtr->next;
	    }
	    ptr= ptr->next;
	}
    }

    fclose(ssj_file);
    bam_header_destroy(header);
    bam_close(bam_input);
    bam_destroy1(b);

    //*********************************************************************************************************************//

    if(ssc_file == NULL || ssc_file_name[0]==0) {
    	current_time = time(NULL);
    	fprintf(log_file,"Completed in %1.0lf seconds\n",difftime(current_time,timestamp));
    	return 0;
    }

    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    for(i = 0; i < header->n_targets; i++) {
        site *qtr = root_site[i];
        while(qtr != NULL) {
	    qtr->count[0] = qtr->count[1] = 0;
            qtr = qtr->next;
        }
	curr_site[i] = &root_site[i];
    }

    b = bam_init1();
    k = 0;
    ref_id = -1;
    n_reads = 0;
    while(bam_read1(bam_input, b)>=0) {
        c   = &b->core;
        if(c->tid < 0 || c->tid >= header->n_targets) continue;

        ref_id_prev = ref_id;
        ref_id = c->tid;

        cigar = bam1_cigar(b);

        s = ((c->flag & BAM_FREVERSE)>0);

        mapped_strand = (c->flag & BAM_FREAD1) ? (s + rev_compl[0]) & 1 : (s + rev_compl[1]) & 1;

        beg = c->pos + 1;
        end = bam_calend(c, cigar);

        while((*curr_site[ref_id])!=NULL && (*curr_site[ref_id])->pos < beg) curr_site[ref_id] = &((*curr_site[ref_id])->next);

        if(ref_id != ref_id_prev  && ref_id_prev >= 0) {
            progressbar(1, 1, header->target_name[ref_id_prev], verbose);
            k=0;
        }

        for(;k<beg;k++) progressbar(k, header->target_len[ref_id], header->target_name[ref_id], verbose);

	flag = 1;
    	pos = beg;
        for(i = 0; i < c->n_cigar; i++) {
            offset = cigar[i] >> 4;
            switch(cigar[i] & 0x0F) {
                case BAM_CMATCH:        pos += offset;  // match to the reference
                                        break;
                case BAM_CINS:          pos += 0;       // insertion to the reference, pointer stays unchanged
                                        break;
                case BAM_CDEL:          pos += offset;  // deletion from the reference (technically the same as 'N') pointer moves
                                        break;
                case BAM_CREF_SKIP:
                case BAM_CSOFT_CLIP:
                case BAM_CHARD_CLIP:
                case BAM_CPAD:           
                default:                flag=0;
					break;
            }
        }
        end = pos;

	if(flag) {
	    site *qtr = (*curr_site[ref_id]);
	    while(qtr != NULL && qtr->pos < end) {
		if(qtr->pos > beg + margin && qtr->pos < end - margin) qtr->count[mapped_strand]++;
		qtr = qtr->next;
	    }
	}

        n_reads++;
        if(n_reads>limit_counts && limit_counts>0) break;
    }
    if(verbose) progressbar(1, 1, header->target_name[ref_id_prev], verbose);

    for(i = 0; i < header->n_targets; i++) {
        site *qtr = root_site[i];
        while(qtr != NULL) {
            fprintf(ssc_file, "%s\t%i\t%i\t%i\n", header->target_name[i], qtr->pos, qtr->count[0], qtr->count[1]);
            qtr = qtr->next;
        }
    }
    fclose(ssc_file);
    bam_header_destroy(header);
    bam_close(bam_input);
    bam_destroy1(b);

    current_time = time(NULL);
    fprintf(log_file,"Completed in %1.0lf seconds\n",difftime(current_time,timestamp));
    return 0;

}
