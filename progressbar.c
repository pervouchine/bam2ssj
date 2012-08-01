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

#include"progressbar.h"

void progressbar(unsigned long current, unsigned long last, char *inscription, int verbose) {
    int i,l,k;
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    unsigned long m;

    m = (int)(w.ws_col*0.8);

    if(last==0) return;

    if(!verbose) return;

    if((m*(current-1))/last < (m*current)/last) {
        k = (m*current)/last;
        fprintf(stderr,"%c%s\t[",13,inscription);
        for(i=0;i<k;i++) fprintf(stderr,"=");
        if(k<m) fprintf(stderr,">");
        for(i++;i<m;i++) fprintf(stderr," ");
        fprintf(stderr,"] %2.1lf%% ",(double)100*current/last);
    }
    if(current==last) fprintf(stderr,"\n");
}
