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

#include "list.h"

void create_node(list_element **pptr, int id, int side) {
    list_element *ptr;
    *pptr = (list_element*) malloc(sizeof(list_element));
    ptr = *pptr;
    ptr->id = id;
    ptr->count53[0] = ptr->count53[1] = 0;
    ptr->next = NULL;
    ptr->count53[side]++;
}

void update_count(list_element **pptr, int id, int side) {
    list_element *ptr, *prev;

    ptr = *pptr;
    if(ptr==NULL) {
	create_node(pptr, id, side);
        return;
    }
    else  {
        while(ptr!=NULL) {
            if(ptr->id == id) {
                ptr->count53[side]++;
                return;
            }
            prev = ptr;
            ptr = ptr->next;
        }
	create_node(&prev->next, id, side);
    }
}

void destroy_list(list_element *ptr) { 
    list_element *next;
    while(ptr!=NULL) {
        next = ptr->next;
	free(ptr);
	ptr = next;
    }
}
