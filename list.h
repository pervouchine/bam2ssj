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

struct list_el {
   int id;
   int count53[2];
   struct list_el *next;
};

typedef struct list_el list_element;


void update_count(list_element **pptr, int id, int side);
void destroy_list(list_element *ptr);
