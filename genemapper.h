//
//  genemapper.h
//  vcfparser
//
//  Created by Joël Spaltenstein on 3/13/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#ifndef vcfparser_genemapper_h
#define vcfparser_genemapper_h

#include <stdio.h>

/* All positions are assumed to be 1-indexed */

typedef struct { // exon range, start and end are inclusive
    int32_t start;
    int32_t end;
} exon_range_t;

static inline exon_range_t exon_range(int32_t start, int32_t end) {exon_range_t exon;exon.start = start; exon.end = end; return exon;}

typedef struct {
    int32_t exonCount;
    exon_range_t* exons;
    int32_t exonsAllocated;
} gene_mapper_t;


gene_mapper_t *gene_mapper_init();

gene_mapper_t *gene_mapper_initWithExons(exon_range_t* exons, int32_t exonCount);
gene_mapper_t *gene_mapper_file_init(FILE *fp);

void gene_mapper_add_exon(gene_mapper_t* geneMapper, exon_range_t exon);
void gene_mapper_destroy(gene_mapper_t* geneMapper);

static inline int32_t gene_mapper_exon_count(gene_mapper_t* geneMapper) {return geneMapper->exonCount;}

int32_t gene_mapper_map_position(gene_mapper_t* geneMapper, int32_t genomePosition); // returns -1 if the position does not map
int32_t gene_mapper_reversemap_position(gene_mapper_t* geneMapper, int32_t genePosition); // returns -1 if the position is out of the range


#endif
