//
//  genemapper.c
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/13/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "genemapper.h"

gene_mapper_t *gene_mapper_init()
{
    gene_mapper_t *newGeneMapper = (gene_mapper_t *)malloc(sizeof(gene_mapper_t));
    memset(newGeneMapper, 0, sizeof(gene_mapper_t));
    newGeneMapper->exonsAllocated = 2;
    newGeneMapper->exons = (exon_range_t *)malloc(sizeof(exon_range_t) * 2);
    memset(newGeneMapper->exons, 0, sizeof(exon_range_t) * 2);
    
    return newGeneMapper;
}

gene_mapper_t *gene_mapper_initWithExons(exon_range_t* exons, int32_t exonCount)
{
    gene_mapper_t *newGeneMapper = (gene_mapper_t *)malloc(sizeof(gene_mapper_t));
    memset(newGeneMapper, 0, sizeof(gene_mapper_t));
    
    newGeneMapper->exonCount = exonCount;
    newGeneMapper->exonsAllocated = exonCount;
    newGeneMapper->exons = (exon_range_t *)malloc(sizeof(exon_range_t) * exonCount);
    memcpy(newGeneMapper->exons, exons, sizeof(exon_range_t) * exonCount);
    
    return newGeneMapper;
}

gene_mapper_t *gene_mapper_file_init(FILE *fp)
{
    gene_mapper_t *newGeneMapper = gene_mapper_init();
    int start;
    int end;
    
    while (fscanf(fp, "%d %d\n", &start, &end) == 2) {
        gene_mapper_add_exon(newGeneMapper, exon_range(start, end));
    }
    
    return newGeneMapper;
}

void gene_mapper_add_exon(gene_mapper_t* geneMapper, exon_range_t exon)
{
    if (geneMapper->exonCount == geneMapper->exonsAllocated) {
        geneMapper->exons = realloc(geneMapper->exons, sizeof(exon_range_t) * (geneMapper->exonsAllocated * 2));
        geneMapper->exonsAllocated *= 2;
        memset(geneMapper->exons + geneMapper->exonCount, 0, sizeof(exon_range_t) * (geneMapper->exonsAllocated - geneMapper->exonCount));
    }
    geneMapper->exons[geneMapper->exonCount] = exon;
    geneMapper->exonCount++;
}

void gene_mapper_destroy(gene_mapper_t* geneMapper)
{
    free(geneMapper->exons);
    free(geneMapper);
}

int32_t gene_mapper_map_position(gene_mapper_t* geneMapper, int32_t genomePosition, exon_range_t* exonRangeOut)
{
    int32_t runLength = 0;
    int32_t i;
    
    for (i = 0; i < geneMapper->exonCount; i++) {
        exon_range_t exon = geneMapper->exons[i];
        char plusExon = exon.end >= exon.start;
        int32_t exonLength;
        if (plusExon) {
            exonLength = (exon.end - exon.start) + 1;
            if (genomePosition <= exon.end && genomePosition >= exon.start) {
                if (exonRangeOut) {
                    *exonRangeOut = exon;
                }
                return runLength + (genomePosition - exon.start);
            } else {
                runLength += exonLength;
            }
        } else {
            exonLength = (exon.start - exon.end) + 1;
            if (genomePosition <= exon.start && genomePosition >= exon.end) {
                if (exonRangeOut) {
                    *exonRangeOut = exon;
                }
                return runLength + (exon.start - genomePosition);
            } else {
                runLength += exonLength;
            }
        }
    }
    
    return -1;
}

int32_t gene_mapper_reversemap_position(gene_mapper_t* geneMapper, int32_t genePosition)
{
    int32_t runPosition = genePosition;
    int32_t i;

    for (i = 0; i < geneMapper->exonCount; i++) {
        exon_range_t exon = geneMapper->exons[i];
        char plusExon = exon.end >= exon.start;
        int32_t exonLength;
        if (plusExon) {
            exonLength = (exon.end - exon.start) + 1;
            if (exonLength >= runPosition) {
                return exon.start + runPosition;
            } else {
                runPosition -= exonLength;
            }
        } else {
            exonLength = (exon.start - exon.end) + 1;
            if (exonLength >= runPosition) {
                return exon.start - runPosition;
            } else {
                runPosition -= exonLength;
            }
        }
    }
    
    return -1;
}

void gene_mapper_print_exons(gene_mapper_t* geneMapper, FILE *fp)
{
    int32_t i;

    fprintf(fp, "There %s %d exons in the exon file.\n", geneMapper->exonCount == 1?"is":"are", (int)geneMapper->exonCount);
    
    for (i = 0; i < geneMapper->exonCount; i++) {
        exon_range_t exon = geneMapper->exons[i];
        fprintf(fp, "    %8d  %8d ", (int)exon.start, (int)exon.end);
        if (exon_range_strand(exon) == plusstrand) {
            fprintf(fp, "  Length: %5d  (%c)strand\n", (int)(exon.end - exon.start) + 1, exon_range_strand(exon));
        } else {
            fprintf(fp, "  Length: %5d  (%c)strand\n", (int)(exon.start - exon.end) + 1, exon_range_strand(exon));
        }
    }
}



















