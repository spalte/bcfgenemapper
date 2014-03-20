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
#include <ctype.h>
#include "genemapper.h"

int32_t exon_range_length(exon_range_t exon) {
    if (exon_range_strand(exon) == plusstrand) {
        return (exon.end - exon.start) + 1;
    } else {
        return (exon.start - exon.end) + 1;
    }
}


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
    
    fscanf(fp, "Sequence:\n");
    
    long sequenceFilePoss = ftell(fp);
    int32_t sequenceLength = 0;
    while (isspace(fgetc(fp)) == 0) {
        sequenceLength++;
    }
    fseek(fp, sequenceFilePoss, SEEK_SET);
    
    newGeneMapper->referenceGenome = (char *)malloc(sequenceLength + 1);

#warning dynamically allocate essential position array length
    newGeneMapper->essentialPositions = (int32_t *)malloc(sizeof(int32_t) * 1000);
    memset(newGeneMapper->essentialPositions, 0, sizeof(int32_t) * 1000);
    
    fscanf(fp, "%s\n", newGeneMapper->referenceGenome);

    int position = 0;
    for (newGeneMapper->essentialPositionCount = 0; fscanf(fp, "%d\n", &position) == 1;) {
        newGeneMapper->essentialPositions[newGeneMapper->essentialPositionCount] = position;
        newGeneMapper->essentialPositionCount++;
    }
    
    int32_t totalExonLength = 0;
    int i = 0;
    for (i = 0; i < newGeneMapper->exonCount; i++) {
        totalExonLength += exon_range_length(newGeneMapper->exons[i]);
    }
    
    if (strlen(newGeneMapper->referenceGenome) != totalExonLength) {
        fprintf(stderr, "***WARNING*** The Genomic Reference sequence has a length of %d when the total length of the exons is %d.\n", (int)strlen(newGeneMapper->referenceGenome), (int)totalExonLength);
    }
    
    return newGeneMapper;
}

void gene_mapper_add_exon(gene_mapper_t* geneMapper, exon_range_t exon)
{
    if (geneMapper->exonCount == geneMapper->exonsAllocated) {
        geneMapper->exons = (exon_range_t *)realloc(geneMapper->exons, sizeof(exon_range_t) * (geneMapper->exonsAllocated * 2));
        geneMapper->exonsAllocated *= 2;
        memset(geneMapper->exons + geneMapper->exonCount, 0, sizeof(exon_range_t) * (geneMapper->exonsAllocated - geneMapper->exonCount));
    }
    geneMapper->exons[geneMapper->exonCount] = exon;
    geneMapper->exonCount++;
}

void gene_mapper_destroy(gene_mapper_t* geneMapper)
{
    free(geneMapper->exons);
    free(geneMapper->referenceGenome);
    free(geneMapper->essentialPositions);
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
    
    int32_t totalLength = 0;
    for (i = 0; i < geneMapper->exonCount; i++) {
        exon_range_t exon = geneMapper->exons[i];
        fprintf(fp, "    %8d  %8d ", (int)exon.start, (int)exon.end);
        fprintf(fp, "  Length: %5d  (%c)strand\n", exon_range_length(exon), exon_range_strand(exon));
        totalLength += exon_range_length(exon);
    }
    fprintf(fp, "Total Length: %d\n", totalLength);
}

















