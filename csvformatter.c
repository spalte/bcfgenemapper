//
//  csvformater.c
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/15/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#include <stdio.h>

#include "csvformatter.h"
#include "main.h"

static const char *emptyString = "";

csv_formatter_sample_t *csv_formatter_sample_init(const char *sampleName, char allele)
{
    csv_formatter_sample_t *newSample = malloc(sizeof(csv_formatter_sample_t));
    memset(newSample, 0, sizeof(csv_formatter_sample_t));
    
    char *newSampleName = malloc(strlen(sampleName) + 1);;
    strcpy(newSampleName, sampleName);
    newSample->sampleName = newSampleName;
    
    newSample->allele = allele;
    
    return newSample;
}


void csv_formatter_sample_destroy(csv_formatter_sample_t* sample)
{
    free((char *)sample->sampleName);
    free(sample);
}

csv_formatter_variation_list_t *csv_formatter_variation_list_init(int32_t sampleCount, int32_t position)
{
    csv_formatter_variation_list_t *newVariations = malloc(sizeof(csv_formatter_variation_list_t));
    
    newVariations->position = position;
    newVariations->variationsCount = sampleCount;
    newVariations->variations = malloc(sizeof(char *) * sampleCount);
    int32_t i;
    for (i = 0; i < newVariations->variationsCount; i++) {
        newVariations->variations[i] = emptyString;
    }
    
    return newVariations;
}

void csv_formatter_variation_list_destroy(csv_formatter_variation_list_t *variationList)
{
    int i;
    for (i = 0; i < variationList->variationsCount; i++) {
        if (variationList->variations[i] != emptyString) {
            free((char *)variationList->variations[i]);
        }
    }
    free(variationList->variations);
    free(variationList);
}

void csv_formatter_variation_list_add(csv_formatter_variation_list_t *variationList, const char * variation, int32_t sampleIndex)
{
    if (variationList->variationsCount <= sampleIndex) {
        fprintf(stderr, "[%s:%d %s] sampleIndex %d out of bounds\n", __FILE__, __LINE__, __FUNCTION__, sampleIndex);
        exit(1);
    }
    
    if (variationList->variations[sampleIndex] != emptyString) {
        free((char *)variationList->variations[sampleIndex]);
    }
    
    char *newVariation = malloc(strlen(variation) + 1);
    strcpy(newVariation, variation);
    variationList->variations[sampleIndex] = newVariation;
}

csv_formatter_t *csv_formatter_init(bcf_hdr_t *bcfHeader)
{
    csv_formatter_t *newFormatter = malloc(sizeof(csv_formatter_t));
    memset(newFormatter, 0, sizeof(csv_formatter_t));

    newFormatter->referenceSample = csv_formatter_sample_init("reference", 0);
    
    newFormatter->sampleCount = bcf_hdr_nsamples(bcfHeader) * 2;
    newFormatter->samples = malloc(sizeof(csv_formatter_sample_t*) * bcf_hdr_nsamples(bcfHeader) * 2);
    
    int i;
    for (i=0; i<bcf_hdr_nsamples(bcfHeader); i++)
    {
        char *name = bcfHeader->samples[i];
        
        newFormatter->samples[i*2] = csv_formatter_sample_init(name, 1);
        newFormatter->samples[i*2+1] = csv_formatter_sample_init(name, 2);
    }
    
    newFormatter->variationListsAllocated = 16;
    newFormatter->variationLists = malloc(sizeof(csv_formatter_variation_list_t *) * 16);
    memset(newFormatter->variationLists, 0, sizeof(csv_formatter_variation_list_t *) * 16);
    
    return newFormatter;
}


void csv_formatter_destroy(csv_formatter_t* csvFormatter)
{
    csv_formatter_sample_destroy(csvFormatter->referenceSample);
    
    int32_t i;
    for (i = 0; i < csvFormatter->sampleCount; i++) {
        csv_formatter_sample_destroy(csvFormatter->samples[i]);
    }
    free(csvFormatter->samples);
    
    for (i = 0; i < csvFormatter->variationListsCount; i++) {
        csv_formatter_variation_list_destroy(csvFormatter->variationLists[i]);
    }
    free(csvFormatter->variationLists);
    
    free(csvFormatter);
}


void csv_formatter_add_record(csv_formatter_t* csvFormatter, bcf_hdr_t *header, bcf1_t *record)
{
    if (bcf_is_snp(record) == 0) { // only handle SNPs for now
        fprintf(stderr, "***WARNING*** The CSV formater can only handle SNPs\n");
        return;
    }
    
    bcf_unpack(record, BCF_UN_ALL);
    
    int32_t *genemapPositionArray = NULL;
    int genemapPositionArrayLength = 0;
    int genemapPositionCount = 0;
    genemapPositionCount = bcf_get_info_int32(header, record, GENEMAP, &genemapPositionArray, &genemapPositionArrayLength);
    if (genemapPositionCount == -1) {
        fprintf(stderr, "***WARNING*** No gene mapping defined in the bcf header\n");
        free(genemapPositionArray);
      return;
    } else if (genemapPositionCount == -2) {
        fprintf(stderr, "***WARNING*** Wrong gene mapping type in the header\n");
        free(genemapPositionArray);
      return;
    } else if (genemapPositionCount == -3) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with no gene mappings\n");
        free(genemapPositionArray);
       return;
    } else if (genemapPositionCount < 0) {
        fprintf(stderr, "***WARNING*** Unknown error occured while reading the gene mappings\n");
        free(genemapPositionArray);
        return;
    }
    
    if (genemapPositionCount > 1) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with multiple gene mappings\n");
        free(genemapPositionArray);
        return;
    }
    if (genemapPositionCount < 1) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with 0 gene mappings\n");
        free(genemapPositionArray);
       return;
    }
    
    int32_t genemapPosition = genemapPositionArray[0];
    free(genemapPositionArray);
    genemapPositionArray = NULL;
   
    char *genemapStrandString = NULL;
    int genemapStrandStringLength = 0;
    int genemapStrandStringCount = 0;
    genemapStrandStringCount = bcf_get_info_string(header, record, GENEMAP_STRAND, &genemapStrandString, &genemapStrandStringLength);
    if (genemapStrandStringCount == -1) {
        fprintf(stderr, "***WARNING*** No gene mapping strand defined in the bcf header\n");
        free(genemapStrandString);
       return;
    } else if (genemapStrandStringCount == -2) {
        fprintf(stderr, "***WARNING*** Wrong gene mapping strand type in the header\n");
        free(genemapStrandString);
       return;
    } else if (genemapStrandStringCount == -3) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with no gene mapping strand\n");
        free(genemapStrandString);
        return;
    } else if (genemapStrandStringCount < 0) {
        fprintf(stderr, "***WARNING*** Unknown error occured while reading the gene mappings strand\n");
        free(genemapStrandString);
        return;
    }
    
    if (genemapStrandStringCount > 1) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with multiple gene mapping strands\n");
        free(genemapStrandString);
        return;
    }
    if (genemapStrandStringCount < 1) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with 0 gene mapping strands\n");
        free(genemapStrandString);
        return;
    }
    
    char genemapStrand = genemapStrandString[0];
    free(genemapStrandString);
    genemapStrandString = NULL;

    if (genemapStrand != plusstrand && genemapStrand != minusstrand) {
        fprintf(stderr, "***WARNING*** Trying to CSV format a variant call with an illegal gene mapping strand\n");
        return;
    }
    
    csv_formatter_variation_list_t *variationList = NULL;
    if (csvFormatter->variationListsCount > 0) {
        variationList = csvFormatter->variationLists[csvFormatter->variationListsCount - 1];
    }
    
    if (variationList == NULL || variationList->position != genemapPosition) {
        if (csvFormatter->variationListsCount == csvFormatter->variationListsAllocated) {
            csvFormatter->variationListsAllocated *= 2;
            csvFormatter->variationLists = realloc(csvFormatter->variationLists, csvFormatter->variationListsAllocated);
            memset(csvFormatter->variationLists + csvFormatter->variationListsCount, 0,
                   sizeof(csv_formatter_variation_list_t *) * (csvFormatter->variationListsAllocated - csvFormatter->variationListsCount));
        }
        csvFormatter->variationListsCount++;
        csvFormatter->variationLists[csvFormatter->variationListsCount - 1] = csv_formatter_variation_list_init(csvFormatter->sampleCount + 1, genemapPosition);
        variationList = csvFormatter->variationLists[csvFormatter->variationListsCount - 1];
    }
    
    int *genotypesArray = NULL;
    int genotypesCount = 0;
    int genotypesArrayLength = 0;
    
    genotypesCount = bcf_get_genotypes(header, record, &genotypesArray, &genotypesArrayLength);
    if (genotypesCount < 0) {
        fprintf(stderr, "Error getting genotypes\n");
        free(genotypesArray);
        exit(1);
    }
    if (genotypesCount != csvFormatter->sampleCount) {
        fprintf(stderr, "***WARNING*** Not diploid\n");
        free(genotypesArray);
        return;
    }
    
    // add reference variant
    char *referenceVariation = record->d.allele[0];
    char *referenceVariationComplement = NULL;
    
    if (genemapStrand == minusstrand) {
        referenceVariationComplement = malloc(strlen(referenceVariation) + 1);
        strcpy(referenceVariationComplement, complement_nucleotide_sequence(referenceVariation));
        referenceVariation = referenceVariationComplement;
    }
    
    csv_formatter_variation_list_add(variationList, referenceVariation, 0);
    free(referenceVariationComplement);
    
    int i;
    for (i = 0; i < csvFormatter->sampleCount / 2; i++) {
        int genotypeIndex1 = genotypesArray[i*2];
        int genotypeIndex2 = genotypesArray[(i*2)+1];
        const char *genotype1 = NULL;
        const char *genotype2 = NULL;
        if (genotypeIndex1 == bcf_gt_missing) {
            genotype1 = "N";
        } else if (genotypeIndex1 == bcf_int32_vector_end) {
            genotype1 = "-";
        } else {
            genotype1 = record->d.allele[bcf_gt_allele(genotypeIndex1)];
        }
        if (genotypeIndex2 == bcf_gt_missing) {
            genotype2 = "N";
        } else if (genotypeIndex2 == bcf_int32_vector_end) {
            genotype1 = "-";
        } else {
            genotype2 = record->d.allele[bcf_gt_allele(genotypeIndex2)];
        }
        
        char *genotype1Complement = NULL;
        char *genotype2Complement = NULL;
        if (genemapStrand == minusstrand) {
            printf("genotype1 = %s\n", genotype1);
            printf("genotype2 = %s\n", genotype2);
            genotype1Complement = malloc(strlen(genotype1) + 1);
            strcpy(genotype1Complement, complement_nucleotide_sequence(genotype1));
            genotype1 = genotype1Complement;
            genotype2Complement = malloc(strlen(genotype2) + 1);
            strcpy(genotype2Complement, complement_nucleotide_sequence(genotype2));
            genotype2 = genotype2Complement;
            printf("complement genotype1 = %s\n", genotype1);
            printf("complement genotype2 = %s\n", genotype2);
        }
        
        char *concat_genotype = NULL;
        if (bcf_gt_is_phased(genotypeIndex2) == 0) {
            concat_genotype = malloc(strlen(genotype1) + strlen(genotype2) + 5);
            sprintf(concat_genotype, "(%s, %s)", genotype1, genotype2);
            genotype1 = concat_genotype;
            genotype2 = concat_genotype;
        }
        
        csv_formatter_variation_list_add(variationList, genotype1, (i*2)+1); // +1 because of reference genome
        csv_formatter_variation_list_add(variationList, genotype2, (i*2)+2); // +1 because of reference genome
        
        free(concat_genotype);
        free(genotype1Complement);
        free(genotype2Complement);
    }
    
    free(genotypesArray);
}

void csv_formatter_print(csv_formatter_t* csvFormatter, FILE *fp)
{
    int i;
    int j;
    
    fprintf(fp, "Sample");
    for (i = 0; i< csvFormatter->variationListsCount; i++) {
        fprintf(fp, "\t%d", (int)csvFormatter->variationLists[i]->position);
    }
    fprintf(fp, "\n");

    fprintf(fp, "%s", csvFormatter->referenceSample->sampleName);
    for (i = 0; i < csvFormatter->variationListsCount; i++) {
        fprintf(fp, "\t%s", csvFormatter->variationLists[i]->variations[0]);
    }
    fprintf(fp, "\n");
    
    for (i = 0; i < csvFormatter->sampleCount; i++) {
        fprintf(fp, "%s (%d)", csvFormatter->samples[i]->sampleName, csvFormatter->samples[i]->allele);
        for (j = 0; j < csvFormatter->variationListsCount; j++) {
            fprintf(fp, "\t%s", csvFormatter->variationLists[j]->variations[i+1]);
        }
        fprintf(fp, "\n");
    }
}





