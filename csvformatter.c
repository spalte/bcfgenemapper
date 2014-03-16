//
//  csvformater.c
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/15/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#include <stdio.h>

#include "csvformatter.h"

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
    if (variationList->variationsCount >= sampleIndex) {
        fprintf(stderr, "[%s:%d %s] sampleIndex %d out of bounds", __FILE__,__LINE__,__FUNCTION__,sampleIndex);
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
    
    newFormatter->sampleCount = bcf_hdr_nsamples(bcfHeader);
    newFormatter->samples = malloc(sizeof(csv_formatter_sample_t*) * bcf_hdr_nsamples(bcfHeader) * 2);
    
    int i;
    for (i=0; i<bcf_hdr_nsamples(bcfHeader); i++)
    {
        char *name = bcfHeader->samples[i];
        
        newFormatter->samples[i*2] = csv_formatter_sample_init(name, 0);
        newFormatter->samples[i*2+1] = csv_formatter_sample_init(name,1);
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









