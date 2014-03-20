//
//  csvformater.h
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/15/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#ifndef bcfgenemapper_csvformater_h
#define bcfgenemapper_csvformater_h

#include <stdio.h>
#include <htslib/vcf.h>

typedef struct {
    const char *sampleName;
    char allele; // 0 or 1
} csv_formatter_sample_t;

typedef struct {
    int32_t position;
    char *sequenceName;;
    
    int32_t variationsCount; // this will be the number of samples
    const char **variations;
} csv_formatter_variation_list_t;

typedef struct {
    csv_formatter_sample_t *referenceSample;
    
    int32_t sampleCount;
    csv_formatter_sample_t **samples;
    
    int32_t variationListsCount;
    int32_t variationListsAllocated;
    csv_formatter_variation_list_t **variationLists;
} csv_formatter_t;

csv_formatter_sample_t *csv_formatter_sample_init(const char *sampleName, char allele);
void csv_formatter_sample_destroy(csv_formatter_sample_t* sample);

csv_formatter_variation_list_t *csv_formatter_variation_list_init(int32_t sampleCount, int32_t position);
void csv_formatter_variation_list_destroy(csv_formatter_variation_list_t *variationList);
void csv_formatter_variation_list_add(csv_formatter_variation_list_t *variationList, const char * variation, int32_t sampleIndex);

csv_formatter_t *csv_formatter_init(bcf_hdr_t *bcfHeader);
void csv_formatter_destroy(csv_formatter_t* csvFormatter);

void csv_formatter_collapse_variant_lists(csv_formatter_t* csvFormatter);
void csv_formatter_sort_variant_lists(csv_formatter_t* csvFormatter);

void csv_formatter_add_record(csv_formatter_t* csvFormatter, bcf_hdr_t *header, bcf1_t *record);
void csv_formatter_add_postition(csv_formatter_t* csvFormatter, int32_t position, const char *referenceSequence);
void csv_formatter_print(csv_formatter_t* csvFormatter, FILE *fp);

#endif
