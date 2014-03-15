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
    int variationsCount;
    const char **variations;
    char allele; // 0 or 1
} csv_formatter_variations_t;

typedef struct {
    csv_formatter_sample_t *referenceSample;
    
    int sampleCount;
    csv_formatter_sample_t *samples;
} csv_formatter_t;

csv_formatter_t *csv_formatter_sample_init(const char *sampleName, char allele);
void csv_formatter_sample_destroy(csv_formatter_t* csvFormatter);

csv_formatter_t *csv_formatter_init(bcf_hdr_t *bcfHeader);
void csv_formatter_destroy(csv_formatter_t* csvFormatter);

void csv_formatter_add_record(csv_formatter_t* csvFormatter, bcf1_t *record);
void csv_formatter_print(csv_formatter_t* csvFormatter, FILE *fp);

#endif
