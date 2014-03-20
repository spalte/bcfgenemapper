//
//  main.h
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/15/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#ifndef bcfgenemapper_main_h
#define bcfgenemapper_main_h

#include <htslib/vcf.h>

#define GENEMAP "GENEMAP"
#define GENEMAP_INFO_HEADER "##INFO=<ID=" GENEMAP ",Number=1,Type=Integer,Description=\"Mapped Gene Location\">"

#define GENEMAP_NAME GENEMAP "NAME"
#define GENEMAP_NAME_INFO_HEADER "##INFO=<ID=" GENEMAP_NAME ",Number=1,Type=String,Description=\"Mapped Gene Name\">"

#define GENEMAP_STRAND GENEMAP "STRAND"
#define GENEMAP_STRAND_INFO_HEADER "##INFO=<ID=" GENEMAP_STRAND ",Number=1,Type=String,Description=\"Mapped Gene Strand\">"

#define GENEMAP_FILE_VERSION 0.1f
#define GENEMAP_FILE_VERSION_STRING "0.1"
#define GENEMAP_VERSION_STRING "genemapperVersion"
#define GENEMAP_VERSION_HEADER "##" GENEMAP_VERSION_STRING "=" GENEMAP_FILE_VERSION_STRING

char complement_nucleotide(char n);
const char *complement_nucleotide_sequence(const char *sequence);

enum _strand_t {
    plusstrand = '+',
    minusstrand = '-'
};
typedef char strand_t;

static const char plusstrandString[] = {plusstrand, 0};
static const char minusstrandString[] = {minusstrand, 0};

// returns 0 on success
int bcf_update_genemapper_info(const bcf_hdr_t *hdr, bcf1_t *line, int32_t index, strand_t strand);

#endif
