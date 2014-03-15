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

#define GENEMAP_STRAND GENEMAP "STRAND"
#define GENEMAP_STRAND_INFO_HEADER "##INFO=<ID=" GENEMAP_STRAND ",Number=1,Type=String,Description=\"Mapped Gene Strand\">"

enum _strand_t {
    plusstrand = '+',
    minusstrand = '-'
};
typedef char strand_t;

// returns 0 on success
int bcf_update_genemapper_info(const bcf_hdr_t *hdr, bcf1_t *line, int32_t index, strand_t strand);

#endif
