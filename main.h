//
//  main.h
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/15/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#ifndef bcfgenemapper_main_h
#define bcfgenemapper_main_h

#define GENEMAP "GENEMAP"
#define GENEMAP_INFO_HEADER "##INFO=<ID=" GENEMAP ",Number=1,Type=Integer,Description=\"Mapped Gene Location\">"

#define GENEMAP_STRAND GENEMAP "STRAND"
#define GENEMAP_STRAND_INFO_HEADER "##INFO=<ID=" GENEMAP_STRAND ",Number=1,Type=Strins,Description=\"Mapped Gene Strand\">"

enum _strand_t {
    plusstrand = '+',
    minusstrand = '-'
};
typedef char strand_t;


#endif
