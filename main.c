//
//  main.c
//  vcfparser
//
//  Created by Joël Spaltenstein on 3/6/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#include <stdio.h>
#include <htslib/vcf.h>
#include <getopt.h>

#include "genemapper.h"
#include "version.h"


/* Flag set by ‘--verbose’. */
static int verbose_flag;

static const char* program_name;

char validate_output_type(const char *type)
{
    if (strcmp(type, "b") == 0) {
        return 1;
    } else if (strcmp(type, "u") == 0) {
        return 1;
    } else if (strcmp(type, "z") == 0) {
        return 1;
    } else if (strcmp(type, "v") == 0) {
        return 1;
    } else {
        return 0;
    }
}

void print_usage (FILE* stream, int exit_code)
{
    fprintf (stream, "Gene Mapper (%s)\n", BCFGENEMAPPER_VERSION);
    fprintf (stream, "Usage:  %s <-e exon_filename> [options] [input filename]\n", program_name);
    fprintf (stream,
             "  -h  --help                 Display this usage information.\n"
             "  -o  --output filename      Write output to file.\n"
             "  -O  --output-type b|u|z|v  Compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v).\n"
             "  -e  --exons filename       Read exon ranges from this file.\n"
             "  -c  --csv filename         Write variants to this file.\n"
             "  -v  --verbose              Print verbose messages.\n\n"
             "If the input or output filenames are ommited %s will use stdin and stdout respectively\n", program_name);
    exit (exit_code);
}

int main(int argc, char * const *argv)
{
    int c;
    
    const char* input_filename = NULL;
    const char* output_filename = NULL;
    const char *output_type = NULL;
    const char *exons_filename = NULL;
    const char *csv_filename = NULL;
    
    program_name = argv[0];
    verbose_flag = 0;
    
    while (1)
    {
        static const char* const short_options = "vho:O:e:c:";
        static struct option long_options[] =
        {
            {"verbose",     no_argument,       NULL, 'v'},
            {"help",        no_argument,       NULL, 'h'},

            {"output",      required_argument, NULL, 'o'},
            {"output-type", required_argument, NULL, 'O'},
            {"exons",       required_argument, NULL, 'e'},
            {"csv",         required_argument, NULL, 'c'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        
        c = getopt_long(argc, argv, short_options, long_options, &option_index);
        
        if (c == -1)
            break;
        
        switch (c)
        {
            case 'h':
                print_usage (stdout, 0);
            
            case 'v':   /* -v or --verbose */
                verbose_flag = 1;
                break;

            case 'o':
                output_filename = optarg;
                break;
                
            case 'O':
                output_type = optarg;
                if (validate_output_type(output_type) == 0) {
                    fprintf(stdout, "Invalid output type\n");
                    print_usage(stdout, 1);
                }
                break;
                
            case 'e':
                exons_filename = optarg;
                break;
                

            case 'c':
                csv_filename = optarg;
                break;
                
            case '?':
                print_usage(stdout, 1);
                break;
                
            default:
                abort();
        }
    }
    
    // read the input file if it is there
    if (optind + 1 == argc) {
        input_filename = argv[optind];
    } else if (optind < argc) {
        printf("Specify only one input file.\n");
        print_usage(stdout, 1);
    }
    
    if (input_filename == NULL && output_filename == NULL && exons_filename == NULL) {
        print_usage(stdout, 1);
    }
    
    if (exons_filename == NULL) {
        printf("No exon file specified.\n");
        print_usage(stdout, 1);
    }
    
    FILE *exonFp = fopen(exons_filename, "r");
    if (exonFp < 0) {
        printf("Unable to open exon file. \"%s\".", exons_filename);
        print_usage(stdout, 1);
    }
    
    gene_mapper_t *geneMapper = gene_mapper_file_init(exonFp);
    if (gene_mapper_exon_count(geneMapper) == 0) {
        printf("Unable to read exons from file \"%s\".", exons_filename);
        print_usage(stdout, 1);
    }
    fclose(exonFp);
    exonFp = NULL;

    
    if (input_filename == NULL) {
        input_filename = "-";
    }
    if (output_filename == NULL) {
        output_filename = "-";
    }
    
    htsFile *htsInFile = hts_open(input_filename, "r");
    if (htsInFile < 0) {
        printf("Unable to open input file \"%s\".", input_filename);
        print_usage(stdout, 1);
    }
    
    char outputFileMode[3] = {'w', 'b', 0};
    if (output_type) {
        outputFileMode[1] = output_type[0];
    }
    htsFile *vcfOutFile = hts_open(output_filename, outputFileMode);
    if (htsInFile < 0) {
        printf("Unable to open output file \"%s\".", output_filename);
        print_usage(stdout, 1);
    }

    bcf_hdr_t *bcf_header = bcf_hdr_read(htsInFile);
    bcf_hdr_t *hdr_out = bcf_hdr_dup(bcf_header);
    
    int error = bcf_hdr_append(hdr_out, "##INFO=<ID=GENEMAP,Number=1,Type=Integer,Description=\"Gene Location\">");
    if (error) {
        printf("bcf_hdr_append error %d", error);
        abort();
    }

    bcf_hdr_write(vcfOutFile, hdr_out);
    
    bcf1_t *bcf_record = bcf_init();
    while (bcf_read(htsInFile, bcf_header, bcf_record)>=0 )
    {
        int32_t geneLocation = gene_mapper_map_position(geneMapper, bcf_record->pos);
        if (geneLocation >= 0) {
            geneLocation++; // we 0 index all locations, but locations in the vcf file are 1 indexed;
            bcf_update_info_int32(hdr_out, bcf_record, "GENEMAP", &geneLocation, 1);
            bcf_write(vcfOutFile, hdr_out, bcf_record);
        }
    }
    
    hts_close(htsInFile);
    htsInFile = NULL;
    hts_close(vcfOutFile);
    vcfOutFile = NULL;

    gene_mapper_destroy(geneMapper);
    geneMapper = NULL;
    
    bcf_hdr_destroy(bcf_header);
    bcf_header = NULL;
    bcf_hdr_destroy(hdr_out);
    hdr_out = NULL;
    bcf_destroy(bcf_record);
    bcf_record = NULL;

    exit (0);
}

