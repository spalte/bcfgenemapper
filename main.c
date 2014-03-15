//
//  main.c
//  bcfgenemapper
//
//  Created by Joël Spaltenstein on 3/6/14.
//  Copyright (c) 2014 Spaltenstein Natural Image. All rights reserved.
//

#include <stdio.h>
#include <htslib/vcf.h>
#include <getopt.h>

#include "genemapper.h"
#include "version.h"
#include "main.h"

static int verbose_flag;
static int strip_flag;

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

// returns 0 on success
int bcf_update_genemapper_info(const bcf_hdr_t *hdr, bcf1_t *line, int32_t index, strand_t strand)
{
    int oneBasedIndex = index+1; // we use 1 base indexing, but vcf uses 0 base indexing
    int error = 0;
    error = bcf_update_info_int32(hdr, line, GENEMAP, &oneBasedIndex, 1);
    if (error) {
        return error;
    }
    char strandStr[] = {0,0};
    strandStr[0] = strand;
    error = bcf_update_info_string(hdr, line, GENEMAP_STRAND, strandStr);
    if (error) {
        return error;
    }
    return 0;
}


void print_usage (FILE* stream, int exit_code)
{
    fprintf(stream, "Gene Mapper (%s, htslib version:%s)\n", BCFGENEMAPPER_VERSION, hts_version());
    fprintf(stream, "Usage:  %s <-e exon_filename> [options] [input_filename]\n", program_name);
    fprintf(stream,
            "  -h  --help                 Display this usage information.\n"
            "  -o  --output filename      Write output to file.\n"
            "  -O  --output-type b|u|z|v  Compressed BCF (b), Uncompressed BCF (u),\n"
            "                             Compressed VCF (z), Uncompressed VCF (v).\n"
            "  -e  --exons filename       Read exon ranges from this file.\n"
            "  -c  --csv filename         Write variants to a cvs file.\n"
            "  -v  --verbose              Print verbose messages.\n\n"
            "  -s  --strip                Don't output variants that are not in exons.\n"
            "  -v  --verbose              Print verbose messages.\n\n"
            "If the input or output filenames are omitted %s will use\n"
            "stdin and stdout respectively\n\n", program_name);
    fprintf(stream,
            "The format of the exon range file is: \"(start_position) (end_position)newLine\"\n"
            "Both start and end are inclusive and 0-indexed (like in NCBI XML files).\n"
            "If the exon is on the (-)strand, the start should be a larger index than\n"
            "the end index.\n\n"
            "Example for the RHCE reference peptide (NP_065231.3) onto GRCh38 Chr1\n"
            "(which happens to be on the (-)strand):\n"
            "25420785 25420638\n"
            "25408868 25408682\n"
            "25402745 25402595\n"
            "25392140 25391993\n"
            "25390914 25390748\n"
            "25389112 25388975\n"
            "25385843 25385710\n"
            "25375427 25375348\n"
            "25370539 25370466\n"
            "25362552 25362526\n");

    exit(exit_code);
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
        static const char* const short_options = "vsho:O:e:c:";
        static struct option long_options[] =
        {
            {"verbose",     no_argument,       NULL, 'v'},
            {"strip",       no_argument,       NULL, 's'},
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
            case 'v':
                verbose_flag = 1;
                break;
            case 's':
                strip_flag = 1;
                break;
            case 'o':
                output_filename = optarg;
                break;
            case 'O':
                output_type = optarg;
                if (validate_output_type(output_type) == 0) {
                    fprintf(stderr, "Invalid output type: \"%s\", legal values are b|u|z|v\n", output_type);
                    print_usage(stderr, 1);
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
                fprintf(stderr, "about to abort because c = %d, %c\n", c, (char)c);
                abort();
        }
    }
    
    // read the input file if it is there
    if (optind + 1 == argc) {
        input_filename = argv[optind];
    } else if (optind < argc) {
        fprintf(stderr, "Specify only one input file.\n");
        print_usage(stderr, 1);
    }
    
    if (input_filename == NULL && output_filename == NULL && exons_filename == NULL) {
        print_usage(stdout, 1);
    }
    
    if (exons_filename == NULL) {
        fprintf(stderr, "No exon file specified. Use -e to specify the exon file.\n");
        print_usage(stderr, 1);
    }
    
    FILE *exonFp = fopen(exons_filename, "r");
    if (exonFp < 0) {
        fprintf(stderr, "Unable to open exon file. \"%s\".", exons_filename);
        print_usage(stderr, 1);
    }
    
    gene_mapper_t *geneMapper = gene_mapper_file_init(exonFp);
    if (gene_mapper_exon_count(geneMapper) == 0) {
        fprintf(stderr, "Unable to read exons from file \"%s\".", exons_filename);
        print_usage(stderr, 1);
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
        fprintf(stderr, "Unable to open input file \"%s\".", input_filename);
        print_usage(stderr, 1);
    }
    
    char outputFileMode[3] = {'w', 'v', 0};
    if (output_type) {
        outputFileMode[1] = output_type[0];
    }
    htsFile *vcfOutFile = hts_open(output_filename, outputFileMode);
    if (htsInFile < 0) {
        fprintf(stderr, "Unable to open output file \"%s\".", output_filename);
        print_usage(stderr, 1);
    }
    
    if (verbose_flag) {
        gene_mapper_print_exons(geneMapper, stdout);
    }

    bcf_hdr_t *bcf_header = bcf_hdr_read(htsInFile);
    bcf_hdr_t *hdr_out = bcf_hdr_dup(bcf_header);
    
    int error = bcf_hdr_append(hdr_out, GENEMAP_INFO_HEADER);
    if (error) {
        fprintf(stderr, "bcf_hdr_append error %d", error);
        abort();
    }
    error = bcf_hdr_append(hdr_out, GENEMAP_STRAND_INFO_HEADER);
    if (error) {
        fprintf(stderr, "bcf_hdr_append error %d", error);
        abort();
    }

    bcf_hdr_write(vcfOutFile, hdr_out);
    
    int32_t keptRecords = 0;
    int32_t updatedRecords = 0;
    int32_t removedRecords = 0;
    
    bcf1_t *bcf_record = bcf_init();
    while (bcf_read(htsInFile, bcf_header, bcf_record)>=0 )
    {
        exon_range_t exon;
        int32_t geneLocation = gene_mapper_map_position(geneMapper, bcf_record->pos, &exon);
        
        if (geneLocation >= 0) {
            bcf_update_genemapper_info(hdr_out, bcf_record, geneLocation, exon_range_strand(exon));
            updatedRecords++;
            bcf_write(vcfOutFile, hdr_out, bcf_record);
            keptRecords++;
        } else if (strip_flag == 0) {
            bcf_write(vcfOutFile, hdr_out, bcf_record);
            keptRecords++;
        } else {
            removedRecords++;
        }
    }
    
    if (verbose_flag) {
        printf("%d record%s kept.\n", (int)keptRecords, keptRecords != 1?"s":"");
        printf("%d record%s updated.\n", (int)updatedRecords, updatedRecords != 1?"s":"");
        printf("%d record%s removed.\n", (int)removedRecords, removedRecords != 1?"s":"");
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

