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


/* Flag set by ‘--verbose’. */
//static int verbose_flag;

int main(int argc, char * const *argv)
{
//    int c;
//    
//    while (1)
//    {
//        static struct option long_options[] =
//        {
//            /* These options set a flag. */
//            {"verbose", no_argument,       &verbose_flag, 1},
//            /* These options don't set a flag.
//             We distinguish them by their indices. */
////            {"add",     no_argument,       0, 'a'},
////            {"append",  no_argument,       0, 'b'},
////            {"delete",  required_argument, 0, 'd'},
////            {"create",  required_argument, 0, 'c'},
//            {"file",    required_argument, 0, 'f'},
//            {0, 0, 0, 0}
//        };
//        /* getopt_long stores the option index here. */
//        int option_index = 0;
//        
//        c = getopt_long (argc, argv, "abc:d:f:",
//                         long_options, &option_index);
//        
//        /* Detect the end of the options. */
//        if (c == -1)
//            break;
//        
//        switch (c)
//        {
//            case 0:
//                /* If this option set a flag, do nothing else now. */
//                if (long_options[option_index].flag != 0)
//                    break;
//                printf ("option %s", long_options[option_index].name);
//                if (optarg)
//                    printf (" with arg %s", optarg);
//                printf ("\n");
//                break;
//                
//            case 'a':
//                puts ("option -a\n");
//                break;
//                
//            case 'b':
//                puts ("option -b\n");
//                break;
//                
//            case 'c':
//                printf ("option -c with value `%s'\n", optarg);
//                break;
//                
//            case 'd':
//                printf ("option -d with value `%s'\n", optarg);
//                break;
//                
//            case 'f':
//                printf ("option -f with value `%s'\n", optarg);
//                break;
//                
//            case '?':
//                /* getopt_long already printed an error message. */
//                break;
//                
//            default:
//                abort ();
//        }
//    }
//    
//    /* Instead of reporting ‘--verbose’
//     and ‘--brief’ as they are encountered,
//     we report the final status resulting from them. */
//    if (verbose_flag)
//        puts ("verbose flag is set");
//    
//    /* Print any remaining command line arguments (not options). */
//    if (optind < argc)
//    {
//        printf ("non-option ARGV-elements: ");
//        while (optind < argc)
//            printf ("%s ", argv[optind++]);
//        putchar ('\n');
//    }
//    
//
  
    gene_mapper_t *geneMapper = gene_mapper_init();
    gene_mapper_add_exon(geneMapper, exon_range(25420785, 25420638));
    gene_mapper_add_exon(geneMapper, exon_range(25408868, 25408682));
    gene_mapper_add_exon(geneMapper, exon_range(25402745, 25402595));
    gene_mapper_add_exon(geneMapper, exon_range(25392140, 25391993));
    gene_mapper_add_exon(geneMapper, exon_range(25390914, 25390748));
    
    gene_mapper_add_exon(geneMapper, exon_range(25389112, 25388975));
    gene_mapper_add_exon(geneMapper, exon_range(25385843, 25385710));
    gene_mapper_add_exon(geneMapper, exon_range(25375427, 25375348));
    gene_mapper_add_exon(geneMapper, exon_range(25370539, 25370466));
    gene_mapper_add_exon(geneMapper, exon_range(25362552, 25362526));
    
    // do some debug reverse mapping
    printf("reverse map\n");
    printf("0 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 0));
    printf("1253 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 1253));
    
    printf("675 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 675));
    printf("817 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 817));
    printf("871 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 871));
    printf("889 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 889));
    printf("907 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 907));
    printf("1056 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 1056));
    printf("1060 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 1060));
    printf("1226 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 1226));
    printf("1249 -> %d\n", (int)gene_mapper_reversemap_position(geneMapper, 1249));
    
    htsFile *htsInFile = hts_open(argv[1], "r");
    htsFile *vcfOutFile = hts_open(argv[2], "w");

    bcf_hdr_t *bcf_header = bcf_hdr_read(htsInFile);
    
    int headerTextLen;
    char *headerText = bcf_hdr_fmt_text(bcf_header, 0, &headerTextLen);
    printf("Header Text:\n%s", headerText);
    

    bcf1_t *bcf_record = bcf_init();
    
    
    
    bcf_hdr_t *hdr_out = bcf_hdr_dup(bcf_header);
    
    int error = bcf_hdr_append(hdr_out, "##INFO=<ID=GENEMAP,Number=1,Type=Integer,Description=\"Gene Location\">");
        if (error) {
            printf("bcf_hdr_append error %d", error);
        }

    
    bcf_hdr_write(vcfOutFile, hdr_out);
    
    while (bcf_read(htsInFile, bcf_header, bcf_record)>=0 )
    {
        int32_t geneLocation = gene_mapper_map_position(geneMapper, bcf_record->pos);
        if (geneLocation >= 0) {
            geneLocation++; // we 0 index all locations, but locations in the vcf file are 1 indexed;
            bcf_update_info_int32(hdr_out, bcf_record, "GENEMAP", &geneLocation, 1);
            bcf_write(vcfOutFile, hdr_out, bcf_record);
        } else {
            printf("position %d not found\n", (int)bcf_record->pos);
        }
    }
    
//    bcf_hdr_destroy(hdr);
//    bcf_hdr_destroy(hdr_out);
    hts_close(htsInFile);
    hts_close(vcfOutFile);

    
    
//    int error = 0;
//    while ( (error = bcf_read(htsFile, bcf_header, bcf_record)) >= 0) {
//    	error = bcf_unpack(bcf_record, BCF_UN_ALL);
//        
//        if (error) {
//            printf("bcf_unpack error %d", error);
//        }
//
//    }
    
    exit (0);
}

