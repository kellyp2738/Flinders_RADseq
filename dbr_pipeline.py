#!/usr/bin/python

#########################################################################
### De novo pipeline for processing complete Adelaide tick dataset    ###
#########################################################################

#### IMPORTS ####
from integrated_denovo_pipeline import *
#from DBR_Parsing import *
from assembled_DBR_filtering import *

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-t', '--threads', help = 'Number of threads to use for parallel operations.')
    opts = parser.parse_args()
    
    threads=int(opts.threads)
    
    #########################################################################
    ### SET UP THREADS TO USE FOR PARALLEL PROCESSES                      ###
    #########################################################################
    
    for i in xrange(threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()
    
    #########################################################################
    ### PATHS TO EXECUTABLES                                              ###
    #########################################################################
    
    #pearPath = 'pear-0.9.6-bin-64'
    qualityFilter = 'fastq_quality_filter'
    trimmer = 'fastx_trimmer'
    demultiplexPath = 'fastx_barcode_splitter.pl'
    denovo_path = 'denovo_map.pl'
    stacks_executables = '/usr/lib/stacks-1.42/bin' # this is where ustacks and cstacks live, which should be all we need. other stacks scripts are in /opt/software/stacks-1.26/scripts
    BWA = 'bwa'
    samtoolsPath = 'samtools'
    bcftoolsPath = 'bcftools'
    
    #########################################################################
    ### INPUTS AND OUTPUTS                                                ###
    #########################################################################
    
    #### PART 1: INITIAL ASSEMBLY AND DBR FILTERING
    # user only needs to specify parent directory; the remaining directories should be automatically generated
    #dataDir = '/home/kpierce/AdelaideTicks/AllData'
    parentDir = '/mnt/HGST4TB/legs_larvae_reference_assembly_output' # THIS CHANGED PART-WAY THROUGH THE ANALYSIS AS DATA WERE MOVED TO A BIGGER DRIVE
    
    # the remaining directories are automatically generated from the parent directory
    #catInDir = dataDir
    #catOutDir = parentDir + '/cat_parallel/'
    
    #filterInDir = catOutDir
    #filterOutDir =  parentDir + '/qual_filtered/'
    
    #mergeLanesInDir = filterOutDir
    #mergeLanesOutDir = parentDir + '/qual_filtered_merged_lanes/'
    
    # demultiplexing done with GBSX using a shell script (no python wrapper yet)
    #demultiplexInDir = mergeLanesOutDir
    #demultiplexOutDir = parentDir + '/demultiplexed/'
    
    # bodies
    #dbrInDir = '/mnt/HGST4TB/demultiplexed_testing_data'
    dbrInDir_bodies = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Bodies' # b. hydro
    dbrInDir_others = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Bodies_nonBH' # non b. hydro
    
    # legs
    dbrInDir_legs = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Legs'
    
    # larvae
    dbrInDir_larvae = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Larvae'
    
    # since samples will be duplicated in the bodies/legs groups and file name won't specify the tissue type, we need different directories for tissue type to ensure the right DBR dict is referenced
    dbrOutDir_legs_larvae = parentDir + '/DBR_dirs_legs_larvae/'
    dbrOutDir_bodies = parentDir + '/DBR_dirs_bodies/'
    
    trimInDir_bodies = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Bodies'
    trimInDir_others = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Bodies_nonBH'
    trimInDir_legs = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Legs'
    trimInDir_larvae = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Larvae'
    
    trimOutDir_BWA = parentDir + '/trimmed_for_BWA' # bodies
    trimOutDir_stacks = parentDir + '/trimmed_for_Stacks/' # legs & larvae
    
    stacksInDir = trimOutDir_stacks
    stacksOutDir = parentDir + '/Stacks_Legs_Larvae_Output/' # stacks doesn't allow an output to be specified
    
    pseudorefInDir = stacksOutDir
    pseudorefOutDir = parentDir + '/Legs_Larvae_pseudoreference.fastq'
    
    BWAoutDir_Legs_Larvae = parentDir + '/BWA_Legs_Larvae/'
    BWAoutDir_Bodies = parentDir + '/BWA_Bodies/'
    
    DBRfilteredseqs_Legs_Larvae = parentDir + '/dbrFiltered_Legs_Larvae/'
    DBRfilteredseqs_Bodies = parentDir + '/dbrFiltered_Bodies/'
    
    #### PART 2: REASSEMBLING THE FILTERED SEQUENCES
    #re_demultiplexInDir = DBRfilteredseqs
    #re_demultiplexOutDir = parentDir + '/dbrFiltered_demultiplexed/'
    #re_BWAinDir = DBRfilteredseqs
    #re_BWAoutDir = parentDir + '/dbrFiltered_BWA2/'
    #finalBCFout = parentDir + '/testing_bodies_pseudorefMapped_genotypes.bcf'
    #finalVCFout = parentDir + '/testing_bodies_pseudorefMapped_genotypes.vcf'
    
    #########################################################################
    ### FUNCTION CALLS TO RUN THE PIPELINE                                ###
    #########################################################################
    
    '''
    # COMPLETED 10/2/2016
    # CONCATENATE READ 1 WITH REVERSE OF READ 2
    parallel_concatenate(in_dir = catInDir, regexR1='R1', regexR2='R2', out_dir = catOutDir)

    # COMPLETED 10/10/2016 (-ish)
    # QUALITY FILTER DATA
    out_name = '.qual_filtered.gz' # gets appended to input file name
    q = 30
    p = 50
    parallel_FASTQ_quality_filter(in_dir = filterInDir, 
                                   out_dir = filterOutDir, 
                                   out_name = out_name, 
                                   q = q, 
                                   p = p, 
                                   qualityFilter = qualityFilter,
                                   read = read
    
    ## this isn't working, but was done manually with zcat
    # zcat IDXn* > IDXn_qual_filtered_fully_concatenated.fq.gz
    # saved in mergeLanesOutDir
    # COMPLETED 10/14/2016
    # PASTE LIBRARIES SPLIT ACROSS LANES INTO A SINGLE FILE
    parallel_merge_lanes(in_dir = mergeLanesInDir,
                         regexLibrary = 'IDX[0-9]',
                         out_dir = mergeLanesOutDir)
    
    #COMPLETED FOR INDEX 1 ON 10/23/2016 -- TOOK APPROXIMATELY 11 HOURS TO PROCESS THAT 18GB FILE
                       
    # A BASH SCRIPT WAS USED TO RUN A BUNCH OF CALLS TO GBSX DEMULTIPLEXER INSTEAD OF USING THIS FASTX BARCODE SPLITTER WRAPPER
    # DEMULTIPLEX
    out_prefix = '/demultiplexed_'
    iterative_Demultiplex2(in_dir = demultiplexInDir, 
                          barcode_dir = '/mnt/HGST4TB/Flinders_RADseq/Barcodes_First5', 
                          out_dir = demultiplexOutDir,
                          regexLibrary = 'IDX\d{1}',
                          demultiplexPath = demultiplexPath,
                          startPoint = 'barcodes')
    
    # COMPLETED FOR 65 SAMPLE TEST SET ON 6/23/2017 (10 MIN RUN TIME)
    # MAKE DBR DICTIONARIES FOR QUAL FILTERED PEAR DATA
    seq_type = 'pear' # even though these aren't pear-merged, this is still the proper argument -- it tells the fxn that the reads are merged/concatenated
    parallel_DBR_dict(in_dir = dbrInDir_bodies, 
                      seqType = seq_type, 
                      dbr_start = -10, 
                      dbr_stop = None, 
                      threads = threads,
                      test_dict = True, 
                      save = dbrOutDir_bodies)
                      
    parallel_DBR_dict(in_dir = dbrInDir_others, 
                      seqType = seq_type, 
                      dbr_start = -10, 
                      dbr_stop = None, 
                      threads = threads,
                      test_dict = True, 
                      save = dbrOutDir_bodies)
                      
    parallel_DBR_dict(in_dir = dbrInDir_legs, 
                      seqType = seq_type, 
                      dbr_start = -10, 
                      dbr_stop = None, 
                      threads = threads,
                      test_dict = True, 
                      save = dbrOutDir_legs_larvae)
                      
    parallel_DBR_dict(in_dir = dbrInDir_larvae, 
                      seqType = seq_type, 
                      dbr_start = -10, 
                      dbr_stop = None, 
                      threads = threads,
                      test_dict = True, 
                      save = dbrOutDir_legs_larvae)
    '''                      
    ## DO NOT REDO THIS STEP FOR LEGS AND BODIES
    ## ONLY TRIM LARVAL SEQUENCES (NOT YET DONE)
    # COMPLETED IN THE WEE HOURS OF 10/25/2016
    # TRIM TO UNIFORM LENGTH
    suffix = '_trimmed.fq'
    first_base = 1 # barcodes and enzyme cut sites already trimmed using GBSX demultiplexer; start at 1 for fastx trimmer
    # this is the number of bases left after the longest barcode + other sequences are removed: 
    #       150bp read length - 9mer barcode (max) - 4mer R1 cut - 4mer R2 cut - 10mer DBR = 123
    last_base = 123 
    #parallel_Trim(in_dir = trimInDir_bodies, 
    #     out_dir = trimOutDir_BWA, 
    #     trimPath = trimmer, 
    #     threads = threads,
    #     first_base = first_base, 
    #     last_base = last_base,
    #     suffix = suffix)
         
    #parallel_Trim(in_dir = trimInDir_others, 
    #     out_dir = trimOutDir_BWA, 
    #     trimPath = trimmer, 
    #     threads = threads,
    #     first_base = first_base, 
    #     last_base = last_base,
    #     suffix = suffix)
         
    #parallel_Trim(in_dir = trimInDir_legs, 
    #     out_dir = trimOutDir_stacks, 
    #     trimPath = trimmer, 
    #     threads = threads,
    #     first_base = first_base, 
    #     last_base = last_base,
    #     suffix = suffix)
         
    #parallel_Trim(in_dir = trimInDir_larvae, 
    #     out_dir = trimOutDir_stacks, 
    #     trimPath = trimmer, 
    #     threads = threads,
    #     first_base = first_base, 
    #     last_base = last_base,
    #     suffix = suffix)
         
    # COMPLETED 10/25/2016, BUT PROGRAM DIDN'T EXIT OR PROCEED TO denovo_Cstacs()
    # I think the function needs a return?
    # RUN USTACKS SIMULTANEOUSLY ON ALL LIBRARIES
    denovo_Ustacks(in_dir = stacksInDir, 
                  denovo_path = denovo_path, 
                  stacks_executables = stacks_executables, 
                  out_dir = stacksOutDir, 
                  m = 10, 
                  n = 2, 
                  num_threads = threads, # -p
                  b = 1, 
                  D = '_initial_assembly',
                  unmatchedName = 'undetermined')

    # COMPLETED 10/26/2016
    # RUN CSTACKS SIMULTANEOUSLY ON ALL LIBRARIES (same args as above)
    denovo_Cstacks(in_dir = stacksInDir, 
                  denovo_path = denovo_path, 
                  stacks_executables = stacks_executables, 
                  out_dir = stacksOutDir, 
                  m = 10, 
                  n = 2,
                  num_threads = threads, # -p
                  b = 1, 
                  D = '_initial_assembly')
                  
    # COMPLETED 10/26/2016          
    # GENERATE THE PSEUDOREFERENCE GENOME
    GeneratePseudoref(in_dir = pseudorefInDir, 
                      out_file = pseudorefOutDir,  
                      BWA_path = BWA) # imported from integrated_denovo_pipeline.py
                      
    # COMPLETED 10/26/2016
    # REFERENCE MAP QUALITY FILTERED/DEMULTIPLEXED MERGED READS TO THE PSEUDOREFERENCE
    parallel_refmap_BWA(in_dir = trimOutDir_BWA, # input demultiplexed, trimmed reads
               out_dir = BWAoutDir_Bodies, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
               
    parallel_refmap_BWA(in_dir = trimOutDir_stacks, # input demultiplexed, trimmed reads
               out_dir = BWAoutDir_Legs_Larvae, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
               
    # FILTER OUT PCR DUPLICATES USING THE DBR SEQUENCES
    parallel_DBR_Filter(assembled_dir = BWAoutDir_Bodies, # the SAM files for the data mapped to pseudoreference
               out_dir = DBRfilteredseqs_Bodies, # the output file, full path, ending with .fasta
               n_expected = 2, # the number of differences to be tolerated
               barcode_dir = None, # the barcodes for individuals in the library referenced in dict_in
               dict_dir = dbrOutDir_bodies, # a single dictionary of DBRs (for one library only)
               sample_regex = r'(.*_psti.R1)',
               sam_list = sample_sam_file_list,
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None)
               
    parallel_DBR_Filter(assembled_dir = BWAoutDir_Legs_Larvae, # the SAM files for the data mapped to pseudoreference
               out_dir = DBRfilteredseqs_Legs_Larvae, # the output file, full path, ending with .fasta
               n_expected = 2, # the number of differences to be tolerated
               barcode_dir = None, # the barcodes for individuals in the library referenced in dict_in
               dict_dir = dbrOutDir_legs_larvae, # a single dictionary of DBRs (for one library only)
               sample_regex = r'(.*_psti.R1)',
               sam_list = sample_sam_file_list,
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None)
    
    # COMPLETED 10/26/2016 
    # REFERENCE MAP DBR FILTERED READS TO THE PSEUDOREFERENCE
    #parallel_refmap_BWA(in_dir = re_BWAinDir, # input demultiplexed, trimmed reads
    #           out_dir = re_BWAoutDir, 
    #           BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
    #           pseudoref_full_path = pseudorefOutDir)
               
    # CALL THE GENOTYPES USING SAMTOOLS MPILEUP; CONVERT OUTPUT TO VCF FILE
    #callGeno(sam_in = BWAoutDir, 
    #         pseudoref = pseudorefOutDir, 
    #         BCFout = finalBCFout, 
    #         VCFout = finalVCFout,
    #         samtoolsPath = samtoolsPath,
    #         bcftoolsPath = bcftoolsPath)

