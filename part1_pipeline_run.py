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
	stacks_executables = '/usr/lib/stacks/bin' # this is where ustacks lives, which should be all we need. other stacks scripts are in /opt/software/stacks-1.26/scripts
	BWA = 'bwa'
	samtoolsPath = 'samtools'
	bcftoolsPath = 'bcftools'
    
    #########################################################################
    ### INPUTS AND OUTPUTS                                                ###
    #########################################################################
    
    #### PART 1: INITIAL ASSEMBLY AND DBR FILTERING
    # user only needs to specify parent directory; the remaining directories should be automatically generated
    dataDir = '/home/kpierce/AdelaideTicks/AllData'
    parentDir = '/mnt/sda1/"
    
    # the remaining directories are automatically generated from the parent directory
    catInDir = dataDir
    catOutDir = parentDir + '/cat_parallel/'
    filterInDir = parentDir
    filterOutDir =  parentDir + '/qual_filtered/'
    dbrInDir = filterOutDir
    dbrOutDir = parentDir + '/DBR_dir/'
    demultiplexInDir = filterOutDir
    demultiplexOutDir = parentDir + '/demultiplexed/'
    trimInDir = demultiplexOutDir
    trimOutDir = parentDir + '/trimmed/'
    
    #stacksInDir = trimOutDir
    #stacksOutDir = parentDir + '/StacksOutput/' # stacks doesn't allow an output to be specified
    #pseudorefInDir = stacksOutDir
    #pseudorefOutDir = parentDir + '/pseudoreference.fastq'
    #BWAoutDir = parentDir + '/BWA/'
    #DBRfilteredseqs = parentDir + '/dbrFiltered/'
    
    #### PART 2: REASSEMBLING THE FILTERED SEQUENCES
    #re_demultiplexInDir = DBRfilteredseqs
    #re_demultiplexOutDir = parentDir + '/dbrFiltered_demultiplexed/'
    #re_BWAinDir = DBRfilteredseqs
    #re_BWAoutDir = parentDir + '/dbrFiltered_BWA2/'
    #finalBCFout = parentDir + '/dbrFiltered_pseudorefMapped_genotypes2.bcf'
    #finalVCFout = parentDir + '/dbrFiltered_pseudorefMapped_genotypes2.vcf'
    
    #########################################################################
    ### FUNCTION CALLS TO RUN THE PIPELINE                                ###
    #########################################################################
    
    # CONCATENATE READ 1 WITH REVERSE OF READ 2
    parallel_concatenate(in_dir = catInDir, regexR1='R1', regexR2='R2', out_dir = catOutDir)
    
    '''
    # QUALITY FILTER DATA
    out_name = '.qual_filtered' # gets appended to input file name
    q = 30
    p = 50
    read = '.assembled.fastq' # extension for pear-assembled reads
    parallel_FASTQ_quality_filter(in_dir = filterInDir, 
                                   out_dir = filterOutDir, 
                                   out_name = out_name, 
                                   q = q, 
                                   p = p, 
                                   qualityFilter = qualityFilter,
                                   read = read)
    
    # MAKE DBR DICTIONARIES FOR QUAL FILTERED PEAR DATA
    seq_type = 'pear'
    parallel_DBR_dict(in_dir = dbrInDir, 
                       seqType = seq_type,
                       save = dbrOutDir,
                       dbr_start = -10,
                       dbr_stop = -2)
    
    # DEMULTIPLEX
    out_prefix = '/demultiplexed_'
    iterative_Demultiplex2(in_dir = demultiplexInDir, 
                          barcode_dir = '/home/kpierce/Flinders_RADseq/Barcodes', 
                          out_dir = demultiplexOutDir,
                          regexLibrary = 'IDX\d{1}',
                          demultiplexPath = demultiplexPath,
                          startPoint = 'barcodes', 
                          out_prefix = out_prefix)
    ##### NEED TO DETERMINE THE FIRST AND LAST BASE TO RETAIN
    ##### ALIGN FROM THE END OF THE STRING AND KEEP SET NUMBER OF BASES
    ##### maybe something like last_base = -10, first_base = last base + 100 (or whatever)?
    
    # TRIM TO UNIFORM LENGTH
    suffix = '_trimmed.fq'
    first_base = 11
    last_base = 196
    Trim(in_dir = trimInDir, 
         out_dir = trimOutDir, 
         suffix = suffix, 
         first_base = first_base, 
         last_base = last_base)
         
    ##### AFTER THE DATA ARE TRIMMED, SORT OUT BARCODE FILES TO GUIDE ASSEMBLY #####     
    
    # RUN USTACKS SIMULTANEOUSLY ON ALL LIBRARIES
    denovo_Ustacks(in_dir = stacksInDir, 
                  denovo_path = denovo_path, 
                  stacks_executables = stacks_executables, 
                  out_dir = stacksOutDir, 
                  m = 10, 
                  n = 2, 
                  b = 1, 
                  D = '_initial_assembly')
    
    # RUN CSTACKS SIMULTANEOUSLY ON ALL LIBRARIES (same args as above)
    denovo_Cstacks(in_dir = stacksInDir, 
                  denovo_path = denovo_path, 
                  stacks_executables = stacks_executables, 
                  out_dir = stacksOutDir, 
                  m = 10, 
                  n = 2, 
                  b = 1, 
                  D = '_initial_assembly')
    
    # GENERATE THE PSEUDOREFERENCE GENOME
    GeneratePseudoref(in_dir = pseudorefInDir, 
                      out_file = pseudorefOutDir,  
                      BWA_path = BWA) # imported from integrated_denovo_pipeline.py
    
    # REFERENCE MAP QUALITY FILTERED/DEMULTIPLEXED MERGED READS TO THE PSEUDOREFERENCE
    parallel_refmap_BWA(in_dir = trimOutDir, # input demultiplexed, trimmed reads
               out_dir = BWAoutDir, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
        
    # FILTER OUT PCR DUPLICATES USING THE DBR SEQUENCES
    DBR_Filter(assembled_dir = BWAoutDir, # the SAM files for the data mapped to pseudoreference
               out_dir = DBRfilteredseqs, # the output file, full path, ending with .fasta
               n_expected = 2, # the number of differences to be tolerated
               barcode_dir = '/home/pierce/CSU_ChronicWasting/RevisedBarcodes', # the barcodes for individuals in the library referenced in dict_in
               dict_dir = dbrOutDir, # a single dictionary of DBRs (for one library only)
               sample_regex = '.*_(\d{1,3}T?)_.*',
               barcode_file=None, # if just a single library is being used, can directly pass the barcode file
               test_dict=True, # optionally print testing info to stdout for checking the dictionary construction
               phred_dict=phred_dict, # dictionary containing ASCII quality filter scores to help with tie breaks
               samMapLen=None)
    
    # REFERENCE MAP DBR FILTERED READS TO THE PSEUDOREFERENCE
    parallel_refmap_BWA(in_dir = re_BWAinDir, # input demultiplexed, trimmed reads
               out_dir = re_BWAoutDir, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
    
    # CALL THE GENOTYPES USING SAMTOOLS MPILEUP; CONVERT OUTPUT TO VCF FILE
    callGeno(sam_in = re_BWAoutDir, 
             pseudoref = pseudorefOutDir, 
             BCFout = finalBCFout, 
             VCFout = finalVCFout,
             samtoolsPath = samtoolsPath,
             bcftoolsPath = bcftoolsPath)
'''