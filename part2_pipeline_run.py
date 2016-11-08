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
	stacks_executables = '/usr/lib/stacks/bin' # this is where ustacks and cstacks live, which should be all we need. other stacks scripts are in /opt/software/stacks-1.26/scripts
	BWA = 'bwa'
	samtoolsPath = 'samtools'
	bcftoolsPath = 'bcftools'
    
    #########################################################################
    ### INPUTS AND OUTPUTS                                                ###
    #########################################################################
    
    #### PART 1: INITIAL ASSEMBLY AND DBR FILTERING
    # user only needs to specify parent directory; the remaining directories should be automatically generated
    dataDir = '/home/kpierce/AdelaideTicks/AllData'
    parentDir = '/mnt/HGST4TB/' # THIS CHANGED PART-WAY THROUGH THE ANALYSIS AS DATA WERE MOVED TO A BIGGER DRIVE
    
    # the remaining directories are automatically generated from the parent directory
    # demultiplexing done with GBSX using a shell script (no python wrapper yet)
    trimInDir = '/mnt/HGST4TB/demultiplexed_gsbx/demultiplexed_Bodies'
    trimOutDir = parentDir + '/trimmed_legs_and_bodies/'
    
    pseudorefOutDir = parentDir + '/Legs_pseudoreference.fastq'
    BWAoutDir = parentDir + '/BWA_Legs_and_Bodies/'

    finalBCFout = parentDir + '/Legs_pseudorefMapped_genotypes.bcf'
    finalVCFout = parentDir + '/Legs_pseudorefMapped_genotypes.vcf'
    
    #########################################################################
    ### FUNCTION CALLS TO RUN THE PIPELINE                                ###
    #########################################################################


    # TRIM TO UNIFORM LENGTH
    suffix = '_trimmed.fq'
    first_base = 1 # barcodes and enzyme cut sites already trimmed using GBSX demultiplexer
    last_base = 123 # this is the number of bases left after the longest barcode + other sequences are removed: 150bp read length - 9mer barcode (max) - 4mer R1 cut - 4mer R2 cut - 10mer DBR = 123
    parallel_Trim(in_dir = trimInDir, 
         out_dir = trimOutDir, 
         trimPath = trimmer, 
         first_base = first_base, 
         last_base = last_base,
         suffix = suffix)
    
    # REFERENCE MAP QUALITY FILTERED/DEMULTIPLEXED MERGED READS TO THE PSEUDOREFERENCE
    parallel_refmap_BWA(in_dir = trimOutDir, # input demultiplexed, trimmed reads
               out_dir = BWAoutDir, 
               BWA_path = BWA, # imported from integrated_denovo_pipeline.py 
               pseudoref_full_path = pseudorefOutDir)
               
    # CALL THE GENOTYPES USING SAMTOOLS MPILEUP; CONVERT OUTPUT TO VCF FILE
    callGeno(sam_in = BWAoutDir, 
             pseudoref = pseudorefOutDir, 
             BCFout = finalBCFout, 
             VCFout = finalVCFout,
             samtoolsPath = samtoolsPath,
             bcftoolsPath = bcftoolsPath,
             threads = threads)