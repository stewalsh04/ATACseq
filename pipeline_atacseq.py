"""=====================================

******This is work in progress******

Code to run a downstream analysis of ATACseq data

The code will first run ATACqc, then count over a master peak file, then run DESeq2 to carry out further QC and differential analysis

I plan to add an option to normalise the reads to a housekeeping gene list using RUVseq before carry out the DEseq analysis.

The pipline will take:

Bam files
A master peaks file in GTF file format to count over (I also plan to integrate a convertor here, so a .bed file can be used instead)
A sample sheet for DESeq2


=============================================================

"""
import sys
import os
from ruffus import *
from CGATCore import Pipeline as P


# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])
    
#Perform ATACqc on .bam files
@transform('*.bam',
           suffix(".bam"),
           "_atacqc.html")
def atac_qc(infile, outfile):
    '''Perform atac qc on bam files.'''
    
    cwd = os.getcwd()
    
    statement = '''R -e "rmarkdown::render('/ifs/projects/proj078/1_atac_seq/SteveW_Rproject/atac_qc.Rmd',
    output_file='%(cwd)s/%(outfile)s')" --args "%(cwd)s/%(infile)s" "%(cwd)s/"''' % locals()
    
    P.run(statement)

#Count over .bam files
@merge(('*.bam',
       suffix(".bam")),
       "_counts.csv")
def atac_counting(infiles, outfile):
     '''Read counting with rsubread'''
     
     cwd = os.getcwd()
     peak_file =  PARAMS['peak_file_gtf']
     infiles2 = []
     
     for x in infiles:
         infiles2.append(x)
     
     length = len(infiles2)
     
     for x in range(length):
         infiles2[x] = cwd + "/" + str(infiles2[x])
         print (infiles2[x])
     
     infiles2. pop()
     infile_str = ','.join(infiles2)                   
     
     statement = '''R -e "rmarkdown::render('/ifs/projects/proj078/1_atac_seq/SteveW_Rproject/rsubread_multiple.Rmd',output_file='%(cwd)s/rsubread.htm')"  
     --args "%(infile_str)s"  "%(peak_file)s" "%(cwd)s/"''' % locals()
     
     P.run(statement)

#Use counts table and sample sheet for DEseq2
@transform('*_counts.csv',
       suffix("_counts.csv"),
       "deseq.html")
def atac_deseq(infile, outfile):
    ''''Differential anlaysis with deseq'''
    
    cwd = os.getcwd()
    sample_sheet = PARAMS['sample_sheet']
    peak_file =  PARAMS['peak_file']
    gene = PARAMS['gene']
    
    statement = '''R -e "rmarkdown::render('/ifs/projects/proj078/1_atac_seq/SteveW_Rproject/DEseq2_run.Rmd',output_file='%(cwd)s/deseq.htm')"  
     --args "%(cwd)s/%(infile)s" "%(sample_sheet)s" "%(peak_file)s" "%(gene)s"''' % locals()
     
    P.run(statement)

@transform('*_counts.csv',
       suffix("_counts.csv"),
       "deseq_pca.html")
def atac_deseq_pca(infile, outfile):
    ''''PCA anlaysis with deseq'''
    
    cwd = os.getcwd()
    sample_sheet = PARAMS['sample_sheet']
    
    statement = '''R -e "rmarkdown::render('/ifs/projects/proj078/1_atac_seq/SteveW_Rproject/DEseq_pca.Rmd',output_file='%(cwd)s/deseq_pca.htm')"  
     --args "%(cwd)s/%(infile)s" "%(sample_sheet)s"''' % locals()
     
    P.run(statement)
    
  # ---------------------------------------------------
  #Generic pipeline tasks
#@follows(atac_qc)
#def full():
#    pass
# 
# 
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)
# 
# 
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
# 
# =============================================================================
