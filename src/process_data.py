
import pandas as pd
import numpy as np
import plotly.express as px
import subprocess as sp


def gather_fnames(data_fp, temp_fp):
    
    # Extracting filepaths for files of interest
    p1_1 = sp.Popen(['ls', fp.replace('\\', '')], stdout=sp.PIPE)
    p1_2 = sp.Popen(['grep', '.vcf.gz'], stdin=p1_1.stdout, stdout=sp.PIPE)
    p1_3 = sp.Popen(['grep', '-v', "tbi"], stdin=p1_2.stdout, stdout=sp.PIPE)

    f_1 = open("{}input.list".format(store_fp), "w")
    p1_4 = sp.Popen(['sed', '-e', 's/^/{}/'.format(data_fp)], stdin=p1_3.stdout, stdout=f_1)
    f_1.close()
    
    return


def concat_vcfs():
    
    # Storing VCF header 
    # TODO: Pipe correctly (had some issues here with quoted commands)
    sp.Popen("zgrep '^\#' \"$(head -1 data/temp/input.list)\" | gzip > data/raw/full_chroms.vcf.gz", 
             shell=True, universal_newlines=True,
            stdout=sp.PIPE, stderr=sp.PIPE).communicate()

    # Concatenating VCF files
    sp.Popen("zgrep -v \"^\#\" $(cat data/temp/input.list) | gzip >> data/raw/full_chroms.vcf.gz", 
             shell=True, universal_newlines=True,
             stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    
    return


def filter_vcf(maf, geno, mind):
    
    # Filtering VCF
    sp.Popen(['plink2', '--vcf', 'data/raw/full_chroms.vcf.gz', '--snps-only', '--maf', str(maf), 
              '--geno', str(geno), '--mind', str(mind), '--recode', '--allow-extra-chr', 
              '--make-bed', '--out', 'temp/chromosomes'], 
             universal_newlines=True, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    
    return 

def pca(num_pca):
    
    # Running PCA
    sp.Popen(['plink2', '--bfile', 'data/temp/chromosomes', '--pca', 
              str(num_pca), '--allow-extra-chr', '--out', 'data/out/chrom_pc'], 
              universal_newlines=True, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    
    return
    
    
# ---------------------------------------------------------------------
# Driver Function
# ---------------------------------------------------------------------

def process_data(data_fp, maf, geno, mind, num_pca):
    return