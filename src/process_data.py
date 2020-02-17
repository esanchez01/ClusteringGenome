"""  Data Processing and Plotting

process_data.py processes VCF file data using PLINK2 and
plots final results.

"""

import pandas as pd
import numpy as np
import plotly.express as px
import subprocess as sp


def gather_fnames(data_fp):
    """
    Gathers all VCF files paths at a specified path and
    stores them in a input.list file at a specified path.
    
    :param data_fp: Path where data is housed
    """
    
    # Extracting filepaths for files of interest
    p1_1 = sp.Popen(['ls', data_fp.replace('\\', '')], stdout=sp.PIPE)
    p1_2 = sp.Popen(['grep', '.vcf.gz'], stdin=p1_1.stdout, stdout=sp.PIPE)
    p1_3 = sp.Popen(['grep', '-v', "tbi"], stdin=p1_2.stdout, stdout=sp.PIPE)

    f_1 = open("data/temp/input.list", "w")
    p1_4 = sp.Popen(['sed', '-e', 's/^/{}/'.format(data_fp)], stdin=p1_3.stdout, stdout=f_1)
    f_1.close()
    
    return


def concat_vcfs():
    """
    Gathers all VCF file paths inside a .list file, accesses 
    those paths, and concatenates all the corresponding VCF 
    files together.
    """
    
    # Storing VCF header 
    # TODO: Pipe safely (had some issues here with quoted commands)
    sp.Popen("zgrep '^\#' \"$(head -1 data/temp/input.list)\" | gzip > data/raw/full_chroms.vcf.gz", 
             shell=True, universal_newlines=True,
            stdout=sp.PIPE, stderr=sp.PIPE).communicate()

    # Concatenating VCF files
    sp.Popen("zgrep -v \"^\#\" $(cat data/temp/input.list) | gzip >> data/raw/full_chroms.vcf.gz", 
             shell=True, universal_newlines=True,
             stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    
    return


def filter_vcf(maf, geno, mind):
    """
    Filters VCF file using PLINK2
    
    :param maf: Minor allele frequency threshold
    :param geno: SNP missing call rate threshold
    :param mind: Sample missing call rate threshold
    """
    
    # Commands
    cmds = ['plink2', '--vcf', 'data/raw/full_chroms.vcf.gz', '--snps-only', 
            '--maf', str(maf), '--geno', str(geno), '--mind', str(mind), '--recode', 
            '--allow-extra-chr', '--make-bed', '--out', 'data/temp/chromosomes']
    
    # Filtering VCF
    sp.Popen(cmds, universal_newlines=True, 
             stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    
    return 


def check_outliers(eig_path):
    """
    Checks if outliers exist in eigenvector data
    produced from PLINK2's PCA output
    
    :param eig_path: Path where eigenvector data exists
    """
    
    # Reading in principal component data
    pcs = (pd.read_csv('{}/chrom_pc.eigenvec'.format(eig_path), sep=' ', header=None)
                   .rename(columns={1:'Sample', 2:'PC1', 3:'PC2', 4:'PC3'}))

    # Calculating z_scores
    z_scores = pcs[['PC1', 'PC2', 'PC3']].apply(lambda x: (x-x.mean())/x.std(), axis=0)

    # Finding outlier samples
    check_outlier = z_scores.apply(lambda x: abs(x) > 3, axis=0)
    outliers = check_outlier[check_outlier.any(axis='columns')].index
    outlier_samps = pcs.loc[outliers][[0, 'Sample']]
    
    # Checks if outliers exists
    if len(outlier_samps) > 0:
        # Writing outliers to text file
        outlier_samps.to_csv('data/temp/outliers.txt', sep= ' ', header=False, index=False)
        return True
    
    return False


def pca(num_pca, outpath, outliers=False):
    """
    Runs PCA on VCF file
    
    :param num_pca: Number of principle components
    :param outliers: Whether to filter out outliers or not
    :param outpath: Path to store eigenvector data
    """
    
    # Commands
    cmds = ['plink2', '--bfile', 'data/temp/chromosomes', '--pca', 
              str(num_pca), '--remove', 'data/temp/outliers.txt', 
            '--allow-extra-chr', '--out', '{}/chrom_pc'.format(outpath)]
    
    # Removes '--remove path' argument from command
    if not outliers: 
        cmds.pop(5)
        cmds.pop(6)
        
    # Running PCA
    sp.Popen(cmds, universal_newlines=True, 
             stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    
    return


def plot(outdir):
    """
    Plots PCA clusters
    
    :param outdir: Path where eigenvector data exists
    """
    
    # Reading in principal component data
    pcs = (pd.read_csv('{}/chrom_pc.eigenvec'.format(eig_path), sep=' ', header=None)
                   .rename(columns={1:'Sample', 2:'PC1', 3:'PC2', 4:'PC3'}))

    # Creating sample-superpopulation pair dictionary and 
    # mapping to samples in pc data above
    samples = pd.read_csv('data/raw/sample_pop.csv')
    samples_dict = samples.set_index('Sample').to_dict()['Population']
    pcs['Population'] = pcs.loc[:,0].apply(lambda x: samples_dict[x])

    # Plotting first three principle components
    fig = px.scatter_3d(pcs, x='PC1', y='PC2', z='PC3',color='Population')
    fig.update_traces(marker=dict(size=2.5))
    fig.update_layout(width=800, 
                      margin=dict(r=20, l=50, b=70, t=30), 
                      legend={'x':.75, 'y':.85, 'itemsizing':'constant'},
                      title={'text': "Clustering Chromosomes", 
                      'y':0.9,'x':0.53},
                      font=dict(size=15))
    
    fig.show()
    
    
# ---------------------------------------------------------------------
# Driver Function
# ---------------------------------------------------------------------

def process_data(inpath, maf, geno, mind, num_pca, outdir):
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # Gathers VCF filepaths
    gather_fnames(inpath)
    
    # Loads and concatentates files together
    concat_vcfs()
    
    # Filters master VCF file
    filter_vcf(maf, geno, mind)
    
    # Runs initial PCA
    pca(3, outdir, False)
    
    # Checks for outlier, reruns PCA if any exist
    if check_outliers(outdir):
        pca(3, outdir, True)
        
    # Plots clusters
    plot(outdir)
    
    
    #return