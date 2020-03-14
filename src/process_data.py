"""  Data Processing and Plotting

process_data.py processes VCF file data using PLINK2 and
plots final results.

"""

# Importing libraries
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.offline as ply
import subprocess as sp
import os

SH_PATH = 'src/process_data.sh'



def gather_fnames(data_fp):
    """
    Gathers all VCF files paths at a specified path and
    stores them in a input.list file at a specified path.
    
    :param data_fp: Path where data is housed
    """
    
    input_fp = "data/temp/input.list"
    f_1 = open(input_fp, "w")
    
    # Extracting filepaths for files of interest
    clean_fp = data_fp.replace('\\', '')
    arg2 = 's/^/{}/'.format(data_fp)
    sp.call(['sh', SH_PATH, 'create_list', clean_fp, arg2, input_fp])
    
    f_1.close()

    return



def concat_vcfs():
    """
    Gathers all VCF file paths inside a .list file, accesses 
    those paths, and concatenates all the corresponding VCF 
    files together.
    """
    
    # Storing VCF header 
    sp.call(['sh', SH_PATH, 'header'])

    # Concatenating VCF files
    sp.call(['sh', SH_PATH, 'concatenate'])
    
    return



def filter_vcf(maf, geno, mind):
    """
    Filters VCF file using PLINK2
    
    :param maf: Minor allele frequency threshold
    :param geno: SNP missing call rate threshold
    :param mind: Sample missing call rate threshold
    """
    
    sp.call(['sh', SH_PATH, 'filter', str(maf), str(geno), str(mind)])
    
    return 



def check_outliers():
    """
    Checks if outliers exist in eigenvector data
    produced from PLINK2's PCA output  
    """

    # Reading in principal component data
    pcs = (pd.read_csv(
             'data/temp/chrom_pc.eigenvec', sep=' ', header=None)
             .rename(columns={1:'Sample', 2:'PC1', 3:'PC2', 4:'PC3'}))

    # Calculating z_scores
    z_scores = (pcs[['PC1', 'PC2', 'PC3']]
                   .apply(lambda x: (x-x.mean())/x.std(), axis=0))

    # Finding outlier samples
    check_outlier = z_scores.apply(lambda x: abs(x) > 3, axis=0)
    outliers = check_outlier[check_outlier.any(axis='columns')].index
    outlier_samps = pcs.loc[outliers][[0, 'Sample']]
    
    # Checks if outliers exists
    if len(outlier_samps) > 0:
        # Writing outliers to text file
        outlier_samps.to_csv('data/temp/outliers.txt', sep= ' ', 
                             header=False, index=False)
        return True

    return False



def pca(num_pca, outliers=False):
    """
    Runs PCA on VCF file
    
    :param num_pca: Number of principle components
    :param outliers: Whether to filter out outliers or not
    """
    
    cmds = ['sh', SH_PATH, 'pca', str(num_pca), '--remove', 
            'data/temp/outliers.txt', 'data/temp/chrom_pc']
    
    # Removes '--remove [path]' argument from command
    if not outliers: 
        cmds[4] = ''
        cmds[5] = ''

    # Running PCA
    sp.call(cmds)

    return



def create_trace(df, pop):
    """
    Helper function for 'plot'. Creates traces 
    for plot.
    
    :param df: DataFrame to create trace for
    :param pop: Population contained in df
    :returns: Trace
    """
    trace = go.Scatter3d(
        x=df['PC1'], y=df['PC2'], 
        z=df['PC3'], legendgroup="group", 
        name=pop, mode="markers",
        hovertext=df['Sample'], hoverinfo='text', 
        marker={'size':3}
    )
    
    return trace



def plot(outdir, test=False):
    """
    Plots PCA clusters
    
    :param outdir: Path to output PCA
    :param test: Boolean whether to generate test plot
    """
    
    # Reading in principal component data
    pcs = (pd.read_csv('data/temp/chrom_pc.eigenvec', sep=' ', header=None)
             .rename(columns={1:'Sample', 2:'PC1', 3:'PC2', 4:'PC3'}))

    # Creating sample-superpopulation pair dictionary and 
    # mapping to samples in pc data above
    samples = pd.read_csv('references/sample_pop.csv')
    samples_dict = samples.set_index('Sample').to_dict()['Population']
    pcs['Population'] = pcs.loc[:,0].apply(lambda x: samples_dict[x])
    pcs['Sample'] = pcs['Sample'].apply(lambda x: 'Sample: '+x)

    # Plotting first three principle components
    fig = go.Figure()

    # Subsetting data into populations
    european = pcs[pcs['Population'] == 'European']
    east_asian = pcs[pcs['Population'] == 'East Asian']
    mixed_american = pcs[pcs['Population'] == 'Ad Mixed American']
    south_asia = pcs[pcs['Population'] == 'South Asian']
    african = pcs[pcs['Population'] == 'African']

    # Adding data to plot
    fig.add_trace(create_trace(european, 'European'))
    fig.add_trace(create_trace(east_asian, 'East Asian'))
    fig.add_trace(create_trace(mixed_american, 'Ad Mixed American'))
    fig.add_trace(create_trace(south_asia, 'South Asian'))
    fig.add_trace(create_trace(african, 'African'))

    if test:
        sample = "Test Sample"
    else:
        sample = "Chromosome 1-22"
        
    # Customization
    fig.update_layout(width = 1000, margin = {'r':50, 'l':80, 'b':10, 't':50}, 
                      legend = {'x':.83, 'y':.85, 'itemsizing':'constant'},
                      title = {'text': ("<b>1000 Genomes SNP Mapping</b>"+
                                        "<br>({})".format(sample)), 
                               'y':0.9,'x':0.53},               
                      titlefont = {"size": 24},
                      font = {'size':12}, showlegend=True,
                      scene = {'xaxis_title': "PC1", 'yaxis_title':"PC2", 
                               'zaxis_title':"PC3"}
                     )
        
    plotname = '1000GenomesPlot.html'
    ply.plot(fig, filename='{}/{}'.format(outdir, plotname))
    print('***process finished, plot saved at {}/{}***'.format(outdir, plotname))
    
    
    
# ---------------------------------------------------------------------
# Driver Function
# ---------------------------------------------------------------------

def process_data(inpath, maf, geno, mind, num_pca, outdir, test=False):
    """
    Processes VCF files to generate PCA plot
    
    :param inpath: Location of VCF files
    :param maf: Minor allele frequency threshold
    :param geno: SNP missing call rate threshold
    :param mind: Sample missing call rate threshold
    :param num_pca: Number of principle components
    :param outdir: Directory to output final plot
    :param test: Boolean whether to generate test plot
    """
    
    # Creating out directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    # Creating temp directory
    tmp = 'data/temp'
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    
    # Gathers VCF filepaths
    gather_fnames(inpath)
    
    # Loads and concatentates files together
    concat_vcfs()
    
    # Filters master VCF file
    filter_vcf(maf, geno, mind)
    
    # Runs initial PCA
    pca(num_pca, False)
    
    # Checks for outlier, reruns PCA if any exist
    if check_outliers():
        pca(num_pca, True)
        
    # Plots clusters
    plot(outdir, test)
    
    return
