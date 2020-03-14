""" File Converter

conversion.py contains functions that allow the
conversion of FASTQ files to BAM, BAM to VCF, or 
FASTQ straight to a VCF file.

"""

# Importing libraries
import subprocess as sp
import shutil
import os

REF_ROOT = 'Homo_sapiens_assembly38.fasta'
#PICARD = 'references/picard.jar'
R_GROUP = "@RG\\tID:group1\\tSM:Sample\\tPL:illumina\\tLB:lib1\\tPU:unit1"



def fastq_to_bam(fps, outdir, full=False, temp_kp=False):
    """
    Converts FASTQ file to BAM
    
    :param fps: List of paths to FASTQ files
    :param outdir: Path to output file
    :param full: Boolean whether it is a full 
                 conversion
    :param temp_kp: Whether to leave files in temp
                    directory. Used for full conversion.
    """
    
    # Prepare temp directory if not full 
    # FASTQ->VCF conversion
    if not full:
        prep_directory(fps, outdir, ftype='fastq')
    
    # Changing to temp directory
    temp = 'data/temp'
    os.chdir(temp)
    
    # Defining new filepaths 
    new_outdir = '../../'+outdir
    #new_picard = '../../'+PICARD
    sh_path = '../../src/conversion.sh'
    
    for file in fps:
        filename = file.split('/')[-1]
        filename_wo_ext = filename.split('.')[0]
        sam_path = filename_wo_ext+'.sam'
        bam_path = filename_wo_ext+'.bam'
        
        # Mapping FASTQs to reference file to create SAM
        sp.call(['sh', sh_path, 'mapper', R_GROUP, REF_ROOT, filename, sam_path])
        
        # Converting SAM to BAM
        sp.call(['sh', sh_path, 'sam_bam', sam_path, bam_path])
    
        # Writes to out directory
        if not full:
            sp.call(['sh', sh_path, 'copier', bam_path, new_outdir])
    
    # Clean out temp directory
    if not temp_kp:
        new_temp = '../temp'
        shutil.rmtree(new_temp, ignore_errors=True)
    
    return 



def bam_to_vcf(fps, outdir, full=False):
    """
    Converts a BAM file to VCF
    
    :param fps: List of paths to BAM files
    :param out_fp: Path to output file
    :param full: Boolean whether it is a full 
                 conversion
    """
    
    # Prepare temp directory if not full 
    # FASTQ->VCF conversion
    if not full:
        prep_directory(fps, outdir, ftype='bam')
        
        # Changing to temp directory
        temp = 'data/temp'
        os.chdir(temp)
        
        outdir = '../../'+outdir
    
    for file in fps:
        filename = file.split('/')[-1]
        filename_wo_ext = filename.split('.')[0]
        vcf_path = filename_wo_ext+'.vcf'
        sh_path = '../../src/conversion.sh'
        
        # Creating index for BAM files
        sp.call(['sh', sh_path, 'index_bam', filename, filename])
    
        # Converting BAM to VCF
        sp.call(['sh', sh_path, 'haplotype', REF_ROOT, filename, vcf_path])
    
        sp.call(['sh', sh_path, 'copier', vcf_path, outdir])
     
    # Clean out temp directory
    temp = '../temp'
    shutil.rmtree(temp, ignore_errors=True)
    
    return
    


def fastq_to_vcf(fps, outdir):
    """
    Converts FASTQ file to VCF
    
    :param fps: List of paths to FASTQ files
    :param outdir: Path to output file
    """
    
    # Prepare temp directory
    prep_directory(fps, outdir, ftype='fastq')
    
    # Converting FASTQ to BAM
    fastq_to_bam(fps, outdir, full=True, temp_kp=True)
    
    # Getting new names for input files
    fps_new = [fp.split('/')[-1].split('.')[0]+'.bam' 
               for fp in fps]   
    
    # Converting BAM files to VCF
    bam_to_vcf(fps_new, '../../'+outdir, True)
    
    return



def prep_directory(fps, outdir, ftype='fastq'):
    """
    Prepares temp directory for file conversion
    
    :param fp: Filepath of file to create directory
    :param outdir: Directory to write out files
    :param ftype: File type to prepare directory for
    """
    
    # Creating temp directory
    tmp = 'data/temp'
    if not os.path.exists(tmp):
        os.makedirs(tmp)
        
    # Filepaths to reference files
    ref_files_fp = ('../../../../datasets/dsc180a-wi20-public/'+
                    'Genome/resources/hg38/')
    
    # Shell script path
    sh_path = 'src/conversion.sh'
    
    # Moving reference to temp directory
    sp.call(['sh', sh_path, 'copier', ref_files_fp+REF_ROOT, tmp])
        
    # Moving input files to temp directory
    for file in fps:
        sp.call(['sh', sh_path, 'copier', file, tmp])
        
    # Creating sequence dictionary for reference
    unrooted = REF_ROOT.split('.')[0]
    sp.call(['sh', sh_path, 'seq_dict', tmp+'/'+REF_ROOT, tmp, unrooted])
    
    # Indexing reference
    sp.call(['sh', sh_path, 'fasta_index', tmp+'/'+REF_ROOT])
        
    if ftype == 'fastq':
        
        # Copying references genome files and input file to temp
        ref_ext = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        for ext in ref_ext:
            sp.call(['sh', sh_path, 'copier', ref_files_fp+REF_ROOT+ext, tmp])
        
    return



# ---------------------------------------------------------------------
# Driver Function
# --------------------------------------------------------------------- 


def convert_data(fastq_bam, fastq_vcf, bam_vcf, outdir, **kwargs):
    """
    Converts genetic data based on configuration file content.
    
    :param fastq_bam: List of filepaths to FASTQ files 
                      to convert to BAM
    :param fastq_vcf: List of filepaths to FASTQ files
                      to convert to VCF
    :param bam_vcf: List of filepaths to BAM file to 
                    convert to VCF
    :param outdir: Directory to write out converted files
    """
    
    # Creating out directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Convert FASTQ to BAM files
    if fastq_bam:
        fastq_to_bam(fastq_bam, outdir)
    
    # Convert FASTQ to VCF files
    if fastq_vcf:
        fastq_to_vcf(fastq_vcf, outdir)
    
    # Convert BAM to VCF files
    if bam_vcf:
        bam_to_vcf(bam_vcf, outdir)
        
    return 
    