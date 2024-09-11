#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="dothap1_HG01946_verkko_hg38"
#SBATCH --output=dothap1_HG01946_verkko_hg38.out
#SBATCH --error=dothap1_HG01946_verkko_hg38.err
##SBATCH --requeue
#SBATCH --cpus-per-task=12
#SBATCH --array=0-12 # job array index
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

####make sure the number of array matches the number of chromosomes.
##Also, run the second part of this script separately, and delete the arrays. 

######Defining references
#Ref="/projects/academic/omergokc/charikleia/t2t_chr.fa"#T2T
Ref="/projects/academic/omergokc/hg38.fa" #human hg38
assembly="/projects/academic/omergokc/Luane/HG01946/verkko/diploid/HG01946_verkko_set5/assembly.haplotype1.fasta" #path to your assembly
rn="hg38" #reference name, just for naming files
#rn="t2t"
outdir="HG01946_hg38_verkko" #name of the out directory

#I am activating my conda environment, because this is where I have all the programs installed.
eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

#make outdirectory
mkdir ${outdir}

#loading the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

#mapping the assembly to the reference and generating a paf file
minimap2 -x asm5 ${Ref} ${assembly} > ${outdir}/${rn}_vs_verkko_HG01946.paf

#this command cuts the sixth column of the paf file, that contains the chromosome to which it is aligned
#then it sorts the names of the chromosomes and while reading each chromosome, greps the entire line and save it to a separate file
#then automatically uses the paf2dotplot script on R to generate dotplots for the paf files that are separated by chromosome.
cut -f6  ${outdir}/${rn}_vs_verkko_HG01946.paf| sort -n | uniq | while read line; do grep $line ${outdir}/${rn}_vs_verkko_HG01946.paf > ${outdir}/${rn}_HG01946_${line}.paf; /projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_HG01946_${line}.paf; done

#this is just to also generate a dotplot for the file containing all the chromosomes.
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_vs_verkko_HG01946.paf

#with this approach, chromosome 1 is not separated, because the command groups all chromosomes that start with 1. 
#because this is my chromosome of interest, I am writing a command that only generates a paf and dotplot for chromosome 1.
#only for chr1, because if not it'll be grouped with all the chr that start with 'chr1'
awk -F'\t' '$6 == "chr1"' ${outdir}/${rn}_vs_verkko_HG01946.paf > ${outdir}/${rn}_HG01946_chr1_only.paf #The
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_HG01946_chr1_only.paf
