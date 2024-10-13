#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="dotOpossum1_hifiasm_haploid"
#SBATCH --output=dothap1_Opossum_hifiasm_haploid.out
#SBATCH --error=dothap1_Opossum_hifiasm_haploid.err
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

####make sure the number of array matches the number of chromosomes.
##Also, run the second part of this script separately, and delete the arrays. 

######Defining references
Ref="/projects/academic/omergokc/Luane/Opossum1/ncbi_dataset_old/data/GCF_027887165.1/GCF_027887165.1_mMonDom1.pri_genomic.fna" #Reference genome to generate paf file
assembly_hap1="/projects/academic/omergokc/Luane/Opossum1/TGS/Opossum1_hap1.scaff_seqs" #path to your assembly
assembly_hap2="/projects/academic/omergokc/Luane/Opossum1/TGS/hap2/Opossum1_hap2.scaff_seqs" #path to your assembly

rn="GCF_027887165.1" #reference name, just for naming files
#rn="t2t"
outdir="Opossum1_hifiasm_GCF_027887165.1" #name of the out directory
id="Opossum1"

#I am activating my conda environment, because this is where I have all the programs installed.
eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate Assembly

#make outdirectory
mkdir -p ${outdir}_hap1
mkdir -p ${outdir}_hap2

#loading the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

#mapping the assembly to the reference and generating a paf file
minimap2 -x asm5 ${Ref} ${assembly_hap1} > ${outdir}_hap1/${rn}_vs_${id}_hap1.paf
minimap2 -x asm5 ${Ref} ${assembly_hap2} > ${outdir}_hap2/${rn}_vs_${id}_hap2.paf

############################## HAP1 #################################################### 
cd ${outdir}_hap1
n="1"

#this command cuts the sixth column of the paf file, that contains the chromosome to which it is aligned
#then it sorts the names of the chromosomes and while reading each chromosome, greps the entire line and save it to a separate file
#then automatically uses the paf2dotplot script on R to generate dotplots for the paf files that are separated by chromosome.
cut -f6  ${rn}_vs_${id}_hap${n}.paf| sort -n | uniq | while read line; do grep $line ${rn}_vs_${id}_hap${n}.paf > ${rn}_vs_${id}_hap1_${line}.paf; /projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_vs_${id}_hap${n}_${line}.paf; done

#this is just to also generate a dotplot for the file containing all the chromosomes.
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_vs_${id}_hap${n}.paf

############################## HAP2 #################################################### 
cd ../${outdir}_hap1
n="2"

#this command cuts the sixth column of the paf file, that contains the chromosome to which it is aligned
#then it sorts the names of the chromosomes and while reading each chromosome, greps the entire line and save it to a separate file
#then automatically uses the paf2dotplot script on R to generate dotplots for the paf files that are separated by chromosome.
cut -f6  ${rn}_vs_${id}_hap${n}.paf| sort -n | uniq | while read line; do grep $line ${rn}_vs_${id}_hap${n}.paf > ${rn}_vs_${id}_hap1_${line}.paf; /projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_vs_${id}_hap${n}_${line}.paf; done

#this is just to also generate a dotplot for the file containing all the chromosomes.
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_vs_${id}_hap${n}.paf

conda deactivate
