#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="amy_HG01946_verkko_hg38"
#SBATCH --output=amy_HG01946_verkko_hg38.out
#SBATCH --error=amy_HG01946_verkko_hg38.err
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#This script uses the paf files you have generated with the dotplots scrips. 
#The purpose is to extract certain regions of the contigs from the assembly fasta file, and
#regionally align and visualize it using a reference genome.
#You will need the coordinates of the region of interest from the reference genome. 

#nucmer and mummerplot are part of this package: https://github.com/mummer4/mummer

######Defining references
#Ref="/projects/academic/omergokc/charikleia/t2t_chr.fa"#T2T
Ref="/projects/academic/omergokc/hg38.fa" #human hg38
assembly="/projects/academic/omergokc/Luane/HG01946/verkko/diploid/HG01946_verkko_set5/assembly.haplotype1.fasta" #path to your assembly
rn="hg38" #reference name, just for naming files
#rn="t2t"
outdir="HG01946_hg38_verkko" #name of the out directory

#I am activating my conda environment, because this is where I have all the programs installed.
eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate herro

#########For extracting and local alignment
#Grab the contigs that align to hg38 from 102mb-104mb (this is my region of interest) -- this is to visualize the names of the contigs
#this will use the paf file generated in the dotplot script. If you dont have that, please refer to previous scripts.
#awk '{if ($8>102000000 && $8<104000000) print $0}'  hg38_HG01946_chr1_only.paf > contigs_amy_region_mat.txt #to get the contigs, and the entire line fo the paf file.
#awk '{if ($8>102000000 && $8<104000000) print $0}'  hg38_HG01946_chr1_only.paf | cut -f1 | sort -n | uniq -c > contigs_amy_region.txt #in case you want to save the names of the contigs
awk '{if ($8>102000000 && $8<104000000) print $0}'  hg38_HG01946_chr1_only.paf | cut -f1,3,4 > contig_names.txt #Here it saves a file with column 1 contig name, 2, start and 3, end. see paf files columns for better understanding. 

#this command will read the file I just generated, and will save the name of the contig as i, the start and the end.
#then it will extract using samtools the exact coordinates of the contig from this assebly, and generate an output file with the name of the contig
while read -r i start end; do
  samtools faidx $assembly "${i}:${start}-${end}" > "${i}_${start}-${end}_amylocus.fa"
done < contig_names.txt 

#####in case you have less contigs and want to do it by hand
##from visualizing on the terminal, get the following:
#contig="pat-0010192" #name of contigs, you get this manually
#start="50829" #start
#edn="103" #end
#
##if you have a second contig, so on and so forth
#contig_2="pat-0010268" #name of contigs
#start_2="84143" #start
#edn_2="93" #end
#
##extracting the contigs only the region we want frim the fasta file
#samtools faidx $assembly ${contig}:${start}-${edn} > ${contig}_amylocus.fa
#samtools faidx $assembly ${contig_2}:${start_2}-${edn_2} > ${contig_2}_amylocus.fa
#
#concatenating both fasta files
cat *amylocus.fa > mat_contigs_amylocus.fa 

#align to the contigs to the assembly
#First we need to create a file with the region of interest from hg38 (extracted it from UCSC genome browser for example)
#then, she we map the specific region we extracted from the contig into the hg38 region of interest
#then we looked at the dotplots
nucmer  --prefix amy_hg38_hap1_verkko hg38_template.fa mat_contigs_amylocus.fa #nucmer is installed in the herro environment
#this command will create a delta file called amy_hg38_hap1_verkko.delta
mummerplot -prefix amy_hg38_hap1_verkko --png --large -R hg38_template.fa -Q mat_contigs_amylocus.fa amy_hg38_hap1_verkko.delta
#this command will create files fplot, rplot, and gp. We use the gp for creating the png 

module load gcccore gnuplot
gnuplot amy_hg38_hap1.gp #This will plot the graph and generate a png file.

conda deactivate
