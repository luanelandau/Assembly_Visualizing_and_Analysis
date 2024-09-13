#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="dot_HG01972_hifiasm_hg38"
#SBATCH --output=dot_HG01972_hifiasm_hg38.out
#SBATCH --error=doth_HG01972_hifiasm_hg38.err
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#This script uses the paf files you have generated with the dotplots scrips. 
#The purpose is to extract certain regions of the contigs from the assembly fasta file, and
#regionally align and visualize it using a reference genome.
#You will need the coordinates of the region of interest from the reference genome. 
#nucmer and mummerplot are part of this package: https://github.com/mummer4/mummer

######Defining variables
#Ref="/projects/academic/omergokc/charikleia/t2t_chr.fa"
Ref="/projects/academic/omergokc/hg38.fa" #Reference genome to generate paf file
Ref_cut="/projects/academic/omergokc/Luane/HG01972/dotplots/hg38_103570000-103760000.fa" #Fasta file with the region of the genome you're aiming to align to
assembly_hap1="/projects/academic/omergokc/Luane/HG01972/hifiasm/HG01972.dip.hap1.p_ctg.gfa"
assembly_hap2="/projects/academic/omergokc/Luane/HG01972/hifiasm/HG01972.dip.hap2.p_ctg.gfa"
rn="hg38"
#rn="t2t"
kb="20"
outdir="HG01972_hg38_hifiasm_hap1"
outdir2="HG01972_hg38_hifiasm_hap2"

#I am activating my conda environment, because this is where I have all the programs installed.
eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

mkdir ${outdir}
mkdir ${outdir2}

minimap2 -x asm5 ${Ref} ${assembly_hap1} > ${outdir}/${rn}_vs_verkko_HG01946_hap1.paf
minimap2 -x asm5 ${Ref} ${assembly_hap2} > ${outdir2}/${rn}_vs_verkko_HG01946_hap2.paf

module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

#To create one dotplot per chromosome
#cut -f6  ${outdir}/${rn}_vs_verkko_HG01946.paf| sort -n | uniq | while read line; do grep $line ${outdir}/${rn}_vs_verkko_HG01946.paf > ${outdir}/${rn}_HG01946_${line}.paf; /projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_HG01946_${line}.paf; done

#to create one haplotype for the entire paf file
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_vs_verkko_HG01946_hap1.paf

/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir2}/${rn}_vs_verkko_HG01946_hap2.paf

#only for chr1, because if not it'll be grouped with all the chr that start with 'chr1'
#I am interested in chr1
awk -F'\t' '$6 == "chr1"' ${outdir}/${rn}_vs_verkko_HG01946_hap1.paf > ${outdir}/${rn}_HG01946_chr1_only.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_HG01946_chr1_only.paf

awk -F'\t' '$6 == "chr1"' ${outdir2}/${rn}_vs_verkko_HG01946_hap2.paf > ${outdir2}/${rn}_HG01946_chr1_only.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir2}/${rn}_HG01946_chr1_only.paf


#########For extracting and local alignment
#Grab the contigs that align to hg38 from 102mb-104mb -- this is to visualize the names of the contigs
#awk '{if ($8>102000000 && $8<104000000) print $0}'  hg38_HG01946_chr1_only.paf > contigs_amy_region_mat.txt
#awk '{if ($8>102000000 && $8<104000000) print $0}'  hg38_HG01946_chr1_only.paf | cut -f1 | sort -n | uniq -c > contigs_amy_region.txt #in case you want to save the names of the contigs
awk '{if ($8>102000000 && $8<104000000) print $0}' ${outdir}/${rn}_HG01946_chr1_only.paf | cut -f1,3,4 > ${outdir}/contig_names.txt #Here it saves a file with column 1 contig name, 2, start and 3, end.
awk '{if ($8>102000000 && $8<104000000) print $0}' ${outdir2}/${rn}_HG01946_chr1_only.paf | cut -f1,3,4 > ${outdir2}/contig_names.txt #Here it saves a file with column 1 contig name, 2, start and 3, end.

#this command will read the file I just generated, and will save the name of the contig as i, the start and the end.
#then it will extract using samtools the exact coordinates of the contig from this assebly, and generate an output file with the name of the contig
while read -r i start end; do
  samtools faidx $assembly_hap1 "${i}:${start}-${end}" > ${outdir}/"${i}_amylocus.fa"
done < ${outdir}/contig_names.txt 
#For the other assembly:
while read -r i start end; do
  samtools faidx $assembly_hap2 "${i}:${start}-${end}" > ${outdir2}/"${i}_amylocus.fa"
done < ${outdir2}/contig_names.txt 

#####in case you have less contigs and want to do it by hand
##from visualizing on the terminal, get the following:
#contig="pat-0010192" #name of contigs
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
cat ${outdir}/*_amylocus.fa > ${outdir}/all_contigs_amylocus.fa
cat ${outdir2}/*_amylocus.fa > ${outdir2}/all_contigs_amylocus.fa

conda activate herro
module load gnuplot #if not load this, mummer wont generate the png

#align to the contigs to the assembly
#First we need to create a file with the region of interest from hg38 (extracted it from UCSC genome browser for example)
#then, she we map the specific region we extracted from the contig into the hg38 region of interest
#then we looked at the dotplots
nucmer  --prefix amy_hg38_hap1_hifiasm $Ref_cut ${outdir}/all_contigs_amylocus.fa #nucmer is installed in the herro environment
nucmer  --prefix amy_hg38_hap2_hifiasm $Ref_cut ${outdir2}/all_contigs_amylocus.fa
#this command will create a delta file called amy_hg38_hap1_verkko.delta

#this command will create files fplot, rplot, and gp. We use the gp for creating the png 
mummerplot -prefix amy_hg38_hap1_hifiasm --png --large -R $Ref_cut -Q ${outdir}/all_contigs_amylocus.fa amy_hg38_hap1_hifiasm.delta
mummerplot -prefix amy_hg38_hap2_hifiasm --png --large -R $Ref_cut -Q ${outdir2}/all_contigs_amylocus.fa amy_hg38_hap2_hifiasm.delta


conda deactivate
