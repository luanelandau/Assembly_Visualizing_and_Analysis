#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="dp1946hap_hg38_hifiasm_"
#SBATCH --output=dotplots_HG01946_hifiasm_haploid.out
#SBATCH --error=dotplots_HG01946_hifiasm_haploid.err
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#This script uses the paf files you have generated with the dotplots scrips. 
#The purpose is to extract certain regions of the contigs from the assembly fasta file, and
#regionally align and visualize it using a reference genome.
#You will need the coordinates of the region of interest from the reference genome. 
#nucmer and mummerplot are part of this package: https://github.com/mummer4/mummer

######Defining variables
Ref="/projects/academic/omergokc/hg38.fa" #Reference genome to generate paf file
Ref_cut="/projects/academic/omergokc/Luane/amy_hap/hg38_103570000-103760000.fa" #Fasta file with the region of the genome you're aiming to align to
assembly_hap1="/projects/academic/omergokc/Luane/HG01946/TGS-GapCloser/hifiasm_haploid/hap1/HG1946_hap1.scaff_seqs"
assembly_alt="/projects/academic/omergokc/Luane/HG01946/TGS-GapCloser/hifiasm_haploid/alternative/"

rn="hg38"
outdir="HG01946_haploid"
id="HG01946"

#I am activating my conda environment, because this is where I have all the programs installed.
eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate herro

#Make sure your output from hifiasm is converted using these commands:
#gfa tools is installed in the herro environment
#you can do this on terminal, cause it takes less than a minute
#gfatools gfa2fa $phased_assembly > HG01972.dip.hap1.p_ctg.fasta

##################Create paf files and dotplots for entire chr and assembly##############
mkdir -p ${outdir}
cd ${outdir}

minimap2 -x asm5 ${Ref} ${assembly_hap1} > ${rn}_vs_${id}_hap1.paf
minimap2 -x asm5 ${Ref} ${assembly_alt} > ${rn}_vs_${id}_alt.paf

module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

#to create one dotplot for the entire paf file
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_vs_${id}_hap1.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_vs_${id}_alt.paf


######################################### HAP1

#only for chr1, because if not it'll be grouped with all the chr that start with 'chr1'
#I am interested in chr1
awk -F'\t' '$6 == "chr1"' ${rn}_vs_${id}_hap1.paf > ${rn}_${id}_chr1_only_hap1.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_${id}_chr1_only_hap1.paf

awk '($8 <= 103760000) && ($9 >= 103570000)' ${rn}_${id}_chr1_only_hap1.paf > all_contigs_103570000-103760000_hap1.paf

#sorting first by the start position, then by contig name. So I have all contigs aligned in which the first ones to show will be the ones that start earlier
sort -k3 all_contigs_103570000-103760000_hap1.paf | sort -k1 > all_contigs_103570000-103760000_hap1_sorted.paf

#getting only the start and end of all contigs and the type of strand alignment (+ or -). This will be important in the next steps.
cut -f1,3,4,5 all_contigs_103570000-103760000_hap1_sorted.paf > all_contigs_103570000-103760000_hap1_sorted_start_end.paf

while read -r i start end strand; do
  samtools faidx $assembly_hap1 "${i}:${start}-${end}" > "${i}_${start}_${end}_hap1_amylocus.fa"
done < all_contigs_103570000-103760000_hap1_sorted_start_end.paf

cat *_hap1_amylocus.fa > allcontigs_hap1_amylocus.fa

nucmer  --prefix amylocus_${id}_hap1 $Ref_cut allcontigs_hap1_amylocus.fa
#The "-q" here filters for best alignment by contig. This means for all the possible alignments, it will filter by best identity and larger lenght of the alignment, giving you only one alignment per contig. The -l 1000 will give you alignments of at least 1000bp.
delta-filter -q -l 5000 amylocus_${id}_hap1.delta > amylocus_${id}_hap1.filter

show-coords amylocus_${id}_hap1.filter > amylocus_${id}_hap1.filter.coords.txt

delta="amylocus_${id}_hap1.filter.coords.txt"
awk 'BEGIN { FS="\\|"; OFS="\t"; } /^\s*$/ || /^={5,}/ || /^NUCMER/ || /^\s*\[S1\]/ || $0 !~ /\|/ || NF < 5 { next; } { for (i = 1; i <= NF; i++) { gsub(/^ +| +$/, "", $i); } split($2, pos2_fields, " +"); S2 = pos2_fields[1] + 0; E2 = pos2_fields[2] + 0; gsub(/^ +| +$/, "", $5); split($5, tag_fields, "\t"); contig_tag = tag_fields[2]; if (contig_tag == "" || contig_tag == "-") { next; } split(contig_tag, contig_parts, ":"); contig_name = contig_parts[1]; coords = contig_parts[2]; split(coords, coord_arr, "-"); contig_start = coord_arr[1] + 0; start = contig_start + (S2 - 1); end = contig_start + (E2 - 1); if (start > end) { temp = start; start = end; end = temp; } print contig_name, start, end; }' "${delta}" > contig_coordinates_alignments_nucmer_hap1.txt

while read -r i start end strand; do
  samtools faidx $assembly_hap1 "${i}:${start}-${end}" > "${i}_${start}_${end}_amylocus_hap1_nucmer.fa"
done < contig_coordinates_alignments_nucmer_hap1.txt

cat *_hap1_nucmer.fa > allcontigs_amylocus_hap1_nucmer.fa

######################################### ALT

#only for chr1, because if not it'll be grouped with all the chr that start with 'chr1'
#I am interested in chr1
awk -F'\t' '$6 == "chr1"' ${rn}_vs_${id}_alt.paf > ${rn}_${id}_chr1_only_alt.paf
/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${rn}_${id}_chr1_only_alt.paf

awk '($8 <= 103760000) && ($9 >= 103570000)' ${rn}_${id}_chr1_only_alt.paf > all_contigs_103570000-103760000_alt.paf

#sorting first by the start position, then by contig name. So I have all contigs aligned in which the first ones to show will be the ones that start earlier
sort -k3 all_contigs_103570000-103760000_alt.paf | sort -k1 > all_contigs_103570000-103760000_alt_sorted.paf

#getting only the start and end of all contigs and the type of strand alignment (+ or -). This will be important in the next steps.
cut -f1,3,4,5 all_contigs_103570000-103760000_alt_sorted.paf > all_contigs_103570000-103760000_alt_sorted_start_end.paf

while read -r i start end strand; do
  samtools faidx $assembly_alt "${i}:${start}-${end}" > "${i}_${start}_${end}_alt_amylocus.fa"
done < all_contigs_103570000-103760000_alt_sorted_start_end.paf

cat *_alt_amylocus.fa > allcontigs_alt_amylocus.fa

nucmer  --prefix amylocus_${id}_alt $Ref_cut allcontigs_alt_amylocus.fa
#The "-q" here filters for best alignment by contig. This means for all the possible alignments, it will filter by best identity and larger lenght of the alignment, giving you only one alignment per contig. The -l 1000 will give you alignments of at least 1000bp.
delta-filter -q -l 5000 amylocus_${id}_alt.delta > amylocus_${id}_alt.filter

show-coords amylocus_${id}_alt.filter > amylocus_${id}_alt.filter.coords.txt

delta="amylocus_${id}_alt.filter.coords.txt"
awk 'BEGIN { FS="\\|"; OFS="\t"; } /^\s*$/ || /^={5,}/ || /^NUCMER/ || /^\s*\[S1\]/ || $0 !~ /\|/ || NF < 5 { next; } { for (i = 1; i <= NF; i++) { gsub(/^ +| +$/, "", $i); } split($2, pos2_fields, " +"); S2 = pos2_fields[1] + 0; E2 = pos2_fields[2] + 0; gsub(/^ +| +$/, "", $5); split($5, tag_fields, "\t"); contig_tag = tag_fields[2]; if (contig_tag == "" || contig_tag == "-") { next; } split(contig_tag, contig_parts, ":"); contig_name = contig_parts[1]; coords = contig_parts[2]; split(coords, coord_arr, "-"); contig_start = coord_arr[1] + 0; start = contig_start + (S2 - 1); end = contig_start + (E2 - 1); if (start > end) { temp = start; start = end; end = temp; } print contig_name, start, end; }' "${delta}" > contig_coordinates_alignments_nucmer_alt.txt

while read -r i start end strand; do
  samtools faidx $assembly_alt "${i}:${start}-${end}" > "${i}_${start}_${end}_amylocus_alt_nucmer.fa"
done < contig_coordinates_alignments_nucmer_alt.txt

cat *_alt_nucmer.fa > allcontigs_amylocus_alt_nucmer.fa

##########################################################################################

cat *_nucmer.fa > allcontigs_amylocus_nucmer.fa

nucmer  --prefix amylocus_${id}_haploid_nucmer $Ref_cut allcontigs_amylocus_nucmer.fa

delta-filter -q -l 10000 amylocus_${id}_haploid_nucmer.delta > amylocus_${id}_haploid_nucmer.filter

module load gcc gnuplot

mummerplot -prefix amylocus_${id}_haploid_nucmer.filter --png --large -R $Ref_cut -Q allcontigs_amylocus_nucmer.fa amylocus_${id}_haploid_nucmer.filter

#########################################################################################
conda deactivate
