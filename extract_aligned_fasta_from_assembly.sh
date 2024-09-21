#!/bin/bash
#SBATCH --qos=omergokc
#SBATCH --partition=omergokc
#SBATCH --cluster=faculty
#SBATCH --account=omergokc
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="petarAssemblies"
#SBATCH --output=petarAssemblies_%A_%a.out
#SBATCH --error=petarAssemblies_%A_%a.err
##SBATCH --requeue
#SBATCH --cpus-per-task=6
#SBATCH --array=0-6 # job array index

#This script is used to look into a specific region of interest in your assembly. 
#First step is to align the assembly into a reference genome in which you know the coordinated you need. 
#Then, you generate a paf file with all the coordinates for the region you need and which contigs map
#to that region. After doing that, you extract from your assembly, only the parts of the contigs
#that you want (based on the mapping).
#This script is optimized for doing so for multiple assemblies at the same time. 
#You just have to create a file called "assemblies.txt" organized in the following way:
#HG01946	/projects/academic/omergokc/Luane/HG01946/ragtag/phased_verkko/ragtag_output/ragtag.scaffold.fasta
#HG02252	/projects/academic/omergokc/Luane/HG02252/ragtag/ragtag_output/ragtag.scaffold.fasta
#...
#meaning you have your assembly name in the first column and the path to the fasta file of the assembly in the second. 
#Then, you need to adjust the arrays to match exactly the number of assemblies you have (remember 
#the zero counts in this case). 
#At the end, it will generate an outdirectoty named as your assembly name
#and a fasta file called assemblyname/assemblyname_referencegenome_coordinates.fasta"
#Please remember to change your coordinates in the script. I need to optimize this.

# Read the assembly name and path from the file
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" assemblies.txt)
outdir=$(echo "$line" | awk '{print $1}')
assembly=$(echo "$line" | awk '{$1=""; print $0}' | sed 's/^ *//g')

#Just to make sure this is correctly reading
echo "Outdir: $outdir"
echo "Assembly: $assembly"

#The path for your reference genome and the name of the ref genome
Ref="/projects/academic/omergokc/hg38.fa"
rn="hg38"

#my conda environment that I have samtools and minimap2 installed
eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

#making out directory
mkdir -p ${outdir}

#mapping assembly to the reference genome and generating a paf file
minimap2 -x asm5 ${Ref} ${assembly} > "${outdir}/${outdir}_${rn}_mapped.paf"

#Regions to extract from hg38: chr1:155,186,123-155,192,868, modify accordingly
awk '$6 == "chr1" && $8 <= 155192868 && $9 >= 155186123' "${outdir}/${outdir}_${rn}_mapped.paf" > "${outdir}/contigs_petar_region_${outdir}_${rn}_155186123-155192868.txt"
awk '$6 == "chr1" && $8 <= 155192868 && $9 >= 155186123' "${outdir}/${outdir}_${rn}_mapped.paf" | cut -f1,3,4 > "${outdir}/contigs_names_petar_region_${outdir}_${rn}_155186123-155192868.txt"

# Extracting the region of interest from the contigs
while read -r i start end; do
  # Adjust start to be 1-based
  start=$((start + 1)) #this is adjusting for if the contig coordinate starts on zero
  samtools faidx "$assembly" "${i}:${start}-${end}" > "${outdir}/${i}_petarlocus.fa"
done < "${outdir}/contigs_names_petar_region_${outdir}_${rn}_155186123-155192868.txt"

# Concatenate all extracted regions into a single FASTA file
cat "${outdir}"/*_petarlocus.fa > "${outdir}/${outdir}_${rn}_155186123-155192868.fasta"

conda deactivate

