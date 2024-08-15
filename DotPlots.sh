####I need to organize and comment this file!







#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="hg38_dot_HG01946_reads"
#SBATCH --output=dot_HG01946_reads_hg38.out
#SBATCH --error=dot_HG01946_reads_hg38.err
##SBATCH --requeue
#SBATCH --cpus-per-task=12
#SBATCH --array=0-12 # job array index
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR
TMPDIR=/tmp
echo "TMPDIR="$TMPDIR

ulimit -s unlimited

#Ref="/projects/academic/omergokc/charikleia/t2t_chr.fa"
Ref="/projects/academic/omergokc/hg38.fa"
assembly="/projects/academic/omergokc/Luane/HG01946/nanoq-flye/results_nanoq_10kb_HG01946_ULR/HG01946_10kbreads.fastq"
rn="hg38"
#rn="t2t"
kb="reads_"
outdir="HG01946_hg38_reads_flye"

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"

conda activate Assembly

mkdir ${outdir}

minimap2 -x asm5 ${Ref} ${assembly} > ${outdir}/${rn}_vs_reads_HG01946_ULR_${kb}kb.paf

#module load gcc/11.2.0
#module load openmpi/4.1.1
#module load r/4.2.0

#cut -f6 ${outdir}/${rn}_vs_reads_HG01946_ULR_${kb}kb.paf| sort -n | uniq | while read line; do grep $line ${outdir}/${rn}_vs_reads_HG01946_ULR_${kb}kb.paf > ${outdir}/${rn}_HG01946_${kb}kb_${line}.paf; /projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_HG01946_${kb}kb_${line}.paf; done

/projects/academic/omergokc/Luane/softwares/paf2dotplot/paf2dotplot.r ${outdir}/${rn}_vs_reads_HG01946_ULR_${kb}kb.paf

conda deactivate
