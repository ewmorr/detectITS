#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsxpress"
#SBATCH --output=ots_map.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
conda activate itsxpressenv
#we use itsxpress env just bc it has vsearch and we know it is already created

workDir=$HOME/$1
cd $workDir

indir=itsxpress_out
outdir=its_map
mkdir $outdir

for i in $indir/*R1*
do(
    dir=${i%/*}
    r1File=${i##*/}
    r2File=$(sed -e "s/_R1_/_R2_/" <<< "$r1File")
    sample_id=$(cut -d_ -f1 <<< "$r1File")
    
    vsearch --fastq_mergepairs $dir/$r1File \
        --reverse $dir/$r2File \
        --threads 24 \
        --fastaout $outdir/${sample_id}.merged.fasta

    vsearch --derep_fulllength $outdir/${sample_id}.merged.fasta \
        --sizeout \
        --relabel $sample_id. \
        --threads 24 \
        --output $outdir/${sample_id}.derep.fasta

    sed -i '' "s/^>.*/&;sample=${sample_id};/g" $outdir/${sample_id}.derep.fasta
)
done

cat $outdir/*.derep.fasta > $outdir/all_seqs.fa

db=~/blast_dbs/sh_general_release_dynamic_07252023/sh_general_allEuk_dynamic_singletons_25072023.ITS2.SPP_OF_CONCERN_CORRECTIONS.fasta

vsearch --usearch_global $outdir/all_seqs.fa \
    --db $db \
    --sizein \
    --id 0.985 \
    --threads 24 \
    --blast6out $outdir/all_seqs.blast6.tsv \
    --otutabout $outdir/all_seqs.otu_tab.tsv

grep -e "^#" -e "k__Fungi" $outdir/all_seqs.otu_tab.tsv |
    sed "s/#OTU ID/#OTUID/" |
    perl -n -e 'print if(/^#/ || /\w+.*?\t.*[1-9].*/)' > $outdir/all_seqs.otu_tab.hits.tsv
