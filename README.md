# vsearch based read mapping for fungal species detection

- reads have been primer-trimmed, ITS-extracted, and quality filtered as performed in DADA2 pipeline
- followed by read merging, dereplication (with counting)
- mapping with --usearch_global against UNITE ITS2 db

vsearch v2.29.2_macos_x86_64
```
conda activate its_sim
cd ~/repo/detectITS/data
```
Merge quality filtered and ITS extracted sequences with vsearch and dereplicate.
```
outdir=its_merged
mkdir $outdir

for i in itsxpress_out/*R1*
do(
    dir=${i%/*}
    r1File=${i##*/}
    r2File=$(sed -e "s/_R1_/_R2_/" <<< "$r1File")
    sample_id=$(cut -d_ -f1 <<< "$r1File")
    
    vsearch --fastq_mergepairs $dir/$r1File \
        --reverse $dir/$r2File \
        --threads 4 \
        --fastaout $outdir/${sample_id}.fasta
)
done
```
Derpelicate sequences with cluster size output, and also label sequences with sample IDs
```
indir=its_merged
outdir=merged_derep
mkdir $outdir

for i in $indir/*.fasta
do(
    sample_id=${i#*/}
    sample_id=${sample_id%.fasta}

    vsearch --derep_fulllength $indir/${sample_id}.fasta \
        --sizeout \
        --relabel $sample_id. \
        --output $outdir/${sample_id}.fasta

    sed -i '' "s/^>.*/&;sample=${sample_id};/g" $outdir/${sample_id}.fasta
)
done
```
Cluster at 100% sim within each sample and compare to derrrepped
```
indir=merged_derep
outdir=cluster_100
mkdir $outdir

for i in $indir/*.fasta
do(
    sample_id=${i#*/}
    sample_id=${sample_id%.fasta}

    vsearch --cluster_size $indir/${sample_id}.fasta \
        --threads 4 \
        --sizein \
        --sizeout \
        --id 1 \
        --centroids $outdir/${sample_id}.fasta
)
done

for i in $indir/*.fasta
do(
    sample_id=${i#*/}
    sample_id=${sample_id%.fasta}

    derepNum=$(echo $(cat $indir/${sample_id}.fasta | grep ">" | wc -l))
    otuNum=$(echo $(cat $outdir/${sample_id}.fasta | grep ">" | wc -l))
    echo -e "$sample_id\t$derepNum\t$otuNum"
)
done

```
Cat merged, derepped, labeled sequences to single file for mapping (and separate for 100%)
```
indir=merged_derep
cat $indir/* > all_seqs.derep.fa
indir=cluster_100
cat $indir/* > all_seqs.100.fa

grep ">" all_seqs.derep.fa | wc -l
grep ">" all_seqs.100.fa | wc -l
```
Use `vsearch --usearch_global` mapping to construct otutab
See https://www.drive5.com/usearch/manual/blast6out.html for an explanation of the columns in blast6out format (this gives the alignment statistics)

```
db=~/blast_dbs/sh_general_release_dynamic_07252023/sh_general_allEuk_dynamic_singletons_25072023.ITS2.SPP_OF_CONCERN_CORRECTIONS.fasta
less $db

vsearch --usearch_global all_seqs.derep.fa \
    --db $db \
    --sizein \
    --id 0.985 \
    --blast6out all_seqs.blast6.tsv \
    --otutabout all_seqs.otu_tab.tsv
    
vsearch --usearch_global all_seqs.100.fa \
    --db $db \
    --sizein \
    --id 0.985 \
    --blast6out all_seqs.100.blast6.tsv \
    --otutabout all_seqs.100.otu_tab.tsv

```
remove plants and metazoa from tab (saves time) and then search for rows with at least one non-zero value
```
grep -e "^#" -e "k__Fungi" all_seqs.otu_tab.tsv | 
    sed "s/#OTU ID/#OTUID/" | 
    perl -n -e 'print if(/^#/ || /\w+.*?\t.*[1-9].*/)' > all_seqs.otu_tab.hits.tsv
wc -l all_seqs.otu_tab.hits.tsv
#2415

grep -e "^#" -e "k__Fungi" all_seqs.100.otu_tab.tsv | 
    sed "s/#OTU ID/#OTUID/" | 
    perl -n -e 'print if(/^#/ || /\w+.*?\t.*[1-9].*/)' > all_seqs.100.otu_tab.hits.tsv
wc -l all_seqs.100.otu_tab.hits.tsv
#2398
```
There are less hits to *different* UNITE entries after 100% clustering. This means that real variation is lost. Likely a substr issue. Use derep

get just bretziella hits
```
grep -e "^#" -e "Bretziella" all_seqs.otu_tab.hits.tsv > BRFA_hits.tsv
grep -e "^#" -e "Bretziella" all_seqs.100.otu_tab.hits.tsv > BRFA_hits.100.tsv
```
