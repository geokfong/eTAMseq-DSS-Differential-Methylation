# # step 1: Build the standard HISAT-3N index (with A to G conversion)
# Build an rRNA-specific index to filter rRNA reads during the analysis, can obtain rRNA fasta online
hisat-3n-build --base-change A,G rRNA.fasta rRNA

# Download and prepare the GRCh38 genome for HISAT-3N alignment
wget https://s3.amazonaws.com/igenomes.illumina.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz
tar -xvzf Homo_sapiens_NCBI_GRCh38.tar.gz
hisat-3n-build --base-change A,G Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa GRCh38

## Step 2: Define variables for the analysis
hisat3n_rep_index_name=rRNA
hisat3n_rep_index=/home/onggf/data/2_eTAM/rRNA/$hisat3n_rep_index_name
hisat3n_index_name=GRCh38
hisat3n_index=/home/onggf/data/2_eTAM/hg38/$hisat3n_index_name
hisat3n_out_dir=/home/onggf/data/2_eTAM/result
rep_fa=/home/onggf/data/2_eTAM/rRNA/rRNA.fasta
genome_fa=/home/onggf/data/2_eTAM/hg38/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
ncpus=24
A2G_percent=0.5
cov=1
strandness="F"
parent_folder='/home/onggf/data/2_eTAM/X401SC23022423-Z01-F003/01.RawData/'

# Step 3: Processing files
mamba activate eTAM

for file in "$parent_folder"*/; 
do
    echo ${file} 
    r1=$(find "$file" -type f -name "*1.fq.gz" | head -n 1)
    r2=$(find "$file" -type f -name "*2.fq.gz" | head -n 1)
    sample=$(basename "$(dirname "$r1")")  
    echo "$r1"
    echo "$r2"
    echo "$sample"

    cd /home/onggf/data/2_eTAM/X401SC23022423-Z01-F003/02_processed/
    mkdir -p "$sample"
    cd "$sample"

    # Step 3.1: Adapter trimming
    echo "adapter trimming: $sample"
    ## trim R1 5' end adapters and R2 3' end adapters
    cutadapt -e 0.1 -n 1 -O 1 -q 6 -m 0:46 -g GACGCTCTTCCGATCT -A AGATCGGAAGAGCGTC \
        -o "${sample}.R1.a5trim.fastq.gz" -p "${sample}.R2.a3trim.fastq.gz" "$r1" "$r2" > "${sample}.a3a5trim.umi.round1.log"
    ## trim R1 3' end adapters and R2 5' end adapters
    cutadapt -e 0.1 -n 1 -O 7 -q 6 -m 0:46 -a NNNNNNAGATCGGAAGAGCACA -G TGTGCTCTTCCGATCT \
        -o "${sample}.R1.a3a5trim.fastq.gz" -p "${sample}.R2.a3a5trim.fastq.gz" "${sample}.R1.a5trim.fastq.gz" "${sample}.R2.a3trim.fastq.gz" > "${sample}.a3a5trim.umi.round2.log"

    # Step 3.2: UMI extraction
    umi_tools extract --random-seed=123 --extract-method=regex --bc-pattern=".*" --bc-pattern2="^(?P<umi_1>.{6}).*" \
        -I "${sample}.R1.a3a5trim.fastq.gz" --read2-in="${sample}.R2.a3a5trim.fastq.gz" \
        --stdout="${sample}.R1.a3a5trim.umi.fastq.gz" --read2-out="${sample}.R2.a3a5trim.umi.fastq.gz" -L "${sample}.a3a5trim.umi.log"
    rm $sample.R1.a5trim.fastq.gz $sample.R2.a3trim.fastq.gz
    rm $sample.R[12].a3a5trim.fastq.gz
    
    mkdir -p ${hisat3n_out_dir}/${sample}
    cd ${hisat3n_out_dir}/${sample}
    
    # Step 3.3: rRNA mapping
    echo -e "\nhisat3n align -- rRNA mapping\n"
    hisat-3n -p $ncpus --time --base-change A,G --no-spliced-alignment --no-softclip --norc --no-unal --rna-strandness $strandness -x $hisat3n_rep_index -U /home/onggf/data/2_eTAM/X401SC23022423-Z01-F003/02_processed/${sample}/${sample}.R2.a3a5trim.umi.fastq.gz --un-gz $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz | \
        samtools sort -@ $ncpus -T $hisat3n_out_dir/$sample -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam -
    samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam
    samtools stats -r $rep_fa $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam > $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.stats
    samtools view -@ $ncpus -hb -e "([Yf]+[Zf]>0) && ([Yf]/([Yf]+[Zf])>=$A2G_percent)" $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam -
    /home/onggf/data/2_eTAM/pileup2var-1.1.0/pileup2var -f 524 -s $strandness -g $rep_fa -b $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.pileup2var.flt.txt

    # Step 3.4: Genome mapping
    echo -e "\nhisat3n align -- genome mapping\n"
    hisat-3n -p $ncpus --time --base-change A,G --repeat --repeat-limit 1000 --bowtie2-dp 0 --no-unal --rna-strandness $strandness -x $hisat3n_index -U $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz | \
        samtools view -@ $ncpus -Shb - -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam
    samtools view -@ $ncpus -q 60 -hb $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam -
        rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam
    samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam

    # Step 3.5: Deduplication
    echo -e "\deduplication\n"
    umi_tools dedup --random-seed=123 --method=unique --spliced-is-unique -I $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam -S $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam \
        --output-stats=$hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup -L $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.log
    samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam
    samtools stats -r $genome_fa $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam > $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.stats
    samtools view -@ $ncpus -hb -e "([Yf]+[Zf]>0) && ([Yf]/([Yf]+[Zf])>=$A2G_percent)" $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam -
    /home/onggf/data/2_eTAM/pileup2var-1.1.0/pileup2var -t $ncpus -f 524 -a A -c $cov -s $strandness -g $genome_fa -b $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.pileup2var.flt.txt

done
