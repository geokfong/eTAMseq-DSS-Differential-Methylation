#Step 3: Data Preparation and Conversion to .bed
parent_folder='/home/onggf/data/2_eTAM/result/'

for file in "$parent_folder"*/; do
    echo ${file} 
    rna_pile=$(find "$file" -type f -name "*rRNA.pileup2var.flt.txt" | head -n 1)
    hg38_pile=$(find "$file" -type f -name "*GRCh38.pileup2var.flt.txt" | head -n 1)
    sample=$(basename "$(dirname "$rna_pile")")

    echo "$rna_pile"
    echo "$hg38_pile"
    echo "$sample"

    cd /home/onggf/data/2_eTAM/result/bed 
    mkdir -p "$sample"
    cd "$sample"

    ##rRNA conversion
    tail -n +2 "$rna_pile" | awk -F"\t" '{OFS="\t"; if(($4 == "A")){if($5 > 9){print $1,$2-1,$2,$1"="$2"="$3"="$5"="$6"="$9,$9/$5,$3}}}' - | awk -F"\t" '{OFS="\t"; print $1,($2-2),($3+2),$4,$5,$6}' | awk -F"\t" '{OFS="\t"; if($2>0){print $0}}' > "$sample.rna.tmp.bed"
    bedtools getfasta -fi /home/onggf/data/2_eTAM/rRNA/rRNA.fasta -bed "$sample.rna.tmp.bed" -nameOnly -tab -s | awk -F"\t" '{OFS="\t"; if(length($2)==5){a=toupper($2);gsub(/T/,"U",a); if(a~/[GAU][AG]AC[ACU]/) {m="DRACH"} else {m="nonDRACH"}; print $0"\t"m}}' | perl -pe 's/\([+-]\)//g'  | awk -F"\t" '{OFS="\t"; split($1, info, "="); print info[1],info[2]-1,info[2],$1"="$2"="$3,info[6]/info[4],info[3]}' > "$sample.rRNA.Conversions.bed"

    ##hg38 conversion
    tail -n +2 "$hg38_pile" | awk -F"\t" '{OFS="\t"; if(($4 == "A")){if($5 > 9){print $1,$2-1,$2,$1"="$2"="$3"="$5"="$6"="$9,$9/$5,$3}}}' - | awk -F"\t" '{OFS="\t"; print $1,($2-2),($3+2),$4,$5,$6}' | awk -F"\t" '{OFS="\t"; if($2>0){print $0}}' > "$sample.hg38.tmp.bed"
    bedtools getfasta -fi /home/onggf/data/2_eTAM/hg38/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa -bed "$sample.hg38.tmp.bed" -nameOnly -tab -s | awk -F"\t" '{OFS="\t"; if(length($2)==5){a=toupper($2);gsub(/T/,"U",a); if(a~/[GAU][AG]AC[ACU]/) {m="DRACH"} else {m="nonDRACH"}; print $0"\t"m}}' | perl -pe 's/\([+-]\)//g'  | awk -F"\t" '{OFS="\t"; split($1, info, "="); print info[1],info[2]-1,info[2],$1"="$2"="$3,info[6]/info[4],info[3]}' > "$sample.hg38.Conversions.bed"

done

## Step 4: Generate eTAM Input Files
parent_folder='/home/onggf/data/2_eTAM/result/bed/'

for gen in WT KO; do 
    for treatment in UT BT2 BT4; do 
        for rep in 1 2; do 
            echo Intersection of ${gen}_${treatment}_${rep} and IVT_${treatment}_${rep}; 
            intersectBed -wo -s -nonamecheck -a ${gen}_${treatment}_${rep}/${gen}_${treatment}_${rep}.hg38.Conversions.bed -b IVT_${treatment}_${rep}/IVT_${treatment}_${rep}.hg38.Conversions.bed |  awk -F '\t' 'BEGIN {OFS="\t";print "pos","motif","type","ftom_G_count","ftom_A_count","ivt_G_count","ivt_A_count"} {split($4,a,"=");split($10,b,"=");print a[1]"_"a[2]"_"a[3],a[8],a[7],a[6],a[5],b[6],b[5]}' > result/${gen}_${treatment}_${rep}.eTAMInput.tsv
        done; 
    done; 
done

#Step 5: Run R script to generate output
cd /home/onggf/data/2_eTAM/result/bed/result/
mkdir final_output

for gen in WT KO; do 
    for treatment in UT BT2 BT4; do 
        for rep in 1 2; do 
            input_file="/home/onggf/data/2_eTAM/result/bed/result/${gen}_${treatment}_${rep}.eTAMInput.tsv"
            output_file="/home/onggf/data/2_eTAM/result/bed/result/final_output/${gen}_${treatment}_${rep}.ivtcount.table.out.txt"
            echo "Running ${gen}_${treatment}_${rep}";
            Rscript run_model_ivt.R "$input_file" "$output_file" 50000 T 123 N 10 0.05 0.1
        done; 
    done; 
done
