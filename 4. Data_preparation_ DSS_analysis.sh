## Step 6: Prepare dss input file
for gen in WT KO; do 
    for treatment in UT BT2 BT4; do 
        for rep in 1 2; do 
            input_1="/path/to/result/bed/result/final_output/${gen}_${treatment}_${rep}.ivtcount.table.out.txt"
            input_2="/path/to/result/bed/result/${gen}_${treatment}_${rep}.eTAMInput.tsv"
            output_file="/path/to/result/bed/result/final_output/dss/${gen}_${treatment}_${rep}.dss.txt"
            echo "Running ${gen}_${treatment}_${rep}"
            sort -k1,1 "$input_1" -o "$input_1"
            sort -k1,1 "$input_2" -o "$input_2"
            join -t $'\t' "$input_1" "$input_2" | awk -F"\t" 'BEGIN{print "chr\tpos\tN\tX"} {OFS="\t"; split($1, info, "_"); chr=info[1]; for(i=2; i<=length(info)-2; i++){chr=chr"_"info[i]}; print chr, info[length(info)-1], $4 + $5, int(($13)*($4+$5)/100)}' > "$output_file"
        done;
    done;
done
