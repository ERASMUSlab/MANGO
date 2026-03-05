#!/bin/bash
DIR=$1
TOTAL=$2
NCORE=$3

STEP=$(( (TOTAL + NCORE - 1) / NCORE ))
completed=0


for ((start=1; start<=TOTAL; start+=STEP)); do
    end=$(( start + STEP - 1 ))
    (( end > TOTAL )) && end=$TOTAL

    Rscript ../../GO_BB/MANGO_PREPROCESSING.R "$DIR" "$start" "$end" &
    while (( $(jobs -p | wc -l) >= NCORE )); do
        wait -n
        ((completed+=STEP))
        progress=$(( completed * 100 / TOTAL ))
        ((progress>100)) && progress=100
        echo -ne "\rProgress: ${progress}%"
    done
done

while (( $(jobs -p | wc -l) > 0 )); do
    wait -n
    ((completed+=STEP))
    progress=$(( completed * 100 / TOTAL ))
    ((progress>100)) && progress=100
    echo -ne "\rProgress: ${progress}%"
done

echo -e "\nDONE"

