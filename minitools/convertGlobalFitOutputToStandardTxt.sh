#!/bin/bash

INFILE=$1
OUTFILE=$2

sed 1d $INFILE >onlyparameters.txt
awk ' { t = $1; $1 = $2; $2 = t; print; } ' onlyparameters.txt > temp.txt #swap first two rows
awk '$1="-"$1 {print}' temp.txt > tmp && mv tmp temp.txt #add minus to first column
awk '$2="-"$2 {print}' temp.txt > tmp && mv tmp temp.txt #add minus to first column
awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] }' temp.txt > tmp && mv tmp temp.txt #inverse order of rows
head -1 $INFILE > $OUTFILE
cat temp.txt >> $OUTFILE
awk ' { t = $1; $1 = t; print; } ' onlyparameters.txt > temp.txt #sync formating (even though it is not perfect)
cat temp.txt >> $OUTFILE


