#!/bin/bash
for file in /mnt/c/Users/samhu/source/repos/huettel-msc/export_folder/gaff_files/*.mol2; do
filename="$(basename $file)"
echo "$filename"
antechamber -i $file -fi mol2 -o /mnt/c/Users/samhu/source/repos/huettel-msc/export_folder/antechamber_output/$filename -fo mol2 -at gaff
done
