#!/bin/bash

# This short script takes in a filename formatted as an EnsemblGeneID.txt
# It parses the file and construct a JSON-formatted sister file.

gene_filename=$1

if [[ -z ${gene_filename} ]]; then
    echo "No gene filename supplied. Expected: EnsemblGeneID.txt-formatted file"
    exit 1
fi

cat ${gene_filename} | jq . > ${gene_filename}.ensembl
rm ${gene_filename}
mv ${gene_filename}.ensembl ${gene_filename}