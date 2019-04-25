#!/usr/bin/bash

# Install python dependencies (only needs doing once - comment out once installed).

pip install pyfaidx

pip install biopython

pip install Pillow

# Check hmmer installation.
hash nhmmer 2>/dev/null || { echo >&2 "Requirement: HMMER doesn't appear to be installed. Make sure it's added to your PATH.  Aborting."; exit 1; }

# Specify input files & output location.
# Assemblies must be in FASTA format (with .fa or .fasta extensions) and the filename must be the strain/organism (e.g. NCYC123.fa).
# Alignments must be in STOCKHOLM format (with .stockholm, .sto, or .stk extensions) and the filename must be the gene name (e.g. EMT1.sto).

# Please edit the 4 variables below to your desired input and output files.
# You can run the examples from the main script directory to see what is produced.

#assemblies="your assemblies here"
assemblies="./example/NCYC1384.FA ./example/NCYC3267.fasta ./example/NCYC3431.FASTA"

#alignments="your alignments here"
alignments="./example/EMT1.stockholm ./example/MAC1.stockholm ./example/MMF1.STO"

#output="your output directory"
output="."

#image="your desired image name (no extension)"
image="TEST"

echo "Running HMMER search."

python runHMMERsearch.py $output $assemblies $alignments

echo "Getting top hit FASTA files for matching sequences."

for assembly in $assemblies; do
    strain=$(basename -- "$assembly")
    strain="${strain%.*}"
    echo "Processing ..... $strain"
    files="$output"/hmmer_outputs/"$strain"/*.out
    for file in $files; do
        echo "Parsing ..... $file"
        python mapCoordinates.py $file $assembly --outputDirectory=$output --numHits=1 --geneNameOn=True
    done
done

echo "Finding clusters of matching sequences (same contig & within 10,000 nt)"

python findClusters.py "$output"/found_cluster_files "$output"/clusters.txt

echo "Drawing PNG image of clusters."

python drawClusters.py "$output"/clusters.txt "$output"/"$image"
