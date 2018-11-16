This set of scripts allows you to find a specified gene cluster (provided in the form of a set of nucleotide alignments - one for each gene in the cluster) in a set of genome assemblies.

There is a shell script (runClusterPipeline.sh) that runs the whole analysis, running each script in sequence for all your input files. You will need to edit this script to specify your input and output preferences.

Before running the analysis for the first time, make sure you have installed HMMER (versions 3.1b2 - 3.2 definitely work). The Python dependencies should be installed automatically but these are the first things you should check if things go wrong.

To run the analysis, in the FindClustersPipeline directory, just type:

    bash runClusterPipeline.sh


You can also run each script individually, if you really want to. They can all be run from the command line, should you need their functionality separately.

runHMMERsearch.py - as the name suggests, this takes your alignments and assemblies and creates profile HMMs for each gene, then searches the assemblies using those HMMs. It produces standard nhmmer output files (without the alignment section, for brevity) that you can use for other things, if needed.

mapCoordinates.py - takes nhmmer output files and retrieves the relevant sequences from the given genome assembly. It produces FASTA files of those sequences and can be configured to retrieve multiple sequences (rather than just the top hit), as well as adding a buffer sequence to the end of hits. It saves sequences to a strain specific directory.

findClusters.py - simply finds groups of genes/sequences that share a contig. Relies on the output of mapCoordinates.py as it just cycles through the strain directories checking the contig name in each filename. Outputs a tab-delimited file of strains that have at least two genes/sequences on the same contig.

drawClusters.py - reads the output file of findClusters.py and draws cluster diagrams for each strain (all one PNG image).