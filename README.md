# FindClusters

This set of scripts allows you to find a specified gene cluster (provided in the form of a set of nucleotide alignments - one for each gene in the cluster) in a set of genome assemblies.

There is a shell script (runClusterPipeline.sh) that runs the whole analysis, running each script in sequence for all your input files. You can run this script from the command line (from the FindCLustersPipeline directory), giving your genome assemblies and target gene alignments as described below. If you want to run the pipeline from anywhere in your system, you will need to assign the various python scripts to environment variables and edit the shell script to reflect this. At the moment it is easiest to just run the pipeline from within its installation directory, as you can specify the full paths to your input and output files anyway.

The output image may be quite large. You may find that you will need to split the image into pieces to fit onto A4 for printing. There are several sites online that can do this. You may also need to reduce the resolution (I've set it quite high) or convert to PDF (can also be done online) to reduce your file sizes.



Please make sure that the assemblies have simple FASTA headers as the parser is not particularly sophisticated. Special characters (e.g. pipes, asterisks, etc.) will likely result in missing results. The best option is to have a header that starts with a unique alphanumeric identifier. Any other info can be separated by whitespace and will be ignored.

e.g. ">GI12345" is fine, while ">gi|1234|ref|NC0200019.1|" is not.



Before running the analysis for the first time, make sure you have installed HMMER (versions 3.1b2 - 3.2 definitely work). The Python dependencies should be installed automatically but these are the first things you should check if things go wrong.

To run the analysis, in the FindClustersPipeline directory, just run:

    bash runClusterPipeline.sh -g [assemblies] -a [alignments] -o [output] -i [image]

You can also access help, or run an example analysis, with the following commands:

    bash runClusterPipeline.sh -h

    bash runClusterPipeline.sh -e


You can also run each script individually, if you really want to. They can all be run from the command line, should you need their functionality separately. Run each script with the [-h] flag to see help messages for each one.

runHMMERsearch.py - as the name suggests, this takes your alignments and assemblies and creates profile HMMs for each gene, then searches the assemblies using those HMMs. It produces standard nhmmer output files (without the alignment section, for brevity) that you can use for other things, if needed.

mapCoordinates.py - takes nhmmer output files and retrieves the relevant sequences from the given genome assembly. It produces FASTA files of those sequences and can be configured to retrieve multiple sequences (rather than just the top hit), as well as adding a buffer sequence to the end of hits. It saves sequences to a strain specific directory.

findClusters.py - simply finds groups of genes/sequences that share a contig. Relies on the output of mapCoordinates.py as it just cycles through the strain directories checking the contig name in each filename. Outputs a tab-delimited file of strains that have at least two genes/sequences on the same contig.

drawClusters.py - reads the output file of findClusters.py and draws cluster diagrams for each strain (all one PNG image). You can optionally omit certain strains from the drawing using the --omitStrains option (for example to remove strains that have failed QC checks from a final figure).



FINAL NOTE
This pipeline will identify copies of gene clusters in large datasets but is not intended to be a final description of the clusters it finds. Gene boundaries and exact locations should be checked more thoroughly once they have been identified.