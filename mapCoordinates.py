#!$HOME/Programs/anaconda2/bin/python

from sys import argv
import argparse
import re
from subprocess import call
from collections import defaultdict
from Bio import SeqIO
from pyfaidx import Fasta

# this line assigns the arguments given at command line to variables that can be used within the script (the first is the name of the script itself)
parser = argparse.ArgumentParser(description='Parse a hmmer output file to take a specified number of hits (sequence coordinates) and extract those sequences (using samtools) from the relevant genome assembly.')
parser.add_argument('inputFile',
                    help='The nhmmer output file you wish to get information from. This should take the form [gene_strain-someothertext.fasta.out]. If not searching for a specific gene (i.e. CYP hunt) then [gene] can be omitted if --geneName is set to False.')
parser.add_argument('inputAssembly',
                    help='The genome assembly from which to extract sequence.')
parser.add_argument('--outputDirectory', dest='outputDirectory', default=".",
                    help='The directory in which to save the extracted sequence (a directory is automatically created for the strain). Default is current directory.')
parser.add_argument('--numHits', dest='numHits', default=10000,
                    help='The number of hits to extract from the nhmmer search (will only extract those above inclusion threshold). The default value is set to 10000, which should retrieve all matches above the threshold except in truly unusual circumstances.')
parser.add_argument('--extraSeq', dest='extraSeq', default=0,
                    help='The amount of extra sequence (on either side) you wish to retrieve with each hit. The default is 0 bases. If you intend to annotate the sequence, I recommend you ask for at least 1000 extra bases.')
parser.add_argument('--geneNameOn', dest='geneNameOn', default=False,
                    help='Toggle flag telling the program whether the strain name is preceeded by a gene name in the input filename. Default is False.')

args = parser.parse_args()
inputFile = args.inputFile
inputAssembly = args.inputAssembly
outputDirectory = args.outputDirectory
numHits = args.numHits
extraSeq = int(args.extraSeq)
geneNameOn = args.geneNameOn

#command1 = "samtools faidx " + inputAssembly
#call(command1, shell=True)
assemb_indexed = Fasta(inputAssembly)

lengths={}

indexFile = inputAssembly + ".fai"
#fhi = open(indexFile)
with open(indexFile) as fhi:
    for line in fhi:
        name = line.split()[0]
        length = line.split()[1]
        lengths[name] = length
    
#fhi.close()

#print "debug 1"
#print inputFile + " " + inputAssembly + " " + outputDirectory + " " + numHits

splitInput = inputFile.split("/")
fileName = splitInput[-1]
splitFilename = re.split("[-_.]", fileName)
if geneNameOn == "True" or geneNameOn == "true":
    #print geneNameOn
    geneName = splitFilename[0] + "_"
    strainName = splitFilename[1]
    #print geneName
    #print strainName
else:
    #print geneNameOn
    geneName = ""
    strainName = splitFilename[0]
    #print geneName
    #print strainName
    
#print "GENE=" + geneName + " STRAIN=" + strainName
    
#print "debug 2"
## mapCoordinates takes an input file (an nhmmer output file) and retrieves the sequence name,
## start coordinate, and end coordinate of sequences above the inclusion threshold (up to number supplied - numHits)

# open the file
#fh = open(inputFile)
with open(inputFile) as fh:
    # initiate list [coords] and dictionary [matches] to store matching sequences
    coords = []
    matches = defaultdict(list)
    
    #print "debug 3"
    
    # loop through lines of file until you get to the matching sequences, then extract the required number of sequences (until inclusion threshold is hit), stored as a list of lists (coords).
    inclusion = False
    numSearched = 0
    for line in fh:
        #print line
        if line.startswith("    ------- ------ -----  --------  -----  -----  -----------") or "".join(line.split()) == "-----------------------------------------------":
            inclusion = True
            continue
        if line.startswith("  ------ inclusion threshold ------"):
            inclusion = False
        if inclusion == True and (line.startswith("Annotation for each hit  (and alignments):") or line == "\n"):
            inclusion = False
            break
        if inclusion == True and numSearched < int(numHits):
            fields = line.split()
            #print "FIELDS= " + str(fields)
            contigName = fields[3]
            start = int(fields[4])
            end = int(fields[5])
            length = int(lengths[contigName])
            
            #print fields
            if start > end:
                # modify start/ end coordinates according to extraSeq argument. max() prevents this from dipping below zero in the event of the sequence being near the start of a contig.
                # samtools does not crash if the end coordinate is past the end of a contig but it is helpful to also ensure this does not happen (file header is affected later).
                end = max(1, (end - extraSeq))
                start = min(length, (start + extraSeq))
                lst = [contigName, end, start, "-"]
            else:
                start = max(1, (start - extraSeq))
                end = min(length, (end + extraSeq))
                lst = [contigName, start, end, "+"]
            #matches[fields[3]] = lst
            coords.append(lst)
            numSearched = numSearched + 1

#fh.close()
#return matches
#print coords

#print "debug 4"
#print mapCoordinates(inputFile)
#matches = mapCoordinates(inputFile)

# store matching sequences as dictionary entries with the contig name (as found in nhhmer file) as the key, and a list comprising the coordinates and direction as the value.
for i, j, k, l in coords:
    matches[i].append((j, k, l))
    
#print "debug 5"
#print matches

# if sequences have been found, make a directory (named for the strain) in which to save them if one does not already exist.
if matches:
    command2 = "mkdir -p " + outputDirectory + "/found_cluster_files/" + strainName
    #print command1
    call(command2, shell=True)
    
#print "debug 6"

# for each key (contig) in 'matches', run 'samtools faidx' (index the fasta file) followed by 'samtools faidx contig:start-end' (pull out sequence) and send that to appropriate file.
for contig in matches:
    #print contig, matches.get(contig) #matches.get(match)[0], matches.get(match)[1]
    
    # another for-loop (through list of lists in case multiple matches for one contig)
    for match in matches.get(contig):
        #print match, type(str(match))
        location = strainName + "_" + contig + ":" + str(match[0]) + "-" + str(match[1])
        outputFile = outputDirectory + "/found_cluster_files/" + strainName + "/" + geneName + strainName + "_contig" + contig + "_pos" + str(match[0]) + "to" + str(match[1]) + str(match[2]) + ".fasta"
        #print "OUTPUTFILE= " + outputFile
        #command3 = "samtools faidx " + inputAssembly + " " + location + " > " + outputFile
        #call(command3, shell=True)
        seq_record = assemb_indexed.get_seq(contig, match[0], match[1])
        
        with open(outputFile, "w") as eFn:
            eFn.write(">" + location + "\n")
            eFn.write(str(seq_record) + "\n")
        
        # editFASTAnames functionality added here (add strain name to FASTA header of outputFile - missed out by samtools)
        #eFn = open(outputFile, "r")
        #with open(outputFile, "r") as eFn:
        #    for seq_record in SeqIO.parse(eFn, "fasta"):
        #        old_header = seq_record.id
        #        #print old_header, old_header.type
        #        splitHeader = old_header.split("_")
        #        if splitHeader[0] != strainName:
        #            new_header = strainName + "_" + old_header
        #            seq_record.id = new_header
        #            seq_record.description = ""
        #eFn.close()

        #eFn = open(outputFile, "w")
        #with open(outputFile, "w") as eFn:
        #    SeqIO.write(seq_record, eFn, "fasta")
        #eFn.close()

        # if the sequence is on the opposite strand (indicated by '-' flag in filename), revComp functionality here to reverse complement the sequence, modify the header, and save to a new file.
        if (str(match[2]) == "-"):
            #rC = open(outputFile, "r")
            with open(outputFile, "r") as rC:
                for seq_record in SeqIO.parse(rC, "fasta"):
                    new_seq_record = seq_record.reverse_complement()
            #rC.close()

            newOutputFile = outputDirectory + "/found_cluster_files/" + strainName + "/" + geneName + strainName + "_contig" + contig + "_pos" + str(match[1]) + "to" + str(match[0]) + "+.fasta"
            #print newOutputFile
            # output filename is modified input filename (coordinates reversed)
            #rCo = open(newOutputFile, "w")
            with open(newOutputFile, "w") as rCo:
                new_header = strainName + "_" + contig + ":" + str(match[1]) + "-" + str(match[0])
                new_seq_record.id = new_header
                new_seq_record.description = ""
                SeqIO.write(new_seq_record, rCo, "fasta")
            #rCo.close()

#print "debug 7"