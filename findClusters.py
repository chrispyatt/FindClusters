#!$HOME/Programs/anaconda2/bin/python

from sys import argv
import argparse
import glob
import re
import itertools
from collections import defaultdict

parser = argparse.ArgumentParser(description='Input a directory containing (in subdirectories for each strain) sequence files for cluster genes.')
parser.add_argument('inputDir',
                    help='Directory containing ONLY MapCoordinates.py output files (FASTA files of gene sequences). These can be in subdirectories if multiple strains.')
parser.add_argument('outputFile',
                    help='Desired name of output .txt file (saves to current working directory unless full path specified).')
parser.add_argument('clusterLength',
                    help='The number of genes in the complete cluster.')

args = parser.parse_args()
inputDir = args.inputDir
outputFile = args.outputFile
clusterLength = int(args.clusterLength)

countFile = []

inDir = inputDir + "/*"
folders = glob.glob(inDir)

for folder in folders:
    counter = 0
    f = folder + "/*pos*.fasta"
    files = sorted(glob.glob(f), key=lambda x: re.split("[0-9.]+", x)[-2])
    
    str_files = ""
    strain = folder.split("/")[-1]
    for file in files:
        counter += 1
        splt = file.split("/" + strain + "/")[1]
        str_files = str_files + splt + " "
        #print(counter)
    toPrint = str(strain) + " : " + str(counter) + " : " + str(str_files) + "\n"
    countFile.append(toPrint)

linesToPrint = []

for line in countFile:
    splt = line.split(":")
    strain = splt[0].split("/")[-1].rstrip()
    count = splt[1].strip()
    genes = splt[2].strip().split(" ")
    
    geneInfo = {}
    
    for gene in genes:
        name = re.split("[_]", gene)[0]
        loc = re.split("[a-z]+", gene.split("_")[2])[1]
        coor = [re.split("[a-z.+-]+" ,gene.split("_")[3])[1], re.split("[a-z.+-]+" ,gene.split("_")[3])[2]]
        length = int(coor[1])-int(coor[0])
        dir = re.split("[0-9a-z.]+" ,gene.split("_")[3])[1]
        geneInfo[name] = [loc, dir, coor[0], coor[1], str(length)]
    
    #print("%s %s %s", (strain, count, geneInfo))

    # check if multiple genes appear on the same contig. Produce list of tuples where each tuple is 
    # a pair of "contig + number of genes on it".
    sameContig = [(k, len(list(v))) for k, v in itertools.groupby(sorted(list(geneInfo.values()), key = lambda x: int(x[0])), lambda y: y[0])]
    #print(sameContig)
    
    sameContigSorted = sorted(sameContig, key=lambda x: x[1])
    #print(sameContigSorted)
    
    # check if the genome contains all or nearly all (up to 2 missing) of the cluster genes, but all on different contigs (this picks up cases where a fragmented assembly prevents cluster detection).
    frag = False
    if len(geneInfo) == clusterLength or len(geneInfo) == clusterLength - 1 or len(geneInfo) == clusterLength - 2:
        frag = True
    
    if sameContigSorted[-1][1] > 1 or frag == True:
        toPrint = str(strain) + "\t" + str(sameContigSorted[-1][1]) + "\t" + str(count) + "\t" + str(len(geneInfo)) + "\t"
        for key,value in geneInfo.items():
            toPrint = toPrint + str(key) + ":" + str(value[0]) + ":" + str(value[1]) + ":" + str(value[2]) + ":" + str(value[3]) + ":" + str(value[4]) + " "
        toPrint.rstrip()
        linesToPrint.append(toPrint)

sortedByMaxCluster = sorted(linesToPrint, key=lambda x: int(x.split("\t")[1]), reverse=True)

titles = "Strain \t MaxCluster \t NumFiles \t NumGenes \t Genes"
sortedByMaxCluster.insert(0, titles)

with open(outputFile, "w") as outFile:
    for line in sortedByMaxCluster:
        outFile.write(line + "\n")
    
