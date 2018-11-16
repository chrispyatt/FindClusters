#!$HOME/Programs/anaconda2/bin/python

from sys import argv
import argparse
import re
from subprocess import call
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This program takes a set of genome assemblies, in FASTA format, and a set of nucleotide alignments of cluster genes, in STOCKHOLM format. It uses HMMER to create HMMs for each gene, then uses those to search for matching sequences in the provided assemblies.')
#parser.add_argument('genomeAssemblies',
#                    help='FASTA assemblies to be searched')
#parser.add_argument('clusterGenes',
#                    help='STOCKHOLM alignments of cluster genes (one alignment per file). Files should be named in the style [gene].stockholm')
parser.add_argument('outDir',
                    help='Directory in which to save output.')
parser.add_argument('filenames', nargs='+',
                    help='Assembly and alignment files (will be separated by extension). Assembly files should be named in the style \'[strain].fa\'. Cluster gene files should be named in the style \'[gene].stockholm\'.')

args = parser.parse_args()

#print(args.filenames)

def getByExtension(lst, ext):
  return [f for f in lst if f.endswith(ext)]

genomeAssemblies = getByExtension(args.filenames, '.fa')
for i in getByExtension(args.filenames, '.fasta'):
    genomeAssemblies.append(i)
for i in getByExtension(args.filenames, '.FA'):
    genomeAssemblies.append(i)
for i in getByExtension(args.filenames, '.FASTA'):
    genomeAssemblies.append(i)

clusterGenes = getByExtension(args.filenames, '.stockholm')
for i in getByExtension(args.filenames, '.STOCKHOLM'):
    clusterGenes.append(i)
for i in getByExtension(args.filenames, '.sto'):
    clusterGenes.append(i)
for i in getByExtension(args.filenames, '.STO'):
    clusterGenes.append(i)
for i in getByExtension(args.filenames, '.stk'):
    clusterGenes.append(i)
for i in getByExtension(args.filenames, '.STK'):
    clusterGenes.append(i)
#genomeAssemblies = args.genomeAssemblies
#clusterGenes = args.clusterGenes
outDir = args.outDir

#print(genomeAssemblies)
#print(clusterGenes)
#print(outDir)

## run hmmbuild on each gene

geneList = []

command = "mkdir -p " + outDir + "/profile_hmms"
print(command)
call(command, shell=True)

for gene in clusterGenes:
    geneName = re.split("[/.]", gene)[-2]
    #print(gene)
    #print(re.split("[/.]", gene))
    geneList.append(geneName)
    command = "hmmbuild " + outDir + "/profile_hmms/" + geneName + ".hmm " + gene
    print(command)
    call(command, shell=True)

#print(geneList)

## run nhmmer on each genome for each gene

for assembly in genomeAssemblies:
    strainName = re.split("[/.]", assembly)[-2]
    #print(strainName)
    command1 = "mkdir -p " + outDir + "/hmmer_outputs/" + strainName
    print(command1)
    call(command1, shell=True)
    for gene in geneList:
        #print(gene)
        command2 = "nhmmer --noali -o " + outDir + "/hmmer_outputs/" + strainName + "/" + gene + "_" + strainName + ".out " + outDir + "/profile_hmms/" + gene + ".hmm " + assembly
        print(command2)
        call(command2, shell=True)


