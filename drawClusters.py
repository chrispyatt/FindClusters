from sys import argv
import argparse
import itertools
import numpy
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from PIL import ImageEnhance
import colorsys

#im = Image.new("RGBA", (800,200), (255,255,255,255))

'''
This script takes a tab-delimited file produced by findClusters.py. The drawMultipleClusters() function reads the input file and creates a canvas on which to draw the gene clusters specified in the input file. It then calls drawCluster() to draw each gene cluster, one underneath the other. drawCluster() itself calls the drawArrow() method to draw each arrow (representing each gene) onto a black line representing the DNA strand.

Font used in image is 'Academic M54', downloaded from https://www.wfonts.com/font/academic-m54. It is licensed as 'Free for personal use', covering this script since it is not exploited for profit. No copyright infringement intended.
'''

parser = argparse.ArgumentParser(description='Input a clusters file. Output a drawing of those clusters.')
parser.add_argument('inputFile',
                    help='Clusters file created by findClusters.py.')
parser.add_argument('outputFile',
                    help='Desired name of PNG cluster image (saves to current working directory unless full path specified).')
parser.add_argument('--omitStrains', dest='omitStrains', default="",
                    help='The names of any strains you want to omit from the final drawing (e.g. strains that have failed subsequent quality control tests).')

args = parser.parse_args()
inputFile = args.inputFile
outputFile = args.outputFile
omitStrains = args.omitStrains

def drawArrow(im, direction, colour, position, label, length):
    draw = ImageDraw.Draw(im)
    fnt = ImageFont.truetype('AcademicM54.ttf', 45)
    if direction == "+":
        draw.polygon([position, (position[0]+length,position[1]), (position[0]+length,position[1]-60), (position[0]+length+100,position[1]+40), (position[0]+length,position[1]+140), (position[0]+length,position[1]+80), (position[0],position[1]+80)], fill=colour)
        ImageDraw.Draw(im).text((position[0]+length/3,position[1]+160), label, font=fnt, fill="black")
    elif direction == "-":
        draw.polygon([(position[0],position[1]+40), (position[0]+100,position[1]-60), (position[0]+100,position[1]), (position[0]+100+length,position[1]), (position[0]+100+length,position[1]+80), (position[0]+100,position[1]+80), (position[0]+100,position[1]+140)], fill=colour)
        ImageDraw.Draw(im).text((position[0]+length/3,position[1]+160), label, font=fnt, fill="black")
    else: print("Invalid direction: Specify either \"+\" or \"-\".")
    del draw

### old version
def drawArrow2(im, direction, colour, position, label):
    draw = ImageDraw.Draw(im)
    if direction == "+":
        draw.polygon([position, (position[0]+50,position[1]), (position[0]+50,position[1]-15), (position[0]+75,position[1]+10), (position[0]+50,position[1]+35), (position[0]+50,position[1]+20), (position[0],position[1]+20)], fill=colour)
        ImageDraw.Draw(im).text((position[0]+30,position[1]+40), label, fill="black")
    elif direction == "-":
        draw.polygon([(position[0],position[1]+10), (position[0]+25,position[1]-15), (position[0]+25,position[1]), (position[0]+75,position[1]), (position[0]+75,position[1]+20), (position[0]+25,position[1]+20), (position[0]+25,position[1]+35)], fill=colour)
        ImageDraw.Draw(im).text((position[0]+20,position[1]+40), label, fill="black")
    else: print("Invalid direction: Specify either \"+\" or \"-\".")
    del draw


def get_N_HexCol(num_cols):
    HSV_tuples = [(x * 1.0 / num_cols, 0.5, 0.5) for x in range(num_cols)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#%02x%02x%02x' % tuple(rgb))
    return hex_out


def drawCluster(im, strain, genes, colours, position):
    #im = Image.new("RGBA", (100*len(genes)+125,100), (255,255,255,255))
    #im2 = Image.resize((100*len(genes)+100,200))
    sorted_genes = sorted(genes.strip().split(" "), key=lambda x: x.split(":")[1])
    sameContig = [list(v) for k, v in itertools.groupby(sorted_genes, lambda y: y.split(":")[1])]
    groupPos = 500
    fnt = ImageFont.truetype('AcademicM54.ttf', 60)
    ImageDraw.Draw(im).text((80,400*position-20), strain, font=fnt, fill="black")
    for group in sameContig:
        sorted_group = sorted(group, key=lambda x: int(x.split(":")[3]))
        groupL = 0
        prevEnd = None
        gapList = [0]
        for gene in sorted_group:
            splt = gene.split(":")
            gap = 0
            if prevEnd:
                gap = int(splt[3])-prevEnd
                gapList.append(gap)
                if gap > 10000:
                    gap = 10000
            prevEnd = int(splt[4])
            groupL += int(splt[5]) + gap
        group_start = int(sorted_group[0].split(":")[3])
        group_end = int(sorted_group[-1].split(":")[4])
        group_length = groupL/25+200
        #ImageDraw.Draw(im).line((groupPos, 100*position, 100*len(group)+groupPos+25, 100*position), fill="black")
        #ImageDraw.Draw(im).line((groupPos, 100*position, group_length+groupPos+25, 100*position), fill="black")
        ImageDraw.Draw(im).line((groupPos, 400*position, group_length+groupPos+len(group)*100, 400*position), fill="black", width=11)
        idx = 0
        #colours = get_N_HexCol(len(genes))
        genePos = 0
        prevPos = group_start
        prevLen = 0
        for gene in sorted_group:
            geneSPLT = gene.split(":")
            dir = geneSPLT[2]
            lab = geneSPLT[0]
            col = colours[lab]
            length = int(geneSPLT[5])/25
            #pos = (groupPos+25+100*(idx), 100*position-10)
            if gapList[idx] > 10000:
                genePos += prevLen +100 + 400
                ImageDraw.Draw(im).polygon([(groupPos+genePos-300, 400*position-40), (groupPos+genePos-100, 400*position-40), (groupPos+genePos-100, 400*position+40), (groupPos+genePos-300, 400*position+40)], fill="red")
                #ImageDraw.Draw(im).line((groupPos+genePos-30, 100*position-20, groupPos+genePos-30, 100*position+20), fill="red")
                ImageDraw.Draw(im).text((groupPos+genePos-280, 400*position-20), str(gapList[idx]), fill="black")
                #ImageDraw.Draw(im).text((groupPos+genePos-40, 100*position-30), str(int(geneSPLT[3])-prevPos), fill="black")
                #genePos += length + 75
            #elif 0 < (int(geneSPLT[3])-group_start)/100-prevPos < 35:
            #    genePos += length + 35
            #elif gene != sorted_group[-1]:
            #    genePos += length+gap/100 +25 +gapList[idx]/100
            else:
                genePos += prevLen+gapList[idx]/25 +100
            pos = (groupPos+genePos, 400*position-40)
            drawArrow(im, dir, col, pos, lab, length)
            prevPos = int(geneSPLT[4])
            prevLen = length
            idx += 1
        groupPos += group_length+len(group)*100+100
    #im.save("test_img.png", "PNG")

### old version
def drawCluster2(im, strain, genes, colours, position):
    #im = Image.new("RGBA", (100*len(genes)+125,100), (255,255,255,255))
    #im2 = Image.resize((100*len(genes)+100,200))
    sorted_genes = sorted(genes.strip().split(" "), key=lambda x: x.split(":")[1])
    sameContig = [list(v) for k, v in itertools.groupby(sorted_genes, lambda y: y.split(":")[1])]
    groupPos = 125
    for group in sameContig:
        ImageDraw.Draw(im).line((groupPos, 100*position, 100*len(group)+groupPos+25, 100*position), fill="black")
        idx = 0
        #colours = get_N_HexCol(len(genes))
        sorted_group = sorted(group, key=lambda x: int(x.split(":")[3]))
        for gene in sorted_group:
            fnt = ImageFont.truetype('Pillow/Tests/fonts/FreeMono.ttf', 15)
            ImageDraw.Draw(im).text((20,100*position-5), strain, font=fnt, fill="black")
            geneSPLT = gene.split(":")
            dir = geneSPLT[2]
            lab = geneSPLT[0]
            col = colours[lab]
            pos = (groupPos+25+100*(idx), 100*position-10)
            idx += 1
            drawArrow(im, dir, col, pos, lab)
        groupPos += 100*len(group)+50
    #im.save("test_img.png", "PNG")

def drawMultipleClusters(clusters):
    with open(clusters, "r") as file:
        numGenes = []
        numClusts = 0
        listGenes = {}
        for line in file:
            if line.startswith("Strain") or line.split("\t")[0] in omitStrains:
                continue
            else:
                numClusts += 1
                numGenes.append(int(line.split("\t")[3]))
                for gene in line.split("\t")[4].split(" "):
                    listGenes[gene.split(":")[0]] = None
        maxGenes = max(numGenes)
        #print(maxGenes)
        #print(numClusts)
        colours = get_N_HexCol(len(listGenes))
        idx = 0
        for gene in listGenes:
            listGenes[gene] = colours[idx]
            idx += 1
        #print(listGenes)
        im = Image.new("RGBA", (560*maxGenes+400,440*numClusts+400), (255,255,255,255))
        clustIdx = 1
        file.seek(0)
        for line in file:
            if line.startswith("Strain") or line.split("\t")[0] in omitStrains:
                continue
            else:
                splt = line.split("\t")
                strain = splt[0]
                cluster = splt[4]
                drawCluster(im, strain, cluster, listGenes, clustIdx)
                clustIdx += 1
        outfile = outputFile + ".png"
        outfile2 = outputFile + ".eps"
        im.save(outfile, format = "PNG", quality = 95)
        
        #png = Image.open(outfile)
        #png.load()
        '''
        eps = Image.new("RGB", im.size, (255, 255, 255))
        eps.paste(im, mask=im.split()[3])
        
        eps.save(outfile2, format = "EPS", quality = 80)
        '''
        #im.save(outfile2, "PDF", save_all=True)


drawMultipleClusters(inputFile)

#drawCluster([["cyp1","l"],["cyp4","r"],["cyp2","r"],["cyp5","l"],["cyp3","r"]])

#draw = ImageDraw.Draw(im)
#draw.line((0, 0) + im.size, fill=128)
#draw.line((0, im.size[1], im.size[0], 0), fill=128)
#del draw

# write to stdout
#im.save("test_img.png", "PNG")
