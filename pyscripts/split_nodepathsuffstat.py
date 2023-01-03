import sys
import random

inname = sys.argv[1]
nrep = int(sys.argv[2])
outname = sys.argv[3]

gene2pathss = dict()

genelist = []
ngene = 0

with open(inname, 'r') as infile:
    ngene = int(infile.readline().rstrip('\n'))
    for line in infile:
        gene = line.rstrip('\n').split()[0]
        genelist.append(gene)
        gene2pathss[gene] = line

if ngene != len(genelist):
    print("error: non matching number of genes")
    sys.exit(1)

random.shuffle(genelist)

npersplit = ngene // nrep

for rep in range(nrep):
    with open("{0}{1}.genenodepathsuffstat".format(outname, rep), 'w') as outfile:
        outfile.write("{0}\n".format(npersplit))
        for i in range(npersplit):
            gene = genelist[npersplit*rep + i]
            outfile.write(gene2pathss[gene])






