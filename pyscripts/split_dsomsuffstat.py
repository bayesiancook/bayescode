import sys
import random

inname = sys.argv[1]
nrep = int(sys.argv[2])
outname = sys.argv[3]

ds = dict()
dsnorm = dict()
dn = dict()
dnnorm = dict()

genelist = []
ngene = 0

with open(inname, 'r') as infile:
    line = infile.readline().rstrip('\n')
    ngene = int(line)
    line = infile.readline().rstrip('\n')
    while line:
        gene = line
        genelist.append(gene)
        ds[gene] = infile.readline()
        dsnorm[gene] = infile.readline()
        dn[gene] = infile.readline()
        dnnorm[gene] = infile.readline()
        line = infile.readline().rstrip('\n')

random.shuffle(genelist)

if ngene != len(genelist):
    print("error: non matching number of genes")
    sys.exit(1)

npersplit = ngene // nrep

for rep in range(nrep):
    with open("{0}{1}.genedsomsuffstat".format(outname, rep), 'w') as outfile:
        outfile.write("{0}\n".format(npersplit))
        for i in range(npersplit):
            gene = genelist[npersplit*rep + i]
            outfile.write("{0}\n{1}{2}{3}{4}".format(gene, ds[gene], dsnorm[gene], dn[gene], dnnorm[gene]))






