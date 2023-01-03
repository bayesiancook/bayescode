import sys

innames = sys.argv[1:len(sys.argv)-1]
outname = sys.argv[len(sys.argv)-1]

ds = dict()
dsnorm = dict()
dn = dict()
dnnorm = dict()

genelist = []
totngene = 0

for inname in innames:
    with open(inname, 'r') as infile:
        line = infile.readline().rstrip('\n')
        ngene = int(line)
        totngene += ngene
        for i in range(ngene):
            gene = infile.readline().rstrip('\n')
            genelist.append(gene)
            ds[gene] = infile.readline()
            dsnorm[gene] = infile.readline()
            dn[gene] = infile.readline()
            dnnorm[gene] = infile.readline()

with open(outname, 'w') as outfile:
    outfile.write("{0}\n".format(totngene))
    for gene in genelist:
        outfile.write("{0}\n{1}{2}{3}{4}".format(gene, ds[gene], dsnorm[gene], dn[gene], dnnorm[gene]))

