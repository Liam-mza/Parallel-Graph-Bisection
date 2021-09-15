#!/usr/bin/env python2.7
import sys
import random

if (len(sys.argv) != 4):
    print >> sys.stderr, "Usage: ./create-graph [num-vertices] [num-edges] [file-name]"
    exit(1)

nVtx = int(sys.argv[1])
nEdge = int(sys.argv[2])
fileName = sys.argv[3]

edgePerVtx = nEdge / nVtx
remEdge = nEdge % nVtx

edges = dict()
for i in range(0, nVtx):
    j = 0
    while j < edgePerVtx:
        edge = random.randrange(0, nVtx)
        if ((i, edge) in edges or (edge, i) in edges):
            continue
        else:
            edges[(i, edge)] = 1
            j = j + 1

for i in range(0, remEdge):
    l = random.randrange(0, nVtx);
    r = random.randrange(0, nVtx);
    while l == r or (l,r) in edges or (r,l) in edges:
        l = random.randrange(0, nVtx)
        r = random.randrange(0, nVtx)
    edges[(l, r)] = 1


f = open(fileName, 'w')
f.write('%d %d\n' % (nVtx, nEdge))
for i in edges:
    f.write('%d %d\n' % (i[0], i[1]))

f.close()

print 'A graph with %d vertices and %d edges is written into the file %s.' % (nVtx, nEdge, fileName)
