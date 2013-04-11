#!/usr/bin/env python
#
# Author: Yu Fu
# python my_alignment.py [query.seq] [template.seq]
# gapo, a non-negative number, is the NEGATIVE of gap opening penalty
# gape, a non-negative number, is the NEGATIVE of gap extension penalty
#

import sys
import re

if len(sys.argv) <3:
    print "Usage: python %s query.seq template.seq" % sys.argv[0]

gapo = 10
gape = 1
try:
    mf = open("blosum62.txt", "r")
except:
    print "Cannot open the scoring matrix file."
    sys.exit(1)

aas = mf.readline()
aas = aas.strip()
aas = re.split(re.compile(r"\s+"), aas)
s = {}
for l in mf.readlines():
    l = l.strip()
    v = re.split(re.compile(r"\s+"), l)
    cur_aa = v[0]
    v = [int(x) for x in v[1:]]
    s[cur_aa] = dict(zip(aas, v))
mf.close()

# read the query sequence
query = ""
qf = open(sys.argv[1], "r")
for line in qf.readlines():
    line = line.strip()
    query += line
qf.close()

# Read the template sequence
template = ""
tf = open(sys.argv[2], "r")
for line in tf:
    line = line.strip()
    template += line
tf.close()

# Global alignment
# Initialize all the dictionaries (as Python sadly does not support
# autovivification):
m = {}
pi ={}
pj = {}

for i in range(0, len(query)):
    m[i] = {}
    pi[i] = {}
    pj[i] = {}
    for j in range(0, len(template)):
        m[i][j] = 0
        pi[i][j] = 0
        pj[i][j] = 0

for i in range(0, len(query)):
    m[i][0] = s[query[i]][template[0]]
    pi[i][0] = 0
    pj[i][0] =0
for j in range(0, len(template)):
    m[0][j] = s[query[0]][template[j]]
    pi[0][j] = 0
    pj[0][j] = 0

for i in range(1, len(query)):
    for j in range(1, len(template)):
        # Match/Mismatch
        maxm = m[i-1][j-1] + s[query[i]][template[j]]
        maxi = i - 1
        maxj = j - 1
        # Insert gap(s) in template sequence
        for k in range(0, i-1):
            max_i = m[k][j-1] + s[query[i]][template[j]]
            max_i -= gapo + gape * (i - k - 2)
            if max_i > maxm:
                maxm = max_i
                maxi = k
                maxj = j - 1
        # Insert gap(s) in the query sequnce
        for k in range(0, j-1):
            max_j = m[i-1][k] + s[query[i]][template[j]]
            max_j -= gapo + gape * (j - k - 2)
            if max_j > maxm:
                maxm = max_j
                maxi = i - 1
                maxj = k
        # Store the best score for the current position
        m[i][j] = maxm
        pi[i][j] = maxi
        pj[i][j] = maxj

# Output to see if these three matrices contain the correct values        
# for i in range(0, len(m)):
#     for j in range(0, len(m[i])):
#         print "%2d(%1d,%1d)" % (m[i][j], pi[i][j], pj[i][j]),
#     print 

# Identify the maximum score
maxm = m[len(query)-1][len(template)-1]
pii = len(query) - 1
pjj = len(template) - 1
for i in range(0, len(query)):
    if m[i][len(template)-1] > maxm:
        maxm = m[i][len(template)-1]
        pii = i
        pjj = len(template) - 1
for j in range(0, len(template)):
    if m[len(query)-1][j] > maxm:
        maxm = m[len(query)-1][j]
        pii = len(query) - 1
        pjj = j

print "Max score: %d" % maxm

# trace back
ti = 0
tj = 0
aln1 = ''
aln2 = ''
# Add the '-'s at the very end of the sequence
for i in range(len(query)-1, pii, -1):
    aln1 += query[i]
    aln2 += '-'
for j in range(len(template)-1, pjj, -1):
    aln1 += '-'
    aln2 += template[j]

aln1 += query[pii]
aln2 += template[pjj]

while pii>0 and pjj>0:
    for i in range(pii-1, pi[pii][pjj], -1):
        aln1 += query[i]
        aln2 += '-'
    for j in range(pjj-1, pj[pii][pjj], -1):
        aln1 += '-'
        aln2 += template[j]
    pii2 = pi[pii][pjj]
    pjj2 = pj[pii][pjj]
    pii = pii2
    pjj = pjj2
    aln1 += query[pii]
    aln2 += template[pjj]

# Add the '-'s at the other end of the sequence
for i in range(pii-1, 0-1, -1):
    aln1 += query[i]
    aln2 += "-"
for j in range(pjj-1, 0-1, -1):
    aln1 += "-"
    aln2 += template[j]
try:
    out = open("alignment.txt", "w")
except:
    print "I cannnot write to 'alignment.txt'"

print aln1[::-1]
print >>out, aln1[::-1]
print aln2[::-1]
print >>out, aln2[::-1]

out.close()