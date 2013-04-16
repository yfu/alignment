#!/usr/bin/env python
#
# Pre-calculate the distances between any pairs of AA and also reads the
# secondary structure for each AA
import re
import math
import sys

fh = open(sys.argv[1], "r")
coor = []
structure = ""
seq = ""
pat_comment = re.compile(r"^#")
for line in fh.readlines():
    line = line.strip()
    if pat_comment.search(line) != None:
        continue
    cols = line.split()
    seq += cols[0]
    structure += cols[1]
    coor.append((float(cols[2]), float(cols[3]), float(cols[4])))

# print seq
# print structure
# for c in coor:
#     print c

# dis(tances):
#    A  T  C  T
# A  0  1  4  8
# T  1  3  5  10
# C  4  5  12 15
# G  8  10 15 20 

dis = []
for i in range(0, len(seq)):
    dis.append([])
    for j in range(0, len(seq)):
        dis[i].append(
            # Calculate the distance
            math.sqrt(sum(map(lambda x1, x2: (x1 - x2) * (x1 - x2),
                              coor[i], coor[j]))))

print dis        