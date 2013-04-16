#!/usr/bin/env python
#
# Author: Yu Fu
# python structure_alignment.py [query.seq] [template.seq]
# gapo, a non-negative number, is the NEGATIVE of gap opening penalty
# gape, a non-negative number, is the NEGATIVE of gap extension penalty
#

import sys
import re
import math

def pre_calc(matrix_file):
    """Take a file containing the matrix of coordinates of all the AAs,
    and return the distance and secondar structure between each other
    
    Arguments:
    - `matrix_file`:
    """
    fh = open(matrix_file, "r")
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

    # for d in dis:
    #     print d
    return dis, structure

def report_indel(s1, s2, dis, structure):
    """Report the insertions and deletions.
    s1 is the query sequence and s2 is the template
    """
    # One need to see the letter around the gap and r defines the number
    # thereof.
    s_len = len(s1)
    i = []
    p = re.compile(r"-+")
    for f in p.finditer(s1):
        i.append([f.start(), f.end()])
    j = []
    for f in p.finditer(s2):
        j.append([f.start(), f.end()])
    print "Deletion(s):"
    for ii in i:
        # See if it is structurally compatible
        if ii[0] == 0:
            print "Skipping the gaps at the beginning..."
            print "#" * 80
            continue
        if ii[1] == len(s1):
            print "Skipping the gaps at the end..."
            print "#" * 80            
            continue
        ori_left_str = s2[0:ii[0]].replace("-", "")
        ori_left = len(ori_left_str) - 1
        ori_right_str = s2[0:ii[1]].replace("-", "")
        ori_right = len(ori_right_str)
        print "Locations of the ends: %d, %d" % (ori_left, ori_right)

        if dis[ori_left][ori_right] <= 5.0:
            print "This DELETION is structurally COMPATIBLE."
        else:
            print "This DELETION is structurally INCOMPATIBLE."
        print "The 3D distance of of the two ends is %f." % \
            dis[ori_left][ori_right]
        print ">query"
        print s1[ii[0]:ii[1]]
        print ">Template"
        print s2[ii[0]:ii[1]]
        print "#" * 80
    print "Insertion(s):"
    for jj in j:
        # See if it is structurally compatible (keep in mind that the template
        # sequence already contains gaps)
        # Query ACTG---TCGA
        # Templ --AGACDACT-
        if jj[0] == 0:
            print "Skip the gaps at the beginning..."
            print "#" * 80
            continue
        if jj[1] == len(s1):
            print "Skipping the gaps at the end..."
            print "#" * 80
            continue
        # From the very beginning to the the start of the gap
        ori_left_str = s2[0:jj[0]].replace("-", "")
        ori_left = len(ori_left_str) - 1
        # From the very beginning to the end of the gap
        ori_right_str = s2[0:jj[1]].replace("-", "")
        ori_right = len(ori_right_str)
        print "Locations of the ends: %d, %d" % (ori_left, ori_right)

        
        if structure[ori_left] == "C" or structure[ori_right] == "C":
            print "This INSERTION is structurally COMPATIBLE."
        else:
            print "This INSERTION is structurally INCOMPATIBLE."
        print "The secondary structures of both ends are %s, %s" % \
            (structure[ori_left], structure[ori_right])
            # Recover the original location
        print ">Query"
        print s1[jj[0]:jj[1]]
        print ">Template"
        print s2[jj[0]:jj[1]]
        print "#" * 80

def ins_is_compatible(start, end, structure):
    if structure[start] == "C" or structure[end] == "C":
        return True
    else:
        return False
def del_is_compatible(start, end, dis):
    if dis[start][end] <= 5.0:
        return True
    else:
        return False

if __name__ == "__main__":
    if len(sys.argv) <3:
        print "Usage: python %s query.seq template.seq" % sys.argv[0]
    gapo = 10
    gape = 1
    template_xyz = "template.xyz"
    # template_xyz = "my_temp.xyz"
    # Additional penalty for structurally incompatible gaps
    gapi = 8

    # Change this option to False if you don't want structual alignment
    check_compatibility = True
    try:
        mf = open("blosum62.txt", "r")
    except:
        print "Cannot open the scoring matrix file."
        sys.exit(1)

    # Get the distances between any pairs of AAs
    dis, structure = pre_calc(template_xyz)
    # print "dis[11][12] is " + str(dis[11][12])
    # print dis
    
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
            # Insert gap(s) in template sequence (INSERTION)
            for k in range(0, i-1):
                max_i = m[k][j-1] + s[query[i]][template[j]]
                max_i -= gapo + gape * (i - k - 2)

                # If this option of "check_compatibility" is turned on...
                if check_compatibility == True:
                    # If this insertion is incompatible...
                    if ins_is_compatible(k, j-1, structure) == False:
                        max_i -= gapi
                
                if max_i > maxm:
                    maxm = max_i
                    maxi = k
                    maxj = j - 1
            # Insert gap(s) in the query sequnce (DELETION)
            for k in range(0, j-1):
                max_j = m[i-1][k] + s[query[i]][template[j]]
                max_j -= gapo + gape * (j - k - 2)
                # If this option of "check_compatibility" is turned on...
                if check_compatibility == True:
                    if del_is_compatible(i-1, k, dis) == False:
                        max_j -= gapi

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
        out = open("alignment.fasta", "w")
    except:
        print "I cannnot write to 'alignment.txt'"

    # Output with a fixed width to STDOUT and output as is in the file
    width = 60
    print ">Query"
    reversed_aln1 = aln1[::-1]
    length = len(reversed_aln1)
    for i in range(0, length):
        if i % width == 0 and i != 0:
            print
        sys.stdout.write(reversed_aln1[i])
    print

    print >>out, ">Query"
    print >>out, aln1[::-1]

    print ">Template"
    reversed_aln2 = aln2[::-1]
    length = len(reversed_aln2)
    for i in range(0, length):
        if i % width == 0 and i != 0:
            print
        sys.stdout.write(reversed_aln2[i])
    print

    print >>out, ">Template"
    print >>out, aln2[::-1]
    print # An empty line

    out.close()
    report_indel(reversed_aln1, reversed_aln2, dis, structure)

        
        
            
            
                
            
    
    




