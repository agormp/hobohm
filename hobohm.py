#!/usr/bin/env python3
import sys, os.path

# Parse commandline, check that file name and cutoff is given and that file exists
# Note: should add support for stdin!!!
if not sys.argv[2:]:
    print("usage: %s CUT-OFF NEIGHBOR-LIST [INCLUDE-LIST]\n" % sys.argv[0])
    print("Hobohm2 algorithm for optimal homology-reduction\n")
    print("Input: a list of sequence similarities and a cutoff value.")
    print("       Sequences that are more similar than the cutoff, are")
    print("       considered to be neighbors.\n")
    print("       It is possible to also provide a list of genes")
    print("       that have to be included in the list\n")
    print("Output: a list of sequences that should be kept in the")
    print("       homology-reduced data set\n")
    print("Format of NEIGHBOR-LIST:")
    print("       seqid1 seqid2 similarity")
    print("       seqid1 seqid2 similarity")
    print("       ...\n")
    print("Format of KEEP-LIST:")
    print("       seqid")
    print("       seqid")
    print("       ...\n")
    sys.exit()

cutoff = float(sys.argv[1])

if not os.path.isfile(sys.argv[2]):
    print("%s: file %s does not exist\n" % (sys.argv[0], sys.argv[2]))
    sys.exit()

#######################################################################################
# Open neighbor-file, iterate over it one line at a time
# Read list of similarities, build list of names and keep track of neighbors.
# Neighbor info kept in dictionary, where key is seqname and value is Set of neighbors.
neighborfile = open(sys.argv[2], "r")
neighborlist = {}

# Each line in file has format: seqid1 seqid2 similarity
for simline in neighborfile:

    words=simline.split()
    if len(words)==3:                # Sanity check: do lines conform to expected format?
        seq1=words[0]
        seq2=words[1]
        similarity=float(words[2])

    # Add sequence names as we go along
    if seq1 not in neighborlist:
        neighborlist[seq1]=set()
    if seq2 not in neighborlist:
        neighborlist[seq2]=set()

    # Build lists of neighbors as we go along.
    # Note: Set.add() method automatically enforces member uniqueness - saves expensive test!
    if similarity > cutoff and seq1 != seq2:
        neighborlist[seq1].add(seq2)
        neighborlist[seq2].add(seq1)

neighborfile.close()

# Read includelist if provided
includes = False
if sys.argv[3:]:
    includes = True
    includelist = set()
    includefile = open(sys.argv[3], "r")
    for line in includefile:
            words = line.split()
            includelist.add(words[0])
    includefile.close()

########################################################################################

# If includelist was provided: preprocess list
if includes:

    # (1) Check if any pair of includelist members are neighbors.
    #     If so, print warning and artificially hide this fact
    for incseq1 in includelist:
        for incseq2 in (includelist - set([incseq1])):
            if incseq2 in neighborlist[incseq1]:
                print("# Includelist warning: %s and %s are neighbors!" % (incseq1, incseq2))
                neighborlist[incseq1].remove(incseq2)
                neighborlist[incseq2].remove(incseq1)

    # (2) Remove all neighbors of includelist members
    for incseq in includelist:
        incseq_neighbors = list(neighborlist[incseq]) # To avoid issues with changing set during iteration

        for remseq in incseq_neighbors:

            # Remove any mention of remseq in other entries
            remseq_neighbors = list(neighborlist[remseq])
            for remseq_neighbor in remseq_neighbors:
                neighborlist[remseq_neighbor].remove(remseq)

            # Now remove remseq entry itself
            del neighborlist[remseq]

########################################################################################

# Build dictionary keeping track of how many neighbors each sequence has
nr_dict = {}
for seq in list(neighborlist.keys()):
    nr_dict[seq]=len(neighborlist[seq])

# Find max number of neighbors
maxneighb = max(nr_dict.values())

# While some sequences in list still have neighbors: remove the one with most neighbors, update counts
# Note: could ties be dealt with intelligently?
while maxneighb > 0:

    # Find an entry that has maxneighb neighbors, and remove it from list
    for remove_seq in list(nr_dict.keys()):
        if nr_dict[remove_seq] == maxneighb: break
    del(nr_dict[remove_seq])

    # Update neighbor counts
    for neighbor in neighborlist[remove_seq]:
        if neighbor in nr_dict:
            nr_dict[neighbor] -= 1

    # Find new maximum number of neighbors
    maxneighb = max(nr_dict.values())

##############################################################################################
# Postprocess: reinstate skipped sequences that now have no neighbors
# Note: order may have effect. Could this be optimized?

allseqs=set(neighborlist.keys())
keepseqs=set(nr_dict.keys())
skipseqs=allseqs - keepseqs

for skipped in skipseqs:
    # if skipped sequence has no neighbors in keeplist
    if not (neighborlist[skipped] & keepseqs):
        keepseqs.add(skipped)

# Print remaining sequences
for seq in keepseqs:
    print(seq)

