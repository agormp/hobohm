#!/usr/bin/env python3
import sys, os.path
from optparse import OptionParser

################################################################################################

def main():
    parser = build_parser()
    options = parse_commandline(parser)
    neighbordict = build_neighbordict(options)
    if options.keepfile:
        neighbordict = handle_keeplist(options, neighbordict)
    neighbor_count_dict = remove_neighbors(neighbordict)
    keepnames = reinstate_neighborless_skipped(neighbordict, neighbor_count_dict)
    print_keepnames(keepnames)

################################################################################################

def build_parser():
    parser = OptionParser(usage="usage: hobohm [-s|-d] FILE -c CUTOFF [-k KEEPFILE]",
                          version="0.0.1")

    parser.add_option("-s", type="string", dest="simfile", metavar="SIMFILE",
                          help="file with pairwise similarities: name1 name2 sim")

    parser.add_option("-d", type="string", dest="distfile", metavar="DISTFILE",
                          help="file with pairwise distances: name1 name2 dist")

    parser.add_option("-c", type="float", dest="cutoff", metavar="CUTOFF",
                          help="cutoff for deciding which pairs are neighbors")

    parser.add_option("-k", type="string", dest="keepfile", metavar="KEEPFILE",
                          help="file with names that must be kept (one name per line)")

    parser.set_defaults(simfile=None, distfile=None, cutoff=None, keepfile=None)
    return parser

################################################################################################

def parse_commandline(parser):

    (options, args) = parser.parse_args()
    if ((options.simfile and options.distfile) or
       ((options.simfile is None) and (options.distfile is None))):
        parser.error("Use either option -s (similarity) or option -d (distance)")
    if options.cutoff is None:
        parser.error("Must provide cutoff (option -c)")
    return(options)

################################################################################################

def build_neighbordict(options):
    neighbordict = {}

    # Could use operator module and assign < or > to variable to avoid repeated code
    # but is less readable...
    if options.simfile:
        with open(options.simfile, "r") as infile:
            for line in infile:
                name1,name2,similarity=line.split()
                if name1 not in neighbordict:
                    neighbordict[name1]=set()
                if name2 not in neighbordict:
                    neighbordict[name2]=set()
                if float(similarity) > options.cutoff and name1 != name2:
                    neighbordict[name1].add(name2)
                    neighbordict[name2].add(name1)
    elif options.distfile:
        with open(options.distfile, "r") as infile:
            for line in infile:
                name1,name2,distance=line.split()
                if name1 not in neighbordict:
                    neighbordict[name1]=set()
                if name2 not in neighbordict:
                    neighbordict[name2]=set()
                if float(distance) < options.cutoff and name1 != name2:
                    neighbordict[name1].add(name2)
                    neighbordict[name2].add(name1)

    return neighbordict

################################################################################################

def handle_keeplist(options, neighbordict):

    # Read list of names to keep
    keepset = set()
    with open(options.keepfile, "r") as infile:
        for line in infile:
            words = line.split()
            keepset.add(words[0])

    # First, check if any pair of keepset members are neighbors.
    # If so, print warning and artificially hide this fact
    for keepname1 in keepset:
        for keepname2 in (keepset - set([keepname1])):
            if keepname2 in neighbordict[keepname1]:
                print("# Keeplist warning: {} and {} are neighbors!".format(keepname1, keepname2))
                neighbordict[keepname1].remove(keepname2)
                neighbordict[keepname2].remove(keepname1)

    # Then, remove all neighbors of keepset members
    for keepname in keepset:
        keepname_neighbors = list(neighbordict[keepname]) # To avoid issues with changing set during iteration

        for remseq in keepname_neighbors:
            # Remove any mention of remseq in other entries
            remseq_neighbors = list(neighbordict[remseq])
            for remseq_neighbor in remseq_neighbors:
                neighbordict[remseq_neighbor].remove(remseq)

            # Now remove remseq entry itself
            del neighbordict[remseq]

    return(neighbordict)

################################################################################################

def remove_neighbors(neighbordict):

    # Build dictionary keeping track of how many neighbors each item has
    neighbor_count_dict = {}
    for name in neighbordict:
        neighbor_count_dict[name]=len(neighbordict[name])

    # Find max number of neighbors
    maxneighb = max(neighbor_count_dict.values())

    # While some items still have neighbors: remove item with most neighbors, update counts
    # Note: could ties be dealt with intelligently?
    while maxneighb > 0:

        # Find an item that has maxneighb neighbors, and remove it from list
        for remove_name, count in neighbor_count_dict.items():
            if count == maxneighb:
                break
        del(neighbor_count_dict[remove_name])

        # Update neighbor counts
        for neighbor in neighbordict[remove_name]:
            if neighbor in neighbor_count_dict:
                neighbor_count_dict[neighbor] -= 1

        # Find new maximum number of neighbors
        maxneighb = max(neighbor_count_dict.values())

    return neighbor_count_dict

################################################################################################

def reinstate_neighborless_skipped(neighbordict, neighbor_count_dict):

    # Postprocess: reinstate skipped sequences that now have no neighbors
    # (This can happpen when ...?)
    # Note: order may have effect. Could this be optimized?

    allseqs=set(neighbordict.keys())
    keepseqs=set(neighbor_count_dict.keys())
    skipseqs=allseqs - keepseqs

    for skipped in skipseqs:
        # if skipped sequence has no neighbors in keeplist
        if not (neighbordict[skipped] & keepseqs):
            keepseqs.add(skipped)

    return keepseqs

################################################################################################

def print_keepnames(keepnames):
    # Print retained items
    for name in keepnames:
        print(name)

################################################################################################

if __name__ == "__main__":
    main()
