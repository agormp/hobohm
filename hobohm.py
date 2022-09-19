#!/usr/bin/env python3
import argparse, sys
from collections import defaultdict
from operator import itemgetter

################################################################################################

def main():
    parser = build_parser()
    args = parse_commandline(parser)
    neighbors = parse_neighbor_info(args)
    if args.keepfile:
        neighbors = handle_keeplist(args, neighbors)
    neighbor_count = remove_neighbors(neighbors)
    keepnames = reinstate_neighborless_skipped(neighbors, neighbor_count)
    print_keepnames(keepnames)

################################################################################################

def build_parser():
    parser = argparse.ArgumentParser(description = "Selects representative subset of data based on" +
                                    " list of pairwise similarities (or distances), such that no" +
                                    " retained items are close neighbors")

    parser.add_argument("pairfile", metavar='PAIRFILE', default="-",
                        help="file containing the similarity (option -s) or distance " +
                             "(option -d) for each pair of items: name1 name2 value")

    distsimgroup = parser.add_mutually_exclusive_group()
    distsimgroup.add_argument('-s', action='store_true', dest="values_are_sim",
                              help="values in PAIRFILE are similarities " +
                                    "(larger values = more similar)")
    distsimgroup.add_argument('-d', action='store_true', dest="values_are_dist",
                              help="values in PAIRFILE are distances " +
                                    "(smaller values = more similar)")

    parser.add_argument("-c",  action="store", type=float, dest="cutoff", metavar="CUTOFF",
                          help="cutoff for deciding which pairs are neighbors")

    parser.add_argument("-k", action="store", dest="keepfile", metavar="KEEPFILE",
                          help="file with names of items that must be kept (one name per line)")

    return parser

################################################################################################

def parse_commandline(parser):

    args = parser.parse_args()
    if ((args.values_are_sim and args.values_are_dist) or
       ((not args.values_are_sim ) and (not args.values_are_dist))):
        parser.error("Must specify either option -s (similarity) or option -d (distance)")
    if args.cutoff is None:
        parser.error("Must provide cutoff (option -c)")
    return(args)

################################################################################################

def parse_neighbor_info(args):

    # neighbors should, for each item, contain a set of its neighbors (possibly empty)
    neighbors = defaultdict(set)

    cutoff = args.cutoff    # Micro optimization: save time looking up dotted attributes
    values_are_sim = args.values_are_sim

    with open(args.pairfile, "r") as infile:

        for line in infile:
            name1,name2,value = line.split()

            if name1 != name2:
                value = float(value)

                if values_are_sim:
                    if value > cutoff:
                        neighbors[name1].add(name2)
                        neighbors[name2].add(name1)
                else:
                    if value < cutoff:
                        neighbors[name1].add(name2)
                        neighbors[name2].add(name1)

    return neighbors

################################################################################################

def handle_keeplist(args, neighbors):

    # Read list of names to keep
    keepset = set()
    with open(args.keepfile, "r") as infile:
        for line in infile:
            word = line.rstrip()
            if word:                    # Maybe overkill to allow blank lines...
                keepset.add(word)

    # First, check if any pair of keepset members are neighbors.
    # If so, print warning on stderr and artificially hide this fact
    for keepname1 in keepset:
        for keepname2 in (keepset - set([keepname1])):
            if keepname2 in neighbors[keepname1]:
                sys.stderr.write("# Keeplist warning: {} and {} are neighbors!".format(keepname1, keepname2))
                neighbors[keepname1].remove(keepname2)
                neighbors[keepname2].remove(keepname1)

    # Then, remove all neighbors of keepset members
    for keepname in keepset:
        keepname_neighbors = list(neighbors[keepname]) # To avoid issues with changing set during iteration

        for remseq in keepname_neighbors:
            # Remove any mention of remseq in other entries
            remseq_neighbors = list(neighbors[remseq])
            for remseq_neighbor in remseq_neighbors:
                neighbors[remseq_neighbor].remove(remseq)

            # Now remove remseq entry itself
            del neighbors[remseq]

    return(neighbors)

################################################################################################

def remove_neighbors(neighbors):

    # Build dictionary keeping track of how many neighbors each item has
    neighbor_count = {}
    for name in neighbors:
        neighbor_count[name] = len(neighbors[name])

    # Find max number of neighbors and corresponding name
    item_with_most_nb, max_num_nb = max(neighbor_count.items(), key=itemgetter(1))

    # While some items still have neighbors: remove an item with most neighbors, update counts
    # Note: could ties be dealt with intelligently?
    while max_num_nb > 0:

        del(neighbor_count[item_with_most_nb])

        # Update neighbor counts
        for item in neighbors[item_with_most_nb]:
            if item in neighbor_count:
                neighbor_count[item] -= 1

        # Find new maximum number of neighbors
        item_with_most_nb, max_num_nb = max(neighbor_count.items(), key=itemgetter(1))

    return neighbor_count

################################################################################################

def reinstate_neighborless_skipped(neighbors, neighbor_count):

    # Postprocess: reinstate skipped sequences that now have no neighbors
    # (This can happpen when ...?)
    # Note: order may have effect. Could this be optimized?

    allseqs=set(neighbors.keys())
    keepseqs=set(neighbor_count.keys())
    skipseqs=allseqs - keepseqs

    for skipped in skipseqs:
        # if skipped sequence has no neighbors in keeplist
        if not (neighbors[skipped] & keepseqs):
            keepseqs.add(skipped)

    return keepseqs

################################################################################################

def print_keepnames(keepnames):
    for name in keepnames:
        print(name)

################################################################################################

if __name__ == "__main__":
    main()
