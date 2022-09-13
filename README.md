# hobohm: command line program for selecting representative, non-redundant data set from larger set, based on list of pairwise similarities (or distances).

[![PyPI downloads](https://static.pepy.tech/personalized-badge/hobohm?period=total&units=none&left_color=black&right_color=blue&left_text=downloads&service=github)](https://pepy.tech/project/hobohm)
![](https://img.shields.io/badge/version-1.0.0-blue)

The "Hobohm" algorithm was originally created with the purpose of selecting representative, non-redundant sets of protein data from a larger data set. Non-redundant here means that the resulting data set should contain no pairs of sequences with high similarity:

["Selection of representative protein data sets", Protein Sci. 1992. 1(3):409-17](https://pubmed.ncbi.nlm.nih.gov/1304348/).

This command-line program implements algorithm 2 from that paper.

The `hobohm` program takes as input (1) a text-file containing a list of pairwise similarities and (2) a cutoff for deciding when two sequences are too similar.

The output (written to stdout) is a list of names that should be kept in the homology-reduced set. The algorithm aims to pick the maximally sized such set, given the cutoff.

It is also possible to use a list of pairwise *distances* instead of similarities. The cutoff is then interpreted as the minimum distance required in the output data set.


## Availability

The `hobohm` source code is available on GitHub: https://github.com/agormp/hobohm. The executable can be installed from PyPI: https://pypi.org/project/hobohm/

## Installation

```
python3 -m pip install hobohm
```

Upgrading to latest version:

```
python3 -m pip install --upgrade hobohm
```

## Dependencies

There are no dependencies (apart from the python standard library).

## Overview

### Input:

#### Option -s: pairwise similarities

(1) A text file containing pairwise *similarities*, one pair per line. All pairs of names must be listed.

```
name1 name2 similarity
name1 name3 similarity
...
```

(2) A cutoff value. Pairs of items that are more similar than this cutoff are taken to be redundant, and one of them will be removed in the final output.

#### Option -d: pairwise distances

(1) A text file containing pairwise *distances*, one pair per line. All pairs of names must be listed.

```
name1 name2 distance
name1 name3 distance
...
```

(2) A cutoff value. Pairs of items that are less distant than this cutoff are taken to be redundant, and one of them will be removed in the final output.

### Output:

A list of names of items that should be kept in the non-redundant set, written to stdout. This set contains no pairs of items that are more similar (less distant) than the cutoff. The algorithm aims at making the set the maximal possible size. This can occassionally fail if there are multiple items with the same number of "neighbors".

## Usage

```
Usage: hobohm [-s|-d] FILE -c CUTOFF [-k KEEPFILE]

Options:
  --version    show program's version number and exit
  -h, --help   show this help message and exit
  -s SIMFILE   file with pairwise similarities: name1 name2 sim
  -d DISTFILE  file with pairwise distances: name1 name2 dist
  -c CUTOFF    cutoff for deciding which pairs are neighbors
  -k KEEPFILE  file with names that must be kept (one name per line)
```