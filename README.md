# hobohm: command line program for selecting representative subset of data, based on list of pairwise similarities (or distances) between items.

[![PyPI downloads](https://static.pepy.tech/personalized-badge/hobohm?period=total&units=international_system&left_color=grey&right_color=blue&left_text=downloads)](https://pepy.tech/project/hobohm)
![](https://img.shields.io/badge/version-1.0.3-blue)

The `hobohm` program aims to select a representative subset from a collection of items for which the pairwise similarities are known.

The program takes as input (1) a text-file containing a list of pairwise similarities between items in a data set, and (2) a cutoff for deciding when two items are too similar (i.e., when they are "neighbors").

The output (written to stdout) is a list of names that should be kept in the subset. No retained items are neighbors, and the algorithm aims to pick the maximally sized such set, given the cutoff.

It is also possible to use a list of pairwise *distances* instead of similarities. The cutoff is then interpreted as the *minimum distance* required in the selected subset.

The "Hobohm" algorithm was originally created with the purpose of selecting homology-reduced sets of protein data from larger datasets. "Homology-reduced" here means that the resulting data set should contain no pairs of sequences with high sequence identity:

["Selection of representative protein data sets", Protein Sci. 1992. 1(3):409-17](https://pubmed.ncbi.nlm.nih.gov/1304348/).

This command-line program implements algorithm 2 from that paper, and can be applied to any type of data for which pairwise similarities (or distances) can be defined.



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

(1) A text file containing pairwise *similarities*, one pair per line. All pairs of names must be listed. The similarity matrix is assumed to be symmetric, and it is only necessary to list one direction for each pair of names.

```
name1 name2 similarity
name1 name3 similarity
...
```

(2) A cutoff value. Pairs of items that are *more similar* than this cutoff are taken to be redundant, and at least one of them will be removed in the final output.

#### Option -d: pairwise distances

(1) A text file containing pairwise *distances*, one pair per line. All pairs of names must be listed. The distance matrix is assumed to be symmetric, and it is only necessary to list one direction for each pair of names.

```
name1 name2 distance
name1 name3 distance
...
```

(2) A cutoff value. Pairs of items that are *less distant* than this cutoff are taken to be redundant, and at least one of them will be removed in the final output.

### Output:

A list of names of items that should be kept in the non-redundant set, written to stdout. This set contains no pairs of items that are more similar (less distant) than the cutoff. The algorithm aims at making the set the maximal possible size. This can occassionally fail if there are multiple items with the same number of "neighbors" and the order of removal of items has an impact.

## Usage

```
usage: hobohm.py [-h] [-s | -d] [-c CUTOFF] [-k KEEPFILE] PAIRFILE

Selects representative subset of data based on list of pairwise similarities (or
distances), such that no retained items are close neighbors

positional arguments:
  PAIRFILE     file containing the similarity (option -s) or distance (option -d) for each
               pair of items: name1 name2 value

optional arguments:
  -h, --help   show this help message and exit
  -s           values in PAIRFILE are similarities (larger values = more similar)
  -d           values in PAIRFILE are distances (smaller values = more similar)
  -c CUTOFF    cutoff for deciding which pairs are neighbors
  -k KEEPFILE  file with names of items that must be kept (one name per line)
  ```

## Usage examples

### Select items such that max pairwise similarity is 0.65

```
hobohm -s -c 0.65 pairsims.txt > nonredundant.txt
```

### Select items such that minimum pairwise distance is 10

```
hobohm -d -c 10 pairdist.txt > nonredundant.txt
```

### Select items such that max pairwise similarity is 0.3, while keeping items in keeplist.txt

```
hobohm -s -c 0.3 -k keeplist.txt pairsims.txt > nonredundant.txt
```
