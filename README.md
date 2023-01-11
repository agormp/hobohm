# hobohm: command line program for selecting representative subset of data, based on list of pairwise similarities (or distances) between items.

[![PyPI downloads](https://static.pepy.tech/personalized-badge/hobohm?period=total&units=international_system&left_color=grey&right_color=blue&left_text=downloads)](https://pepy.tech/project/hobohm)
![](https://img.shields.io/badge/version-2.0.1-blue)

---------------------------------------------------------------- 

**Note:** The [greedysub program](https://github.com/agormp/greedysub) implements a better algorithm (typically giving larger subsets), and should be used instead. The `hobohm` program works, but is no longer maintained).

----------------------------------------------------------------

The `hobohm` program aims to select a non-redundant subset of DNA- or protein-sequences, such that all pairwise sequence identities
are below a given threshold.

The program takes as input (1) a text-file containing a list of pairwise similarities between sequences (`name1 name2 similarity`), and (2) a cutoff for deciding when two sequences are too similar (i.e., when they are "neighbors").

The output (written to file) is a list of names that should be kept in the subset. No retained items are neighbors, and the algorithm aims to pick the maximally sized such set, given the cutoff. (Note that this is a hard problem, and this heuristic is not optimal. See  notes on computational intractibility of the problem and performance of heuristics in the [greedysub README](https://github.com/agormp/greedysub)).

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

`hobohm` relies on the [pandas package](https://pandas.pydata.org), which is automatically included when using pip to install.


## Usage

```
usage: hobohm [-h] [--val VALUETYPE] [-c CUTOFF] [-k KEEPFILE] INFILE OUTFILE

Select non-redundant subset of DNA or protein-sequences, such that all pairwise
sequence identities are below threshold.

positional arguments:
  INFILE           input file containing similarity or distance for each pair of
                   items: name1 name2 value
  OUTFILE          output file contatining neighborless subset of items (one name per
                   line)

options:
  -h, --help       show this help message and exit
  --val VALUETYPE  specify whether values in INFILE are distances (--val dist) or
                   similarities (--val sim)
  -c CUTOFF        cutoff value for deciding which pairs are neighbors
  -k KEEPFILE      (optional) file with names of items that must be kept (one name per
                   line)
```
