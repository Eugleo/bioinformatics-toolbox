# Evžen's bioinformatic bag of tricks

You'll need to have [poetry](https://python-poetry.org) and [msms](http://mgltools.scripps.edu/packages/MSMS) installed and available on your PATH. Then you can simply run:

```sh
> git clone https://github.com/Eugleo/bioinformatics-toolbox
> cd bioinformatics-toolbox
> poetry install
> poetry run jupyter notebook samples.ipynb
```

A jupyter notebook will open, where you can play around with the package.

## Functionality

The first 7 assignments are implemented. Most of the implementations rely heavily on the `biopython` library. The functionality is split into five modules:

1. `toolbox.fasta` with the fasta parser
2. `toolbox.hamming` with the Hamming distance calulation
3. `toolbox.alignment` with sequence alignment implemented by edit distance
4. `toolbox.msa` with Clustal-files parser and utilities
5. `toolbox.pdb` with PDB-files parser and utilities

Each feature is showcased in the `samples.ipynb` notebook, which also serves as a kind of documentation.