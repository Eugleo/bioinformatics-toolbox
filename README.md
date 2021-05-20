# Evžen's bioinformatic bag of tricks

You can see the API and most of the functionality in the [samples notebook](https://github.com/Eugleo/bioinformatics-toolbox/blob/main/samples.ipynb).

If you want to play around with the code, you'll need to have [poetry](https://python-poetry.org) and [msms](http://mgltools.scripps.edu/packages/MSMS) installed and available on your PATH. You'll also need Python 3.9.0 or later. Then you can simply run:

```
> git clone https://github.com/Eugleo/bioinformatics-toolbox
> cd bioinformatics-toolbox
> poetry install
> poetry run jupyter notebook samples.ipynb
```

Which will open the mentioned notebook in edit mode.

## Functionality

The first 7 assignments are implemented. Most of the implementations rely heavily on the `biopython` library. The functionality is split into five modules:

1. `toolbox.fasta` with the fasta parser
2. `toolbox.hamming` with the Hamming distance calulation
3. `toolbox.alignment` with sequence alignment implemented by edit distance
4. `toolbox.msa` with Clustal-files parser and utilities
5. `toolbox.pdb` with PDB-files parser and utilities

Each feature is showcased in the `samples.ipynb` notebook, which also serves as a kind of documentation.