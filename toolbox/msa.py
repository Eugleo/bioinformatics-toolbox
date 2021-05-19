from Bio import AlignIO
import itertools


class Clustal:
    def __init__(self, path):
        self.alignments = AlignIO.read(path, "clustal")

    # Works with numpy-style indexing to retrieve whole columns
    def __getitem__(self, index):
        if isinstance(index, str):
            seq = next((a for a in self.alignments if a.id == index), None)
            if seq is None:
                raise KeyError(f"Sequence with id {index} doesn't exist")
            return seq
        else:
            return self.alignments[index]

    def score(self, matrix, cols=None):
        if cols is None:
            cols = range(self.alignments.get_alignment_length())
        return sum(Clustal._score_one_column(self[:, i], matrix) for i in cols)

    def _score_one_column(seqs, matrix):
        return sum(
            matrix[p]  # if p in matrix else matrix[tuple(reversed(p))]
            for p in itertools.combinations((s for s in seqs if s[0] != "-"), 2)
        )
