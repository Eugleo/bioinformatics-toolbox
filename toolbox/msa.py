from Bio import AlignIO
import itertools
from collections import Counter


class Clustal:
    def __init__(self, path):
        self.alignments = AlignIO.read(path, "clustal")
        self._cols = self.alignments.get_alignment_length()
        self._col_scores = [None] * self._cols
        self._col_conservations = [None] * self._cols

    # Works with numpy-style indexing to retrieve whole columns
    def __getitem__(self, index):
        if isinstance(index, str):
            seq = next((a for a in self.alignments if a.id == index), None)
            if seq is None:
                raise KeyError(f"Sequence with id {index} doesn't exist")
            return seq
        else:
            return self.alignments[index]

    def score(self, matrix, indel, cols=None):
        total = 0
        for i in cols if cols is not None else range(self._cols):
            if self._col_scores[i] is None:
                self._col_scores[i] = Clustal._score_one_column(
                    self[:, i], matrix, indel
                )
            total += self._col_scores[i]
        return total

    def _score_one_column(seqs, matrix, indel):
        def score_pair(s1, s2):
            if s1 == "-" and s2 == "-":
                return 0
            elif s1 == "-" or s2 == "-":
                return indel
            else:
                return matrix[(s1, s2)]

        return sum(score_pair(*p) for p in itertools.combinations(seqs, 2))

    def score_conservation(self, position):
        if self._col_conservations[position] is None:
            count = max(v for k, v in Counter(self[:, position]).items() if k != "-")
            self._col_conservations[position] = count / len(self.alignments)
        return self._col_conservations[position]

    def k_most_conserved(self, k):
        return sorted(
            ((i, self.score_conservation(i)) for i in range(self._cols)),
            key=lambda t: t[1],
            reverse=True,
        )[:k]
