from toolbox.fasta import Sequence


def distance(s1: Sequence, s2: Sequence):
    if len(s1) != len(s2):
        raise ValueError(
            f"The sequences should have identical lengths, but the lengths are {len(s1)}, {len(s2)}"
        )
    return sum(1 for a, b in zip(s1, s2) if a != b)
