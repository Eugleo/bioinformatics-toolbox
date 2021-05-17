import dataclasses
from typing import Iterator

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


@dataclasses.dataclass(init=False)
class Sequence:
    description: str
    name: str
    sequence: str

    def __init__(self, bioseq: SeqRecord):
        self.description = bioseq.description
        self.name = bioseq.name
        self.sequence = str(bioseq.seq)

    def __len__(self):
        return len(self.sequence)

    def __iter__(self):
        return self.sequence.__iter__()

    def __getitem__(self, item):
        return self.sequence.__getitem__(item)


def parse(path) -> Iterator[Sequence]:
    return (Sequence(s) for s in SeqIO.parse(path, "fasta"))
