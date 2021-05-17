import dataclasses

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse(path):
    return (Sequence(s) for s in SeqIO.parse(path, "fasta"))


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
