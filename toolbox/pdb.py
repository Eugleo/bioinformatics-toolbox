from typing import Counter
import Bio.PDB as P
from Bio.PDB.ResidueDepth import get_surface, min_dist
from scipy.spatial.distance import pdist
from numpy import isclose


def is_polar(aa):
    return aa in [
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLU",
        "GLN",
        "HIS",
        "LYS",
        "SER",
        "THR",
        "TYR",
    ]


class Structure:
    def __init__(self, id: str, path: str):
        p = P.PDBParser()
        self.structure = p.get_structure(id, path)
        self._width = None

    @property
    def summary(self):
        chains = 0
        residues = 0
        atoms = 0
        for model in self.structure:
            chains += len(model)
            for chain in model:
                residues += len(chain)
                for residue in chain:
                    atoms += len(residue)
        return {
            "models:": len(self.structure),
            "chains": chains,
            "residues": residues,
            "atoms": atoms,
        }

    @property
    def width(self):
        if self._width is None:
            self._width = max(
                pdist([a.get_coord() for a in self.structure.get_atoms()])
            )
        return self._width

    def search_around(self, ligand, radius: float, level: str):
        ns = P.NeighborSearch(P.Selection.unfold_entities(self.structure, "A"))
        return ns.search(ligand.get_coord(), radius, level)

    def get_residue_exposure(self, only_if=lambda x: True, msms_exec="msms"):
        exposed = Counter()
        buried = Counter()
        surface = get_surface(self.structure, MSMS=msms_exec)
        for r in self.structure.get_residues():
            aa = r.get_resname()
            if only_if(aa):
                mindist = min(min_dist(a.get_coord(), surface) for a in r.get_atoms())
                if mindist < 1.4:  # Anybody know why 1.399 is the minimum mindist?
                    exposed[aa] += 1
                else:
                    buried[aa] += 1
        total = sum(buried.values()) + sum(exposed.values())
        return {
            "buried": {
                "ratio": sum(buried.values()) / total,
                "distribution": list(buried.items()),
            },
            "exposed": {
                "ratio": sum(exposed.values()) / total,
                "distribution": list(exposed.items()),
            },
        }
