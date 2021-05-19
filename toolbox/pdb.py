import Bio.PDB as P
from scipy.spatial.distance import pdist


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
