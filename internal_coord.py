"""Read and write CSV files and parse protein structures

Author: Adriaan Lategan
"""

from dataclasses import dataclass
from typing import Callable, Generator, Iterable, Type

from Bio import PDB
from Bio.PDB import Chain, Residue
from Bio.PDB.ic_data import ic_data_backbone, ic_data_sidechains, \
    ic_data_sidechain_extras
from Bio.PDB.internal_coords import IC_Residue, AtomKey

DIHEDRAL_ANGLE_KEYS = ['phi',
                       'psi',
                       'omega',
                       'chi1',
                       'chi2',
                       'chi3',
                       'chi4',
                       'chi5'
                       ]


@dataclass
class CoordinateData:
    """ Object for protein structure coordinates
    Attributes
    ----------
    coordinate_type : str
        type of coordinate, for example cartesian, dihedral angle, bond angle
        or bond length
    coordinate_id : str
        identifier for a specific coordinate value, such as x, y, z, CA:CB,
        N:C, psi, phi
    coordinate_value : float
        floating point value of the coordinate
    """
    coordinate_type: str
    coordinate_id: str
    coordinate_value: float


@dataclass
class ResidueData:
    """ Object identifying an amino acid residue and recording its protein
    structure coordinates
    
    Attributes
    ----------
    protein : str
        pdb entity ID of the protein structure
    model : int
        pdb model number
    chain : str
        pdb polypeptide instance ID
    position : int
        residue position index in the protein chain
    residue_name : str
        amino acid type as its 3-letter name, e.g. ALA or PRO
    coordinates : Generator[CoordinateData, None, None]
        protein structure coordinates describing the residue as atomic
        positions or internal angles
    """
    protein: str
    model: int
    chain: str
    position: int
    residue_name: str
    coordinates: Generator[CoordinateData, None, None]


def get_triads(ic_data: Iterable) -> list[tuple[str, str, str]]:
    """Extract groups of 3 atoms from the biopython internal coordinate data

    Parameters
    ----------
    ic_data : Iterable
        Iterable containing groups of bonded atoms that describe internal
        coordinates


    Returns
    -------
    list[tuple[str, str, str]]
        list of groups of 3 bonded atoms
    """
    return [hedron
            for hedron in ic_data
            if 'H' not in str(hedron) and len(hedron) == 3
            ]


def get_bonds(hedra: list[tuple]) -> list[str]:
    """Extract groups of 2 atoms from the biopython internal coordinate data

        Parameters
        ----------
        hedra : list[tuple]
            List of groups of atoms that describe internal
            coordinates


        Returns
        -------
        list[str]
            list of bonded atom pairs
        """
    return list({bond
                 for hedron in hedra
                 for bond in [f'{hedron[0]}:{hedron[1]}',
                              f'{hedron[1]}:{hedron[2]}'
                              ]
                 }
                )


def get_sidechain_triads(sidechain_atom_groups: dict[str, list[tuple]]
                         ) -> dict[str, list[tuple]]:
    """Extract groups of 3 atoms from the biopython internal coordinate
    sidechain data

    Parameters
    ----------
    sidechain_atom_groups : dict[str, list[tuple]]
       list of bonded atom groups for the sidechain of each standard amino acid

    Returns
    -------
    dict[str, list[tuple]]
        list of groups of 3 bonded atoms for the sidechain of each standard
        amino acid
    """
    for residue, groups in ic_data_sidechain_extras.items():
        for group in groups:
            sidechain_atom_groups[residue] += tuple(group)
    return {residue: get_triads(groups)
            for residue, groups in sidechain_atom_groups.items()
            }


FUNCTION_TYPE: Type = Callable[
    [tuple[AtomKey, AtomKey]
     | tuple[AtomKey, AtomKey, AtomKey]
     | tuple[AtomKey, AtomKey, AtomKey, AtomKey]
     | str
     ], str | float | None
]


def get_coordinate(residue_letter: str,
                   func: FUNCTION_TYPE,
                   backbone_hedra: list[str],
                   sidechain_hedra: dict[str, str]
                   ) -> dict[str, float]:
    """internal coordinate values for the specified atom chains

    Parameters
    ----------
    residue_letter : str
        single letter amino acid name
    func : Callable
        a method of Bio.PDB.internal_coords.IC_Residue to retrieve the desired
        internal coordinate
    backbone_hedra : list[str]
        list of atom chains in the residue's backbone
    sidechain_hedra
        list of atom chains in the residue's sidechain
    Returns
    -------
    dict[str, float]
        dictionary with the coordinate value for each atom chain string

    """
    moieties = backbone_hedra
    side_chain_moieties = sidechain_hedra.get(residue_letter, [])
    if isinstance(side_chain_moieties, str):
        moieties += [side_chain_moieties, ]
    else:
        moieties += side_chain_moieties
    coordinates = {}
    for moiety_key in moieties:
        value = func(moiety_key)
        if value:
            coordinates[moiety_key] = value
    return coordinates


def split_coordinates(coordinate_type: str,
                      coordinates: dict[str, float]
                      ) -> Generator[CoordinateData, None, None]:
    """generate CoordinateData object of each coordinate value

    Parameters
    ----------
    coordinate_type : str
        type of coordinate, for example cartesian, dihedral angle, bond angle
        or bond length
    coordinates : dict[str, float]
        dictionary of coordinate identifiers and values
    Yields
    ------
    CoordinateData
        Object for protein structure coordinates
    """
    for coordinate_id, coordinate_value in coordinates.items():
        yield CoordinateData(coordinate_type, coordinate_id, coordinate_value)


class InternalCoordinates:
    def __init__(
        self,
        backbone_atom_chains: tuple[tuple] = ic_data_backbone,
        sidechain_atom_chains: dict[str, tuple] = ic_data_sidechains.copy()
                ) -> None:
        """Object for accessing residue internal coordinates
        
        Parameters
        ----------
        backbone_atom_chains: tuple[tuple]
            backbone atom chains in a residue
        sidechain_atom_chains: dict[str, tuple]
            sidechain atom chains in each standard amino acid
        
        Attributes
        ----------
        backbone_bonds : list[str]
            backbone atom pair identifiers
        backbone_angle_keys : list[str]
            identifiers for chains of 3 backbone atoms
        sidechain_bonds : dict[str, list[str]]
            sidechain atom pair identifiers for each standard aminoacid
        sidechain_angle_keys : dict[str, list[str]]
            identifiers for chains of 3 sidechain atoms for each standard amino
            acid
        """
        
        backbone_atom_triads = get_triads(backbone_atom_chains)
        self.backbone_bonds = get_bonds(backbone_atom_triads)
        self.backbone_angle_keys = [':'.join(triad) for triad in
                                    backbone_atom_triads]
        sidechain_atom_triads = get_sidechain_triads(sidechain_atom_chains)
        self.sidechain_bonds = {residue: get_bonds(hedra)
                                for residue, hedra in
                                sidechain_atom_triads.items()
                                }
        self.sidechain_angle_keys = {residue: ':'.join(triad)
                                     for residue, triads in
                                     sidechain_atom_triads.items()
                                     for triad in triads
                                     }
    
    def get_coordinates(self,
                        chain: Chain,
                        residue: Residue
                        ) -> ResidueData:
        """ internal coordinates of the residue

        Parameters
        ----------
        chain : Chain
            Bio.PDB chain object
        residue : Residue
            Bio.PDB residue object
        Returns
        -------
        ResidueData
            Object identifying an amino acid residue and recording its protein
            structure coordinates
        """
        chain.atom_to_internal_coordinates(verbose=False)
        ic_residue: IC_Residue = chain[residue.id].internal_coord
        residue_name = residue.resname
        residue_letter = PDB.Polypeptide.protein_letters_3to1[residue_name]
        
        bond_angles = get_coordinate(residue_letter,
                                     ic_residue.get_angle,
                                     self.backbone_angle_keys,
                                     self.sidechain_angle_keys
                                     )
        dihedral_angles = get_coordinate(residue_letter,
                                         ic_residue.get_angle,
                                         DIHEDRAL_ANGLE_KEYS,
                                         {}
                                         )
        bond_lengths = get_coordinate(residue_letter,
                                      ic_residue.get_length,
                                      self.backbone_bonds,
                                      self.sidechain_bonds)
        
        coordinate_dictionary = {'bond_angles': bond_angles,
                                 'dihedral_angles': dihedral_angles,
                                 'bond_lengths': bond_lengths
                                 }
        coordinates = (coordinate_data
                       for name, coordinate in coordinate_dictionary.items()
                       for coordinate_data in
                       split_coordinates(name, coordinate)
                       )
        
        residue_data = ResidueData(residue.full_id[0],
                                   residue.full_id[1],
                                   residue.full_id[2],
                                   residue.full_id[3][1],
                                   residue_name,
                                   coordinates
                                   )
        return residue_data
