"""Read and write CSV files and parse protein structures

Author: Adriaan Lategan
"""

import gzip
import pathlib
import warnings
from dataclasses import dataclass, field
from typing import Iterator, Generator, Iterable

from Bio import PDB
from Bio.PDB import Structure, Chain, Residue


def custom_format(message: Warning | str,
                  category: type[Warning],
                  filename: str,
                  lineno: int,
                  *args, **kwargs
                  ) -> str:
    """ Overwrite the format of Warnings
    
    Returns
    -------
    str
    """
    return f'{filename}:{lineno}: {category.__name__}: {message}\n'


warnings.formatwarning = custom_format


@dataclass
class PolypeptideEntry:
    """ Object representing a contiguous polymer of amino acids, with methods
    to extract amino acid sequences from it
    
    Attributes
    ----------
    pdb_id : str
        pdb entity ID of the protein structure
    model : int
        pdb model number
    chain : PDB.Chain.Chain
         protein chain entity object
    polypeptide : PDB.Polypeptide
         a contiguous chain of amino acids from the protein chain entry
    """
    pdb_id: str
    model: int
    chain: PDB.Chain.Chain
    polypeptide: PDB.Polypeptide
    
    @property
    def sequence(self) -> str:
        """
        
        Returns
        -------
        str
            amino acid sequence of the polypeptide instance
        """
        return self.polypeptide.get_sequence()
    
    def find_motif(self,
                   motif: str
                   ) -> Generator[list[PDB.Residue.Residue], None, None]:
        """find each instance of the sequence motif in the polypeptides

        Parameters
        ----------
        motif : str

        Yields
        -------
        list[PDB.Residue.Residue]
            residues matching the motif from each polypeptide

        """
        end = 0
        start = self.sequence.find(motif, end)
        if start == -1:
            return
        while start != -1:
            end = start + len(motif)
            yield self.polypeptide[start:end]
            start = self.sequence.find(motif, end)


class PdbReader:
    def __init__(self,
                 pdb_directory: str,
                 pdb_type: str = 'cif',
                 gzipped: bool = True
                 ):
        """Object for parsing protein structures in pdb or cif format from a
        specific directory
        
        Parameters
        ----------
        pdb_directory : str
             a directory containing protein structure files
        pdb_type : str
            the format of the protein structure files. Default "cif"
        gzipped : bool
            true if the protein structure files are compressed with gzip, false
            if uncompressed
        
        Attributes
        ----------
        directory_queries
        pdb_directory : pathlib.Path
            path to the directory containing protein structure files
        text_handler : Callable
            method for opening text stream
        parser : PDB.MMCIFParser or PDB.PDBParser
        
        """
        self.pdb_directory = pathlib.Path(pdb_directory)
        if not self.pdb_directory.exists():
            error = f"The directory {self.pdb_directory} does not exist."
            raise FileNotFoundError(error)
        self.text_handler = gzip.open if gzipped else open
        if pdb_type == 'cif':
            self.parser = PDB.MMCIFParser(QUIET=True, auth_chains=False)
        elif pdb_type == 'pdb':
            self.parser = PDB.PDBParser(QUIET=True)
        else:
            error = f'Invalid structure format {pdb_type}. The structure ' \
                    f'format should either be "cif" or "pdb"'
            raise ValueError(error)
    
    def read_file(self,
                  protein_id: str,
                  file_name: str
                  ) -> PDB.Structure.Structure | None:
        """ Parse a protein structure file and return a biopython PDB structure

        Parameters
        ----------
        protein_id : str
            pdb entity ID of the protein structure
        file_name : str
            name of the protein structure file

        Returns
        -------
        PDB.Structure.Structure or None
            Biopython structure entity
        """
        file_path = self.pdb_directory / file_name
        structure = None
        if not file_path.exists():
            warnings.warn(f"The pdb file {file_path} does not exist.")
            
            return structure
        
        with self.text_handler(file_path, "rt") as file:
            try:
                structure = self.parser.get_structure(protein_id, file)
            except (EOFError, TypeError, ValueError) as error:
                warnings.warn(f'Encountered error "{error}" while parsing '
                              f'file "{file_path}"')
        
        return structure
    
    @property
    def directory_queries(self):
        """read each file in the directory and get the protein chain IDs and
        each protein structure file

        Returns
        -------
        PdbQueries
            the path to each protein structure file and the IDs of the
            polypeptide chains in that protein
        """
        queries = PdbQueries()
        for file in self.pdb_directory.iterdir():
            try:
                if not file.is_file():
                    continue
                path = file.name
                protein_id = path[:path.index(".")]
                structure = self.read_file(protein_id, path)
                if not structure:
                    continue
                for chain in structure.get_chains():
                    queries.add_query(protein_id, path, chain.id)
            except OSError as error:
                warnings.warn(f"Encountered OSError {error} for file {file}")
        return queries


@dataclass(slots=True)
class PdbFileQuery:
    """ Object representing a specific proteins structure file and the amino
    protein chain entities that should be parsed from it
    
    Attributes
    ----------
    protein_id: str
        pdb entity ID of the protein structure
    file_path: str
        path to the protein structure file
    chain_ids: list[str]
        list of protein chain entities to parse
    """
    protein_id: str
    file_path: str
    chain_ids: list[str]
    
    def __iter__(self):
        return self.chain_ids.__iter__()
    
    def add_chain(self, chain_id: str | list[str]):
        self.chain_ids += chain_id
    
    def get_structure(self,
                      pdb_reader: PdbReader
                      ) -> PDB.Structure.Structure | None:
        """ Parse a protein structure file and return a biopython PDB structure

        Parameters
        ----------
        pdb_reader : str
            protein structure file parser

        Returns
        -------
        PDB.Structure.Structure or None
            biopython structure entity
        """
        return pdb_reader.read_file(self.protein_id, self.file_path)
    
    def get_polypeptides(self,
                         pdb_reader: PdbReader,
                         builder: PDB.Polypeptide.PPBuilder |
                                  PDB.Polypeptide.CaPPBuilder
                         ) -> Generator[PolypeptideEntry, None, None]:
        """ Identify contiguous chains of amino acids in the protein chain
        entity

        Parameters
        ----------
        pdb_reader : PdbReader
            protein structure file parser
            
        builder :  PDB.Polypeptide.PPBuilder or PDB.Polypeptide.CaPPBuilder
            polypeptide constructor

        Yields
        -------
        PolypeptideEntry
            contiguous chains of amino acids from the protein chain entry
        """
        structure = self.get_structure(pdb_reader)
        if not structure:
            return
        
        for model in structure:
            for chain_id in self.chain_ids:
                if chain_id not in model:
                    warnings.warn(f'Model {model.id} of structure '
                                  f'{structure.id} '
                                  f'does not contain Chain {chain_id}.'
                                  )
                    continue
                chain = model[chain_id]
                for polypeptide in builder.build_peptides(chain):
                    entry = PolypeptideEntry(structure.id,
                                             model.id,
                                             chain,
                                             polypeptide
                                             )
                    yield entry


@dataclass(slots=True)
class PdbQueries:
    """
    List with unique protein structure file query entries
    
    Attributes
    ----------
    query_list : list[PdbFileQuery]
    """
    query_list: list = field(default_factory=list)
    
    def add_query(self, protein_id: str, path_string: str, chain: str):
        """ Creates a new PdbFileQuery object, or appends a protein chain
        entity ID to an existing query

        Parameters
        ----------
        protein_id : str
            pdb entity ID of the protein structure
        path_string : str
            name of the protein structure file
        chain : str
            pdb polypeptide instance ID
        Returns
        -------
        None
        """
        if protein_id not in self:
            query = PdbFileQuery(protein_id, path_string, [chain])
            self.query_list.append(query)
            return
        if chain not in self[protein_id]:
            self[protein_id].add_chain(chain)
        return
    
    def __contains__(self, item):
        return any(query.protein_id == item for query in self)
    
    def __getitem__(self, item: str) -> PdbFileQuery:
        for query in self:
            if query.protein_id == item:
                return query
        raise KeyError
    
    def __iter__(self) -> Iterator[PdbFileQuery]:
        return self.query_list.__iter__()


class PdbQueryCsv:
    def __init__(self, chain_list: str, has_header: bool = True) -> None:
        """Object for parsing csv files listing protein entity IDs, paths to
        protein structure files, and protein chain IDs
        
        Parameters
        ----------
        chain_list : str
            path to a csv file listing protein structures to parse
        has_header: bool
            flag to indicate whether the csv file has column headings.
            Default: True
            
        Attributes
        ----------
        read
        path : pathlib.Path
            path to a csv file listing protein structures to parse
        has_header
            flag to indicate whether the csv file has column headings.
        
        """
        self.path = pathlib.Path(chain_list)
        if not self.path.exists():
            raise FileNotFoundError(f"Chain list file at {self.path} not "
                                    f"found.")
        self.has_header = has_header
    
    @property
    def read(self) -> PdbQueries:
        """Read a csv file listing PDB_IDs, paths, and chains

        Returns
        -------
        PdbQueries
        """
        discard = self.has_header
        queries = PdbQueries()
        with open(self.path, 'r') as file:
            for line in file:
                if not discard:
                    protein_id, path, chain = line.rstrip("\n").split(',')
                    queries.add_query(protein_id, path, chain)
                discard = False
        return queries


class CsvWriter:
    
    def __init__(self, path_string: str, fields: Iterable[str]):
        """Object for writing values to defined fields in a csv file
        
        Parameters
        ----------
        path_string : str
            path of the output file to write
        fields : Iterable[str]
            column names of csv file
            
        Attributes
        ---------
        output_handle : IO
            text stream for writing output file
        fields : Iterable[str]
            column names of csv file
        """
        output_path = pathlib.Path(path_string)
        self.output_handle = open(output_path, 'w')
        self.fields = fields
    
    def write_headings(self):
        header = ','.join(self.fields) + '\n'
        self.output_handle.write(header)
    
    def write_line(self, field_values: Iterable[str]):
        """add a row to the csv format table and write it to the output file
        
        Parameters
        ----------
        field_values : Iterable[str]
            list containing a string value for each column
        """
        self.output_handle.write(f'{",".join(field_values)}\n')
    
    def close(self):
        """
        Close the TextIO stream
        """
        self.output_handle.close()
