"""Calculate the internal coordinates of each instance of a sequence motif in
protein structures

Author: Adriaan Lategan
"""

import argparse

from Bio.PDB.Polypeptide import PPBuilder

from internal_coord import InternalCoordinates
from pdb_io import PdbQueries, PdbReader, PdbQueryCsv, CsvWriter

HAS_HEADER = True


def get_arguments() -> argparse.Namespace:
    """Fetch command-line arguments

    Returns
    -------
    argparse.Namespace
    """
    argument_parser = argparse.ArgumentParser("Calculate the internal "
                                              "coordinates of each instance "
                                              "of a sequence motif in protein "
                                              "structures")
    argument_parser.add_argument('motif',
                                 type=str,
                                 help='the amino acid sequence to find and '
                                      'extract from the protein structures'
                                 )
    argument_parser.add_argument('-q', '--query_list',
                                 type=str,
                                 help='path to the file specifying the names, '
                                      'file paths, and chain ids of the '
                                      'protein structures to search for the '
                                      'motif. If no query_list is specified,'
                                      'it will read all files in '
                                      'structure_directory.'
                                 )
    argument_parser.add_argument('structure_directory',
                                 type=str,
                                 help='path to the directory containing '
                                      'protein structure files'
                                 )
    argument_parser.add_argument('structure_format',
                                 choices=['pdb', 'cif'],
                                 type=str,
                                 default='pdb',
                                 help='format of the protein structure files: '
                                      'either "pdb" or "cif"'
                                 )
    argument_parser.add_argument('-g', '--gzipped',
                                 action='store_true',
                                 help='structure files are compressed using '
                                      'gzip'
                                 )
    argument_parser.add_argument('output_file',
                                 type=str,
                                 help='file path for the output csv file'
                                 )
    arguments = argument_parser.parse_args()
    return arguments


def motif_to_coordinates(
    motif: str,
    query_list: str,
    pdb_directory: str,
    pdb_format: str,
    gzipped: bool,
    output_file: str
    ) -> None:
    """ Read the protein structure files in the specified directory,
    search for the motif sequence in each structure and output the internal
    coordinates of the motifs as a csv file
    
    Parameters
    ----------
    motif : str
        an amino acid sequence to find in the protein structures
    query_list : str
        a csv file specifying protein file paths and polymer instances
    pdb_directory : str
        a directory containing protein structure files
    pdb_format : str
        the format of the protein structure files. Either "pdb" or "cif"
    gzipped : bool
        true if the protein structure files are compressed with gzip, false if
        uncompressed
    output_file : str
        name of the file to which the internal coordinates are written
    """
    
    has_header = HAS_HEADER
    pdb_reader = PdbReader(pdb_directory, pdb_format, gzipped)
    internal = InternalCoordinates()
    
    if query_list:
        chain_file = PdbQueryCsv(query_list, has_header)
        pdb_queries: PdbQueries = chain_file.read
    else:
        pdb_queries: PdbQueries = pdb_reader.directory_queries
    
    residue_data = (internal.get_coordinates(polypeptide.chain, residue)
                    for query in pdb_queries
                    for polypeptide in
                    query.get_polypeptides(pdb_reader, PPBuilder())
                    for fragment in polypeptide.find_motif(motif)
                    for residue in fragment
                    )
    fields = ['Protein',
              'Model',
              'Chain',
              'Position',
              'Residue Name',
              'Coordinate Type',
              'Coordinate ID',
              'Coordinate Value'
              ]
    csv = CsvWriter(output_file, fields)
    csv.write_headings()
    for residue in residue_data:
        for coordinate in residue.coordinates:
            field_values = [residue.protein,
                            f'{residue.model}',
                            residue.chain,
                            f'{residue.position}',
                            residue.residue_name,
                            coordinate.coordinate_type,
                            coordinate.coordinate_id,
                            f'{coordinate.coordinate_value:f}'
                            ]
            csv.write_line(field_values)
    csv.close()


if __name__ == "__main__":
    arguments = get_arguments()
    motif_to_coordinates(arguments.motif,
                         arguments.query_list,
                         arguments.structure_directory,
                         arguments.structure_format,
                         arguments.gzipped,
                         arguments.output_file
                         )
