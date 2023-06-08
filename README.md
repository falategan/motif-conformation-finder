

<div align = center>
  
# motif conformation finder
find all instances of an amino acid sequence motif in protein structures and parse their internal coordinates
  
</div>
 
<br>
 
<div align = right>
 
  
[![Python](https://img.shields.io/badge/Python-3.11-informational?style=flat-square&logo=appveyor)](https://www.python.org/)
[![License](https://img.shields.io/github/license/falategan/motif-conformation-finder?style=flat-square)](/LICENCE)
<br>
<br>
[API Documentation](https://motif-conformation-finder.readthedocs.io/en/latest/)
</div>

## Dependencies
- ![Python 3.11](https://www.python.org/)
- ![Biopython](biopython.org/)

## Usage:

 1. Install dependencies

 2. Download the protein structures from the PDB

 3. Write chain query file (optional)
    - A list specifying which chains to search, and which protein structure files to search in.<br>
      This file is a long-format csv file containing the protein entry ids, the path to each <br>
      protein file from the root directory that contains all the protein structure files, and <br>
      the chain entry identifiers for each protein file, e.g.
      
            ``` 
            Protein ID, File Path, Chain
            1i6w, 1i6w.pdb.ent.gz, A
            1i6w, 1i6w.pdb.ent.gz, B
            1gci, 1gci.pdb.gz, A
             ```
          

 4. Run motif-conformations.py from the command-line:
 
    ```
    python [-h] [-q QUERY_LIST] [-g] motif structure_directory {pdb,cif} output_file
    ```
