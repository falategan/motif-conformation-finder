

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
      the chain entry identifiers for each protein file, e.g. <br>
      <br>
     ``` 
     Protein ID, File Path, Chain
     1i6w, 1i6w.pdb.ent.gz, A
     1i6w, 1i6w.pdb.ent.gz, B
     1gci, 1gci.pdb.gz, A
     ...        
     ```
     
 4. Run motif-conformations.py from the command-line:
 
    ```
    ./motif-conformations.py [-h] [-q QUERY_LIST] [-g] motif structure_directory {pdb,cif} output_file
    ```
    - run `./motif-conformations.py -h` or `./motif-conformations.py --help` for command line <br>
      usage help.
    - use `-q` to specify the path to the chain query file you generated above in 3. If no file <br> 
      is specified, then the program will attempt to parse all files in the `structure_directory` <br>
      and seek `motif` in all the protein chains it finds
    - use the `-g` flag if the protein structure files are gzipped 
    - use `-o` to indicate the output file format 
    - `motif` is the amino acid sequence to seek. Note that the program currently only seeks exact <br>
      matches, and will not interpret promotif or regex motifs
    - `structure_directory` is the path of the the root directory that contains all the protein <br
       structure files
    - choose `cif` if the protein structure files are in the PDBx/MMCIF (.cif) format
    - choose `pdb` if the protein structure files are in the .pdb format
                                                                                                    
## Dive deeper :diving_mask::
   
Read the [API Documentation](https://motif-conformation-finder.readthedocs.io/en/latest/) to take a deeper look at the code.
