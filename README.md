# BindSite: Protein Binding Site Analysis Tool

BindSite is a Python tool that analyzes protein structures to identify and characterize potential binding sites using fpocket integration.

## Features

- Automated binding site detection using fpocket
- Analysis of protein flexibility through:
  - B-factor analysis
  - Secondary structure identification
- Basic binding site characterization including volume and pocket scores
- Command-line interface for easy analysis

## Prerequisites

- Python 3.8+
- fpocket must be installed on your system (https://github.com/Discngine/fpocket)

## Installation

### Windows Users
BindSite requires fpocket, which is a Linux tool. There are two ways to run it on Windows:

1. Using WSL (Windows Subsystem for Linux) - Recommended:
```bash
# Install WSL if you haven't already
wsl --install

# After WSL installation and restart, open Ubuntu terminal and run:
sudo apt-get update
sudo apt-get install fpocket

# Clone and install BindSite
git clone https://github.com/adam-davidson128/bindsite.git
cd bindsite
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

2. Using Conda (Alternative method):
```bash
# Create new conda environment with fpocket
conda create -n bindsite
conda activate bindsite
conda install -c bioconda fpocket

# Install BindSite
git clone https://github.com/yourusername/bindsite.git
cd bindsite
pip install -r requirements.txt
```

### Linux
```bash
# Install fpocket
sudo apt-get install fpocket  # Ubuntu/Debian

# Then install BindSite
git clone https://github.com/yourusername/bindsite.git
cd bindsite

# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

BindSite can be run from the command line:

```bash
python bindsite.py --pdb protein.pdb --ligand ligand.mol2 --pocket 0
```

Arguments:
- `--pdb`: Path to protein structure file (default: protein.pdb)
- `--ligand`: Path to ligand file (default: ligand.mol2)
- `--pocket`: ID of pocket to analyze (default: auto-select best druggability score)

## Output

The program provides:
1. Pocket detection results from fpocket
2. Flexibility analysis including:
   - B-factors for each atom
   - Secondary structure assignments
3. Binding site analysis including:
   - Pocket volume
   - Pocket score
   - Binding site rank
   - Total number of pockets found

## Example

```bash
$ python bindsite.py --pdb 2c6w.pdb --ligand amoxicillin.mol2

Loading protein structure from example.pdb...
Running fpocket analysis...

Looking for fpocket output in: 2c6w_out/2c6w_info.txt
Found fpocket output file, parsing pockets...
Found 17 pockets
Analyzing protein flexibility...

Analyzing compatibility with ligand amoxicillin.mol2 at pocket 3...
(Auto-selected pocket 4 based on highest druggability score)

Compatibility analysis results:
pocket_volume: 3.829
pocket_score: 0.094
druggability_score: 0.665
binding_site_rank: 4
total_pockets_found: 17
```

## Error Handling

BindSite includes comprehensive error checking for:
- Missing input files
- fpocket installation and execution
- Invalid pocket IDs
- Structure parsing errors

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

Adam Davidson - adam.davidson128@gmail.com

Project Link: https://github.com/adam-davidson128/bindsite