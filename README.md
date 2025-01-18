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

```bash
# First install fpocket
## On Ubuntu/Debian:
sudo apt-get install fpocket

## On macOS with Homebrew:
brew install fpocket

# Then install BindSite
git clone https://github.com/adam-davidson128/bindsite.git
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
- `--pocket`: ID of pocket to analyze (default: 0)

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
$ python bindsite.py --pdb example.pdb --ligand drug.mol2 --pocket 0

Loading protein structure from example.pdb...
Running fpocket analysis...
Analyzing protein flexibility...

Analyzing compatibility with ligand drug.mol2 at pocket 0...

Compatibility analysis results:
pocket_volume: 235.6
pocket_score: 28.7
binding_site_rank: 1
total_pockets_found: 12
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

[Your Name] - your.email@example.com

Project Link: https://github.com/adam-davidson128/bindsite