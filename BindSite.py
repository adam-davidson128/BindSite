import os
import numpy as np
from Bio import PDB
from Bio.PDB import *
import subprocess
import pandas as pd
from scipy import spatial
from scipy.spatial import ConvexHull
from scipy.linalg import eigh
import prody
import mdtraj as md

class ProteinStructureAnalyzer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.structure = None
        self.pockets = []
        self.flexibility_data = None
        
    def load_structure(self):
        """Load PDB structure using BioPython"""
        try:
            parser = PDB.PDBParser(QUIET=True)
            self.structure = parser.get_structure('protein', self.pdb_file)
            return self.structure
        except Exception as e:
            raise ValueError(f"Failed to parse PDB file: {e}")
        
    def run_fpocket(self, min_alpha_sphere=3.0, max_alpha_sphere=6.0):
        """Run fpocket for binding site detection"""
        try:
            # Check if fpocket is installed
            if subprocess.run(['which', 'fpocket'], capture_output=True).returncode != 0:
                raise RuntimeError("fpocket is not installed")
            
            # Get absolute path of input file
            abs_pdb_path = os.path.abspath(self.pdb_file)
            print(f"Running fpocket on: {abs_pdb_path}")
            
            # Run fpocket with full output
            cmd = f"fpocket -f {abs_pdb_path}"
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            print("fpocket stdout:", result.stdout)
            if result.stderr:
                print("fpocket stderr:", result.stderr)
                
            self.parse_fpocket_results()
        except subprocess.CalledProcessError as e:
            print(f"Error running fpocket: {e}")
            if e.output:
                print("fpocket output:", e.output)
            if e.stderr:
                print("fpocket error output:", e.stderr)
            
    def parse_fpocket_results(self):
        """Parse fpocket output to get pocket information"""
        base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
        fpocket_dir = base_name + "_out"
        info_file = os.path.join(fpocket_dir, f"{base_name}_info.txt")
        
        print(f"Looking for fpocket output in: {info_file}")
        
        if not os.path.exists(info_file):
            raise FileNotFoundError(f"fpocket output not found at: {info_file}")
        
        self.pockets = []
        with open(info_file, 'r') as f:
            print("Found fpocket output file, parsing pockets...")
            lines = f.readlines()
            current_pocket = None
            for line in lines:
                if line.startswith("Pocket"):
                    if current_pocket:
                        self.pockets.append(current_pocket)
                    current_pocket = {
                        'atoms': [],
                        'score': None,
                        'volume': None
                    }
                elif current_pocket and line.strip():
                    # Parse pocket details
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        if "Score" in line:
                            try:
                                current_pocket['score'] = float(parts[-1])
                            except ValueError:
                                pass
                        elif "Volume" in line:
                            try:
                                current_pocket['volume'] = float(parts[-1])
                            except ValueError:
                                pass
            
            if current_pocket:
                self.pockets.append(current_pocket)
        
        print(f"Found {len(self.pockets)} pockets")
        
    def analyze_pocket(self, pocket_id):
        """Analyze specific pocket characteristics"""
        if not self.pockets or pocket_id >= len(self.pockets):
            raise ValueError(f"Invalid pocket ID: {pocket_id}")
        
        pocket = self.pockets[pocket_id]
        return {
            'volume': pocket.get('volume', None),
            'score': pocket.get('score', None)
        }
        
    def extract_b_factors(self):
        """Extract and normalize B-factors"""
        b_factors = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        b_factors.append(atom.get_bfactor())
        return np.array(b_factors)
        
    def analyze_secondary_structure(self):
        """Analyze secondary structure elements"""
        try:
            dssp = PDB.DSSP(self.structure[0], self.pdb_file)
            ss_data = {}
            for residue in dssp:
                chain_id = residue[0]
                ss_type = residue[2]
                if chain_id not in ss_data:
                    ss_data[chain_id] = []
                ss_data[chain_id].append(ss_type)
            return ss_data
        except Exception as e:
            print(f"Warning: DSSP analysis failed: {e}")
            return None

    def analyze_flexibility(self):
        """Analyze protein flexibility using multiple methods"""
        self.flexibility_data = {
            'b_factors': self.extract_b_factors(),
            'secondary_structure': self.analyze_secondary_structure()
        }
        return self.flexibility_data

class ChemicalCompatibilityAnalyzer:
    def __init__(self, protein_analyzer):
        self.protein_analyzer = protein_analyzer
        
    def analyze_ligand_compatibility(self, ligand, pocket_id):
        """Analyze chemical compatibility between ligand and pocket"""
        if not os.path.exists(ligand):
            raise FileNotFoundError(f"Ligand file not found: {ligand}")
            
        pocket_analysis = self.protein_analyzer.analyze_pocket(pocket_id)
        
        try:
            # Basic compatibility check
            compatibility = {
                'pocket_volume': pocket_analysis['volume'],
                'pocket_score': pocket_analysis['score'],
                'binding_site_rank': pocket_id + 1,  # Convert to 1-based indexing for display
                'total_pockets_found': len(self.protein_analyzer.pockets)
            }
            
            return compatibility
            
        except Exception as e:
            raise ValueError(f"Failed to analyze compatibility: {e}")

def main():
    import argparse
    
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Analyze protein-ligand binding compatibility')
    parser.add_argument('--pdb', type=str, default="protein.pdb",
                      help='PDB file for protein structure (default: protein.pdb)')
    parser.add_argument('--ligand', type=str, default="ligand.mol2",
                      help='MOL2 file for ligand structure (default: ligand.mol2)')
    parser.add_argument('--pocket', type=int, default=0,
                      help='Pocket ID to analyze (default: 0)')
    
    args = parser.parse_args()
    
    try:
        print(f"Loading protein structure from {args.pdb}...")
        if not os.path.exists(args.pdb):
            print(f"Error: PDB file {args.pdb} not found")
            return
        protein_analyzer = ProteinStructureAnalyzer(args.pdb)
        protein_analyzer.load_structure()
        
        print("Running fpocket analysis...")
        protein_analyzer.run_fpocket()
        
        print("Analyzing protein flexibility...")
        protein_analyzer.analyze_flexibility()
        
        # Initialize chemical compatibility analysis
        compatibility_analyzer = ChemicalCompatibilityAnalyzer(protein_analyzer)
        
        if not os.path.exists(args.ligand):
            raise FileNotFoundError(f"Ligand file {args.ligand} not found")
            
        if args.pocket >= len(protein_analyzer.pockets):
            raise ValueError(f"Invalid pocket ID: {args.pocket}")
            
        print(f"\nAnalyzing compatibility with ligand {args.ligand} at pocket {args.pocket}...")
        compatibility = compatibility_analyzer.analyze_ligand_compatibility(args.ligand, args.pocket)
        print("\nCompatibility analysis results:")
        for key, value in compatibility.items():
            print(f"{key}: {value}")
        
    except Exception as e:
        print(f"Error in analysis: {e}")

if __name__ == "__main__":
    main()