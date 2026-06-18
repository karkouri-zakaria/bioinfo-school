#!/usr/bin/env python3
"""Exercise B: Structure Prediction Confidence Analysis.

Downloads real AlphaFold PDB files from the AlphaFold Database (e.g. Ubiquitin,
Alpha-synuclein), parses the B-factor column to extract pLDDT confidence scores,
identifies structured vs disordered regions, and saves a line plot visualization.
"""

import os
import urllib.request
import matplotlib.pyplot as plt
import numpy as np

import json

def download_pdb(uniprot_id, output_path):
    # Query AlphaFold API to get PDB URL
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"Querying AlphaFold API for {uniprot_id}...")
    try:
        req = urllib.request.Request(api_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
        
        if not data or 'pdbUrl' not in data[0]:
            print(f"No entry found in AlphaFold DB for {uniprot_id}")
            return False
            
        pdb_url = data[0]['pdbUrl']
        print(f"Downloading PDB from {pdb_url}...")
        pdb_req = urllib.request.Request(pdb_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(pdb_req) as response:
            with open(output_path, 'wb') as f:
                f.write(response.read())
        print(f"Saved PDB to {output_path}")
        return True
    except Exception as e:
        print(f"Error downloading {uniprot_id}: {e}")
        return False

def parse_plddt_from_pdb(pdb_path):
    """Parses pLDDT from the B-factor column of CA atoms in the PDB file."""
    plddt_scores = {}
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA": # Alpha carbon representation
                    res_seq = int(line[22:26].strip())
                    plddt = float(line[60:66].strip())
                    plddt_scores[res_seq] = plddt
    # Sort by residue index
    sorted_residues = sorted(plddt_scores.keys())
    scores = [plddt_scores[r] for r in sorted_residues]
    return sorted_residues, scores

def analyze_protein(uniprot_id, name):
    pdb_path = f"/root/bioinfo-school/exercises/week3/{uniprot_id}.pdb"
    if not os.path.exists(pdb_path):
        success = download_pdb(uniprot_id, pdb_path)
        if not success:
            return None
            
    residues, scores = parse_plddt_from_pdb(pdb_path)
    if not scores:
        print(f"Error: No CA atoms found in PDB file {pdb_path}")
        return None
        
    scores = np.array(scores)
    mean_val = np.mean(scores)
    min_val = np.min(scores)
    max_val = np.max(scores)
    
    # Classify confidence
    very_high = np.sum(scores > 90) / len(scores) * 100
    confident = np.sum((scores > 70) & (scores <= 90)) / len(scores) * 100
    low = np.sum((scores > 50) & (scores <= 70)) / len(scores) * 100
    very_low = np.sum(scores <= 50) / len(scores) * 100
    
    print(f"\n--- AlphaFold Structure Confidence Report for {name} ({uniprot_id}) ---")
    print(f"Total residues: {len(residues)}")
    print(f"pLDDT Metrics: Mean = {mean_val:.2f}, Min = {min_val:.2f}, Max = {max_val:.2f}")
    print(f"Confidence Categories:")
    print(f"  Very High (>90 pLDDT): {very_high:.1f}%")
    print(f"  Confident (70-90 pLDDT): {confident:.1f}%")
    print(f"  Low (50-70 pLDDT): {low:.1f}%")
    print(f"  Very Low (<50 pLDDT, likely disordered): {very_low:.1f}%")
    
    # Identify disordered regions (continuous stretches of pLDDT <= 50)
    disordered_regions = []
    current_region = []
    for r, s in zip(residues, scores):
        if s <= 50:
            current_region.append(r)
        else:
            if len(current_region) >= 5: # report regions of length >= 5
                disordered_regions.append((current_region[0], current_region[-1]))
            current_region = []
    if len(current_region) >= 5:
        disordered_regions.append((current_region[0], current_region[-1]))
        
    if disordered_regions:
        print("Disordered / flexible regions identified (residues):")
        for start, end in disordered_regions:
            print(f"  Residues {start}-{end}")
    else:
        print("No significant disordered regions found.")
        
    return residues, scores, mean_val

def main():
    os.makedirs("/root/bioinfo-school/exercises/week3", exist_ok=True)
    
    proteins = [
        ("P0CG48", "Ubiquitin (highly structured core)"),
        ("P37840", "Alpha-synuclein (intrinsically disordered protein)")
    ]
    
    plt.figure(figsize=(10, 6))
    
    for uniprot_id, name in proteins:
        res = analyze_protein(uniprot_id, name)
        if res:
            residues, scores, mean_val = res
            plt.plot(residues, scores, label=f"{name} (Mean={mean_val:.1f})", linewidth=2)
            
    plt.axhline(90, color="green", linestyle="--", alpha=0.5, label="Very High (90)")
    plt.axhline(70, color="orange", linestyle="--", alpha=0.5, label="Confident (70)")
    plt.axhline(50, color="red", linestyle="--", alpha=0.5, label="Very Low (50)")
    
    plt.title("AlphaFold pLDDT Score Profile Comparison")
    plt.xlabel("Residue Position")
    plt.ylabel("pLDDT (Model Confidence)")
    plt.ylim(0, 105)
    plt.legend(loc="lower left")
    plt.grid(True, linestyle=":", alpha=0.5)
    
    plot_path = "/root/bioinfo-school/exercises/week3/structure_plddt.png"
    plt.savefig(plot_path, dpi=300)
    print(f"\npLDDT plot successfully saved at: {plot_path}")

if __name__ == "__main__":
    main()
