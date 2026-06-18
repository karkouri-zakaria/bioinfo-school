#!/usr/bin/env python3
"""Exercise A: Protein Embeddings with ESM2.

Loads proteins from proteins.fasta, extracts embeddings using a small ESM2 model,
compares mean pooling vs CLS-token pooling, and plots pairwise cosine similarity matrices.
"""

import os
import torch
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModel
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA

def main():
    os.makedirs("/root/bioinfo-school/exercises/week3", exist_ok=True)
    
    fasta_path = "/root/bioinfo-school/exercises/week3/proteins.fasta"
    if not os.path.exists(fasta_path):
        print(f"Error: fasta file {fasta_path} does not exist. Run get_proteins.py first.")
        return
        
    # Load sequences and labels
    records = list(SeqIO.parse(fasta_path, "fasta"))
    sequences = []
    labels = []
    accessions = []
    
    for rec in records:
        sequences.append(str(rec.seq))
        accessions.append(rec.id)
        # Extract family from header, e.g. [Family=GPCR]
        family = "Unknown"
        if "[Family=" in rec.description:
            family = rec.description.split("[Family=")[1].split("]")[0]
        labels.append(family)
        
    print(f"Loaded {len(sequences)} sequences from FASTA.")
    
    # Load ESM2 model (T6, 8M params, lightweight for CPU)
    model_name = "facebook/esm2_t6_8M_UR50D"
    print(f"Loading tokenizer and model: {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    
    mean_embeddings = []
    cls_embeddings = []
    
    print("Extracting embeddings...")
    with torch.no_grad():
        for i, seq in enumerate(sequences):
            inputs = tokenizer(seq, return_tensors="pt")
            outputs = model(**inputs)
            # outputs.last_hidden_state shape: (1, seq_len + 2, embed_dim) (includes CLS and EOS tokens)
            hidden_states = outputs.last_hidden_state[0] # (seq_len + 2, embed_dim)
            
            # 1. Mean Pooling (excluding CLS/EOS)
            # The CLS is at index 0, EOS is at index seq_len+1
            residue_states = hidden_states[1:-1, :]
            mean_emb = residue_states.mean(dim=0).numpy()
            mean_embeddings.append(mean_emb)
            
            # 2. CLS Token Pooling (suboptimal representation)
            cls_emb = hidden_states[0, :].numpy()
            cls_embeddings.append(cls_emb)
            
            if (i + 1) % 10 == 0:
                print(f"  Embedded {i + 1}/{len(sequences)} sequences...")
                
    mean_embeddings = np.array(mean_embeddings)
    cls_embeddings = np.array(cls_embeddings)
    
    # Sort by labels to group families together in similarity matrices
    sorted_indices = np.argsort(labels)
    sorted_labels = [labels[idx] for idx in sorted_indices]
    sorted_accessions = [accessions[idx] for idx in sorted_indices]
    
    mean_embeddings_sorted = mean_embeddings[sorted_indices]
    cls_embeddings_sorted = cls_embeddings[sorted_indices]
    
    # Pairwise similarity
    similarity_mean = cosine_similarity(mean_embeddings_sorted)
    similarity_cls = cosine_similarity(cls_embeddings_sorted)
    
    # Calculate quantitative cluster separation metric:
    # average similarity of same-family pairs vs different-family pairs
    def evaluate_clustering(sim_matrix, sorted_lbls):
        same_family = []
        diff_family = []
        n = len(sorted_lbls)
        for r in range(n):
            for c in range(r + 1, n):
                if sorted_lbls[r] == sorted_lbls[c]:
                    same_family.append(sim_matrix[r, c])
                else:
                    diff_family.append(sim_matrix[r, c])
        return np.mean(same_family), np.mean(diff_family)
        
    mean_same, mean_diff = evaluate_clustering(similarity_mean, sorted_labels)
    cls_same, cls_diff = evaluate_clustering(similarity_cls, sorted_labels)
    
    print("\n--- Embedding Metric Evaluation ---")
    print(f"Mean Pooling: Within-family similarity = {mean_same:.4f}, Between-family similarity = {mean_diff:.4f}")
    print(f"  Ratio (higher is better cluster separation) = {mean_same / mean_diff:.4f}")
    print(f"CLS Pooling:  Within-family similarity = {cls_same:.4f}, Between-family similarity = {cls_diff:.4f}")
    print(f"  Ratio = {cls_same / cls_diff:.4f}")
    
    # Plotting Heatmaps and PCA plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Heatmap - Mean Pooling
    im1 = axes[0, 0].imshow(similarity_mean, cmap="coolwarm", vmin=0.4, vmax=1.0)
    axes[0, 0].set_title(f"Mean Pooling Cosine Similarity\nWithin/Between Family Ratio: {mean_same/mean_diff:.3f}")
    fig.colorbar(im1, ax=axes[0, 0])
    
    # Heatmap - CLS Pooling
    im2 = axes[0, 1].imshow(similarity_cls, cmap="coolwarm", vmin=0.4, vmax=1.0)
    axes[0, 1].set_title(f"CLS Pooling Cosine Similarity\nWithin/Between Family Ratio: {cls_same/cls_diff:.3f}")
    fig.colorbar(im2, ax=axes[0, 1])
    
    # PCA - Mean Pooling
    pca_mean = PCA(n_components=2).fit_transform(mean_embeddings)
    # PCA - CLS Pooling
    pca_cls = PCA(n_components=2).fit_transform(cls_embeddings)
    
    unique_labels = sorted(list(set(labels)))
    colors = ['#ef4444', '#3b82f6', '#10b981'] # red, blue, green
    
    for label, color in zip(unique_labels, colors):
        mask = np.array(labels) == label
        # PCA Mean
        axes[1, 0].scatter(pca_mean[mask, 0], pca_mean[mask, 1], label=label, c=color, alpha=0.8, edgecolors='none')
        # PCA CLS
        axes[1, 1].scatter(pca_cls[mask, 0], pca_cls[mask, 1], label=label, c=color, alpha=0.8, edgecolors='none')
        
    axes[1, 0].set_title("PCA: Mean Pooling")
    axes[1, 0].legend()
    axes[1, 0].grid(True, linestyle="--", alpha=0.3)
    
    axes[1, 1].set_title("PCA: CLS Pooling")
    axes[1, 1].legend()
    axes[1, 1].grid(True, linestyle="--", alpha=0.3)
    
    plt.tight_layout()
    plot_path = "/root/bioinfo-school/exercises/week3/family_clustering.png"
    plt.savefig(plot_path, dpi=300)
    print(f"\nClustering plot successfully saved at: {plot_path}")

if __name__ == "__main__":
    main()
