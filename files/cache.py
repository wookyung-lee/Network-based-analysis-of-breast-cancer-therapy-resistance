"""
cache.py — Lightweight persistence helpers for the pipeline.

Each step saves its outputs to a subfolder under CACHE_DIR (default: .cache/).
The next step loads from there instead of re-running everything upstream.

Formats used
------------
DataFrames      → Parquet  (.parquet)   fast + preserves dtypes
NetworkX graphs → GraphML  (.graphml)   human-readable, lossless
Plain dicts/lists → JSON   (.json)
"""

import json
import os
import pickle

import networkx as nx
import pandas as pd

CACHE_DIR = ".cache"


def _path(filename: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, filename)


# ---------------------------------------------------------------------------
# DataFrames
# ---------------------------------------------------------------------------

def save_df(df: pd.DataFrame, filename: str) -> None:
    p = _path(filename)
    df.to_parquet(p, index=True)
    print(f"[cache] saved  → {p}")


def load_df(filename: str) -> pd.DataFrame:
    p = _path(filename)
    df = pd.read_parquet(p)
    print(f"[cache] loaded ← {p}")
    return df


# ---------------------------------------------------------------------------
# NetworkX graphs
# ---------------------------------------------------------------------------

def save_graph(G: nx.Graph, filename: str) -> None:
    """Save as pickle (preserves all edge/node attributes including floats)."""
    p = _path(filename)
    with open(p, "wb") as f:
        pickle.dump(G, f)
    print(f"[cache] saved  → {p}")


def load_graph(filename: str) -> nx.Graph:
    p = _path(filename)
    with open(p, "rb") as f:
        G = pickle.load(f)
    print(f"[cache] loaded ← {p}")
    return G


# ---------------------------------------------------------------------------
# JSON (dicts, lists)
# ---------------------------------------------------------------------------

def save_json(obj, filename: str) -> None:
    p = _path(filename)
    with open(p, "w") as f:
        json.dump(obj, f)
    print(f"[cache] saved  → {p}")


def load_json(filename: str):
    p = _path(filename)
    with open(p) as f:
        obj = json.load(f)
    print(f"[cache] loaded ← {p}")
    return obj


# ---------------------------------------------------------------------------
# Existence check
# ---------------------------------------------------------------------------

def exists(filename: str) -> bool:
    return os.path.exists(_path(filename))
