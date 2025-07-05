"""
Enviro-Mol – core risk calculator
================================
A modular class (`GreenAssess`) that computes the Enviro-Risk Index (ERI)
for an arbitrary molecule supplied as SMILES.  All heavy lifting—descriptor
calculation, parameter-level scoring and ERI aggregation—lives here so that
Streamlit (or any other frontend) can import it as a pure Python library.

File layout expectation
-----------------------
project_root/
├─ green_core/
│  ├─ green_assess.py   ← *you are here*
│  ├─ data/
│  │   ├─ fragments.json      # BIOWIN fragment weights
│  │   ├─ ecosar_coeff.yml    # ECOSAR (a,b) pairs per class
│  │   ├─ oh_groups.pkl       # ∑ k_i group rate constants
│  │   └─ gwp_ref.pkl         # IPCC α,τ lookup
│  └─ __init__.py             # optional
└─ app.py                     # Streamlit UI (separate file)
"""
from __future__ import annotations

import json
import math
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List
import yaml  # dosyanın başına

def _yaml_load(path: Path):
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


# ---------------------------------------------------------------------------
# Helper / utility functions
# ---------------------------------------------------------------------------

def _scale(value: float, low: float, high: float) -> float:
    """Linear scale *value* to 0–1 given a [low, high] window and clip."""
    return max(0.0, min(1.0, (value - low) / (high - low)))


def _json_load(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _pickle_load(path: Path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

@dataclass
class GreenAssess:
    """Compute parameter-level environmental scores and overall ERI."""

    smiles: str
    data_dir: Path = Path(__file__).with_suffix("").parent / "data"
    weights: Dict[str, float] = field(
        default_factory=lambda: {
            "biodeg": 0.25,
            "bcf": 0.25,
            "lc50": 0.20,
            "oh": 0.15,
            "koc": 0.10,
            "gwp": 0.05,
        }
    )

    # Populated after init
    mol: Chem.Mol = field(init=False)
    descriptors: Dict[str, float] = field(init=False)
    scores: Dict[str, float] = field(init=False)
    eri: float = field(init=False)

    # ---------------------------------------------------------------------
    # Life-cycle
    # ---------------------------------------------------------------------

    def __post_init__(self):
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError("Invalid SMILES string provided.")

        self.descriptors = self._compute_descriptors()
        self.scores = self._compute_param_scores()
        self.eri = self._compute_eri()

    # ------------------------------------------------------------------
    # Descriptor block (fast, RDKit)
    # ------------------------------------------------------------------

    def _compute_descriptors(self) -> Dict[str, float]:
        d: Dict[str, float] = {}
        d["logKow"] = Crippen.MolLogP(self.mol)
        d["mol_wt"] = Descriptors.MolWt(self.mol)
        # TODO: add fragment counters, TPSA, etc.
        return d

    # ------------------------------------------------------------------
    # Parameter scores (normalised 0–1)
    # ------------------------------------------------------------------

    # 1. Biodegradability ------------------------------------------------

    def _biodeg_score(self) -> float:
        frag_coeffs = _json_load(self.data_dir / "fragments.json")
        S = 0.0
        # --- count fragments present in molecule
        for frag_smarts, weight in frag_coeffs.items():
            patt = Chem.MolFromSmarts(frag_smarts)
            if self.mol.HasSubstructMatch(patt):
                S += weight * len(self.mol.GetSubstructMatches(patt))
        # BIOWIN-3 correlation
        t_half = 10 ** (2.17 - 0.48 * S)  # days
        return _scale(t_half, 15, 180)

    # 2. BCF --------------------------------------------------------------

    def _bcf_score(self) -> float:
        logKow = self.descriptors["logKow"]
        kow = 10 ** logKow
        k1 = 0.02 * kow ** 0.7
        k2 = 0.003 * kow ** -0.5
        logBCF = math.log10(k1 / k2)
        return _scale(logBCF, 3, 4.5)

    # 3. LC50 -------------------------------------------------------------

    def _lc50_score(self) -> float:
        ecosar = _json_load(self.data_dir / "ecosar_coeff.yml")  # class → (a,b)
        chem_class = self._assign_ecosar_class(ecosar.keys())
        a, b = ecosar[chem_class]
        logKow = self.descriptors["logKow"]
        log_lc50 = a * logKow + b  # mg/L
        return _scale(log_lc50, math.log10(100), math.log10(1))  # invert

    # 4. OH half-life ----------------------------------------------------

    def _oh_score(self) -> float:
        group_k = _pickle_load(self.data_dir / "oh_groups.pkl")  # dict
        k_total = 0.0
        for smarts, k_i in group_k.items():
            patt = Chem.MolFromSmarts(smarts)
            k_total += k_i * len(self.mol.GetSubstructMatches(patt))
        oh_conc = 1e6  # cm-3
        tau_h = 1 / (k_total * oh_conc) / 3600  # seconds→hours
        return _scale(tau_h, 12, 96)

    # 5. Koc --------------------------------------------------------------

    def _koc_score(self) -> float:
        logKow = self.descriptors["logKow"]
        log_koc = 0.54 * logKow + 1.47
        return _scale(log_koc, 2, 4)

    # 6. GWP penalty ------------------------------------------------------

    def _gwp_score(self) -> float:
        gwp_data = _pickle_load(self.data_dir / "gwp_ref.pkl")
        if any(atom.GetSymbol() in {"F", "Cl", "Br"} for atom in self.mol.GetAtoms()):
            # crude — refine with α,τ lookup if needed
            gwp = 1500  # placeholder high value
        else:
            gwp = 0
        return min(1.0, gwp / 1500)

    # ------------------------------------------------------------------
    # Aggregation utilities
    # ------------------------------------------------------------------

    def _compute_param_scores(self) -> Dict[str, float]:
        return {
            "biodeg": self._biodeg_score(),
            "bcf": self._bcf_score(),
            "lc50": self._lc50_score(),
            "oh": self._oh_score(),
            "koc": self._koc_score(),
            "gwp": self._gwp_score(),
        }

    def _compute_eri(self) -> float:
        s = self.scores
        d_pos = math.sqrt(sum((val - 1) ** 2 for val in s.values()))
        d_neg = math.sqrt(sum(val ** 2 for val in s.values()))
        return 100 * (1 - d_neg / (d_pos + d_neg))

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    def to_dataframe(self):
        """Return a pandas DataFrame (Param, Raw, Score)."""
        import pandas as pd
        rows: List[Dict[str, float]] = []
        for param, score in self.scores.items():
            rows.append({
                "Parameter": param,
                "Score (0-1)": round(score, 3),
            })
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------

    # --- stub for ECOSAR classification – to implement properly --------
    def _assign_ecosar_class(self, classes):
        # naive: pick first until we implement SMARTS mapping
        return next(iter(classes))


__all__ = ["GreenAssess"]
