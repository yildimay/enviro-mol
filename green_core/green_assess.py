"""
Enviro-Mol – core risk calculator
=================================
* Güncel sürüm: YAML loader, ECOSAR katsayılarını float’a çevirme, daha güvenli
  SMILES/fragment tarama ve minimal yorum satırları eklendi.
"""

from __future__ import annotations

import json
import math
import pickle
import yaml
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

# ---------------------------------------------------------------------------
# Helper / utility functions
# ---------------------------------------------------------------------------


def _scale(value: float, low: float, high: float) -> float:
    """Lineer ölçek + clip: [low, high] aralığını 0–1’e indirger."""
    return max(0.0, min(1.0, (value - low) / (high - low)))


def _json_load(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _yaml_load(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


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

    mol: Chem.Mol = field(init=False)
    descriptors: Dict[str, float] = field(init=False)
    scores: Dict[str, float] = field(init=False)
    eri: float = field(init=False)

    # ------------------------------------------------------------------ #
    # Life-cycle                                                         #
    # ------------------------------------------------------------------ #

    def __post_init__(self):
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError("Invalid SMILES string provided.")

        self.descriptors = self._compute_descriptors()
        self.scores = self._compute_param_scores()
        self.eri = self._compute_eri()

    # ------------------------------------------------------------------ #
    # Descriptor block                                                   #
    # ------------------------------------------------------------------ #

    def _compute_descriptors(self) -> Dict[str, float]:
        d: Dict[str, float] = {}
        d["logKow"] = Crippen.MolLogP(self.mol)
        d["mol_wt"] = Descriptors.MolWt(self.mol)
        # Gerekirse ek deskriptor ekle (TPSA, fragman sayıları vb.)
        return d

    # ------------------------------------------------------------------ #
    # Individual parameter scores                                        #
    # ------------------------------------------------------------------ #

    # 1 | Biodegradability
    def _biodeg_score(self) -> float:
        frag_coeffs = _json_load(self.data_dir / "fragments.json")
        S = 0.0
        for frag_smarts, weight in frag_coeffs.items():
            patt = Chem.MolFromSmarts(frag_smarts)
            if patt and self.mol.HasSubstructMatch(patt):
                S += weight * len(self.mol.GetSubstructMatches(patt))
        t_half = 10 ** (2.17 - 0.48 * S)  # gün
        return _scale(t_half, 15, 180)

    # 2 | BCF
    def _bcf_score(self) -> float:
        logKow = self.descriptors["logKow"]
        kow = 10 ** logKow
        k1 = 0.02 * kow ** 0.7
        k2 = 0.003 * kow ** -0.5
        logBCF = math.log10(k1 / k2)
        return _scale(logBCF, 3, 4.5)

    # 3 | Balık LC50
    def _lc50_score(self) -> float:
        ecosar = _yaml_load(self.data_dir / "ecosar_coeff.yml")
        chem_class = self._assign_ecosar_class(ecosar.keys())
        coeff = ecosar.get(chem_class)
        if coeff is None:
            raise ValueError(f"ECOSAR class '{chem_class}' not found in coeff file")
        a = float(coeff["a"])
        b = float(coeff["b"])
        logKow = self.descriptors["logKow"]
        log_lc50 = a * logKow + b
        return _scale(log_lc50, math.log10(100), math.log10(1))  # ters skala

    # 4 | Atmosferik OH reaktivitesi
    def _oh_score(self) -> float:
        group_k = _pickle_load(self.data_dir / "oh_groups.pkl")
        k_total = 0.0
        for smarts, k_i in group_k.items():
            patt = Chem.MolFromSmarts(smarts)
            if patt:
                k_total += k_i * len(self.mol.GetSubstructMatches(patt))
        oh_conc = 1e6  # cm⁻³
        tau_h = 1 / (k_total * oh_conc) / 3600  # saniye → saat
        return _scale(tau_h, 12, 96)

    # 5 | Toprak mobilitesi (Koc)
    def _koc_score(self) -> float:
        logKow = self.descriptors["logKow"]
        log_koc = 0.54 * logKow + 1.47
        return _scale(log_koc, 2, 4)

    # 6 | GWP / Halojen cezası
    def _gwp_score(self) -> float:
        # TODO: tam α-τ tablosu entegrasyonu
        halogen = any(atom.GetSymbol() in {"F", "Cl", "Br"} for atom in self.mol.GetAtoms())
        gwp = 1500 if halogen else 0
        return min(1.0, gwp / 1500)

    # ------------------------------------------------------------------ #
    # Aggregation                                                        #
    # ------------------------------------------------------------------ #

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
        d_pos = math.sqrt(sum((v - 1) ** 2 for v in self.scores.values()))
        d_neg = math.sqrt(sum(v ** 2 for v in self.scores.values()))
        return 100 * (1 - d_neg / (d_pos + d_neg))

    # ------------------------------------------------------------------ #
    # Public helper                                                      #
    # ------------------------------------------------------------------ #

    def to_dataframe(self):
        import pandas as pd

        return pd.DataFrame(
            {
                "Parameter": list(self.scores.keys()),
                "Score (0-1)": [round(val, 3) for val in self.scores.values()],
            }
        )

    # ------------------------------------------------------------------ #
    # Internal util                                                      #
    # ------------------------------------------------------------------ #

    def _assign_ecosar_class(self, classes):
        # TODO SMARTS-eşleştirme; şimdilik ilk sınıf
        return next(iter(classes))


__all__ = ["GreenAssess"]
