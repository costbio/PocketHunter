"""
discriminate.py — Pharmacophore-based active/decoy discrimination.

Given a pocket_clusters directory (output of cluster_pockets) and two SDF files
(actives and decoys), scores each cluster representative by its ability to
discriminate actives from decoys using residue-based pharmacophore vectors and
cosine similarity. Outputs ranked conformations with ROC-AUC, EF1%, and EF5%.
"""

import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures

FEATURE_KEYS = ['donor', 'acceptor', 'hydrophobic', 'aromatic', 'positive', 'negative']

_DONOR_RESNAMES = {'SER', 'THR', 'TYR', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'TRP'}
_ACCEPTOR_RESNAMES = {'ASP', 'GLU', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'HIS'}
_HYDROPHOBIC_RESNAMES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}
_AROMATIC_RESNAMES = {'PHE', 'TRP', 'TYR', 'HIS'}
_POSITIVE_RESNAMES = {'LYS', 'ARG', 'HIS'}
_NEGATIVE_RESNAMES = {'ASP', 'GLU'}

_RDKIT_FEATURE_MAP = {
    'Donor': 'donor',
    'Acceptor': 'acceptor',
    'Hydrophobe': 'hydrophobic',
    'LumpedHydrophobe': 'hydrophobic',
    'Aromatic': 'aromatic',
    'PosIonizable': 'positive',
    'NegIonizable': 'negative',
}

_FDEF_PATH = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
_FACTORY = ChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)


def _extract_resname(residue_id: str) -> str:
    parts = residue_id.strip().split('_')
    for part in reversed(parts):
        clean = part.strip().upper()
        if len(clean) == 3 and clean.isalpha():
            return clean
    alpha = ''.join(c for c in residue_id if c.isalpha())
    return alpha[:3].upper() if len(alpha) >= 3 else ''


def get_pocket_features(residue_ids: list) -> dict:
    features = {k: 0 for k in FEATURE_KEYS}
    for rid in residue_ids:
        resname = _extract_resname(rid)
        if not resname:
            continue
        if resname in _DONOR_RESNAMES:
            features['donor'] += 1
        if resname in _ACCEPTOR_RESNAMES:
            features['acceptor'] += 1
        if resname in _HYDROPHOBIC_RESNAMES:
            features['hydrophobic'] += 1
        if resname in _AROMATIC_RESNAMES:
            features['aromatic'] += 1
        if resname in _POSITIVE_RESNAMES:
            features['positive'] += 1
        if resname in _NEGATIVE_RESNAMES:
            features['negative'] += 1
    return features


def get_ligand_features(mol) -> dict:
    features = {k: 0 for k in FEATURE_KEYS}
    if mol is None:
        return features
    try:
        feats = _FACTORY.GetFeaturesForMol(mol)
        for feat in feats:
            family = feat.GetFamily()
            key = _RDKIT_FEATURE_MAP.get(family)
            if key:
                features[key] += 1
    except Exception:
        pass
    return features


def score_complementarity(pocket_features: dict, ligand_features: dict) -> float:
    p = np.array([pocket_features[k] for k in FEATURE_KEYS], dtype=float)
    l = np.array([ligand_features[k] for k in FEATURE_KEYS], dtype=float)
    norm_p = np.linalg.norm(p)
    norm_l = np.linalg.norm(l)
    if norm_p == 0 or norm_l == 0:
        return 0.0
    return float(np.dot(p, l) / (norm_p * norm_l))


def _enrichment_factor(scores_actives: list, scores_decoys: list, fraction: float) -> float:
    n_actives = len(scores_actives)
    n_total = n_actives + len(scores_decoys)
    n_top = max(1, int(fraction * n_total))
    all_scores = [(s, 1) for s in scores_actives] + [(s, 0) for s in scores_decoys]
    all_scores.sort(key=lambda x: -x[0])
    actives_in_top = sum(label for _, label in all_scores[:n_top])
    expected = fraction * n_actives
    if expected == 0:
        return 0.0
    return float(actives_in_top / expected)


def compute_discrimination_metrics(scores_actives: list, scores_decoys: list) -> dict:
    labels = [1] * len(scores_actives) + [0] * len(scores_decoys)
    scores = list(scores_actives) + list(scores_decoys)
    try:
        roc_auc = float(roc_auc_score(labels, scores))
    except ValueError:
        roc_auc = 0.5
    ef1 = _enrichment_factor(scores_actives, scores_decoys, 0.01)
    ef5 = _enrichment_factor(scores_actives, scores_decoys, 0.05)
    return {'roc_auc': roc_auc, 'ef1': ef1, 'ef5': ef5}


def _load_mols_from_sdf(sdf_path: str) -> list:
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    return [mol for mol in supplier if mol is not None]


def run_discrimination(cluster_dir: str, actives_sdf: str, decoys_sdf: str,
                       outfolder: str) -> pd.DataFrame:
    os.makedirs(outfolder, exist_ok=True)

    reps_csv = os.path.join(cluster_dir, 'cluster_representatives.csv')
    if not os.path.exists(reps_csv):
        raise FileNotFoundError(f"cluster_representatives.csv not found in {cluster_dir}")

    df_reps = pd.read_csv(reps_csv)
    required_cols = {'residues', 'cluster', 'Frame'}
    missing_cols = required_cols - set(df_reps.columns)
    if missing_cols:
        raise ValueError(
            f"cluster_representatives.csv is missing required columns: {missing_cols}. "
            f"Found columns: {list(df_reps.columns)}"
        )

    actives = _load_mols_from_sdf(actives_sdf)
    decoys = _load_mols_from_sdf(decoys_sdf)

    if not actives:
        raise ValueError(f"No valid molecules loaded from actives SDF: {actives_sdf}")
    if not decoys:
        raise ValueError(f"No valid molecules loaded from decoys SDF: {decoys_sdf}")

    if len(actives) > 200:
        raise ValueError(f"Too many actives: {len(actives)} (max 200)")
    if len(decoys) > 2000:
        raise ValueError(f"Too many decoys: {len(decoys)} (max 2000)")

    active_features = [get_ligand_features(mol) for mol in actives]
    decoy_features = [get_ligand_features(mol) for mol in decoys]

    rows = []
    for _, rep in df_reps.iterrows():
        residue_ids = [r for r in str(rep['residues']).split() if r]
        pocket_feats = get_pocket_features(residue_ids)
        scores_actives = [score_complementarity(pocket_feats, lf) for lf in active_features]
        scores_decoys = [score_complementarity(pocket_feats, df) for df in decoy_features]
        metrics = compute_discrimination_metrics(scores_actives, scores_decoys)
        rows.append({
            'cluster_id': rep.get('cluster', ''),
            'frame': rep.get('Frame', ''),
            'roc_auc': round(metrics['roc_auc'], 4),
            'ef1': round(metrics['ef1'], 4),
            'ef5': round(metrics['ef5'], 4),
            'residues': rep['residues'],
            'probability': rep.get('probability', ''),
        })

    df_results = pd.DataFrame(rows).sort_values('roc_auc', ascending=False).reset_index(drop=True)
    df_results.to_csv(os.path.join(outfolder, 'discrimination_results.csv'), index=False)
    return df_results
