"""
discriminate.py — Pharmacophore-based active/decoy discrimination.

Given a pocket_clusters directory (output of cluster_pockets) and two SDF files
(actives and decoys), scores each cluster representative by its ability to
discriminate actives from decoys using residue-based pharmacophore vectors and
cosine similarity. Outputs ranked conformations with ROC-AUC, EF1%, and EF5%.
"""

import os
import logging
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
from prody import parsePDB

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
_factory_instance = None

_logger = logging.getLogger(__name__)


def _get_factory():
    global _factory_instance
    if _factory_instance is None:
        if not os.path.exists(_FDEF_PATH):
            raise RuntimeError(
                f"RDKit feature definition file not found: {_FDEF_PATH}. "
                "Ensure rdkit data files are installed (e.g. conda install -c conda-forge rdkit)."
            )
        _factory_instance = ChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)
    return _factory_instance


def _extract_resname(residue_id: str, resname_map: dict = None) -> str:
    """
    Extract the 3-letter residue name from a p2rank residue ID.

    Tries lookup in resname_map first (for 'CHAIN_RESNUM' format like 'A_1019'),
    then falls back to parsing the ID string (for 'CHAIN_NUM_RESNAME' format like 'A_1_ALA').
    """
    rid = residue_id.strip()

    # Try direct lookup in map (handles 'A_1019' format from p2rank)
    if resname_map and rid in resname_map:
        return resname_map[rid]

    # Fall back to string parsing (handles 'A_1_ALA' format)
    parts = rid.split('_')
    for part in reversed(parts):
        clean = part.strip().upper()
        if len(clean) == 3 and clean.isalpha():
            return clean

    # Last resort: first 3 alpha chars
    alpha = ''.join(c for c in rid if c.isalpha())
    result = alpha[:3].upper() if len(alpha) >= 3 else ''
    if result:
        _logger.warning(
            "Residue ID '%s' did not match expected format; extracted '%s' via fallback heuristic.",
            rid, result
        )
    return result


def build_resname_map(pdb_path: str) -> dict:
    """
    Build a mapping from 'CHAIN_RESNUM' to 3-letter residue name by parsing a PDB file.

    Uses prody to read the structure. Returns a dict like {'A_1019': 'ARG', 'A_1022': 'SER', ...}.
    Returns an empty dict if the file cannot be parsed.
    """
    try:
        structure = parsePDB(pdb_path, subset='ca', verbosity='none')
        if structure is None:
            _logger.warning("prody could not parse PDB file: %s", pdb_path)
            return {}
        resname_map = {}
        for res in structure.iterResidues():
            key = f"{res.getChid()}_{res.getResnum()}"
            resname_map[key] = res.getResname()
        return resname_map
    except Exception as exc:
        _logger.warning("Failed to build resname map from %s: %s", pdb_path, exc)
        return {}


def get_pocket_features(residue_ids: list, resname_map: dict = None) -> dict:
    """
    Compute pharmacophore feature counts for a pocket defined by residue IDs.

    Args:
        residue_ids: List of p2rank residue ID strings (e.g. ['A_1019', 'A_1022']
                     or ['A_1_ALA', 'A_2_ARG']).
        resname_map: Optional dict mapping 'CHAIN_RESNUM' to 3-letter residue name,
                     built from the representative PDB via build_resname_map().
                     Required when residue IDs use 'CHAIN_RESNUM' format without
                     the residue name suffix.

    Returns:
        Dict with keys matching FEATURE_KEYS and integer counts.
    """
    features = {k: 0 for k in FEATURE_KEYS}
    for rid in residue_ids:
        resname = _extract_resname(rid, resname_map)
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
        feats = _get_factory().GetFeaturesForMol(mol)
        for feat in feats:
            family = feat.GetFamily()
            key = _RDKIT_FEATURE_MAP.get(family)
            if key:
                features[key] += 1
    except Exception as exc:
        _logger.warning("Failed to extract features from molecule: %s", exc)
    return features


def score_complementarity(pocket_features: dict, ligand_features: dict) -> float:
    """
    Compute pharmacophoric complementarity between a pocket and a ligand.

    H-bond and charge features are cross-mapped (pocket donor pairs with ligand
    acceptor, pocket acceptor pairs with ligand donor, etc.). Hydrophobic and
    aromatic features are self-complementary (like pairs with like).

    Returns a float in [0, 1]: 1 means perfect complementarity, 0 means no
    complementary features. Returns 0.0 if either vector is all zeros.
    """
    # Remap pocket features to their complementary counterparts
    complementary_pocket = {
        'donor':       pocket_features['acceptor'],    # pocket acceptor pairs with ligand donor
        'acceptor':    pocket_features['donor'],       # pocket donor pairs with ligand acceptor
        'hydrophobic': pocket_features['hydrophobic'], # self-complementary
        'aromatic':    pocket_features['aromatic'],    # self-complementary (pi-stacking)
        'positive':    pocket_features['negative'],    # pocket negative pairs with ligand positive
        'negative':    pocket_features['positive'],    # pocket positive pairs with ligand negative
    }
    p = np.array([complementary_pocket[k] for k in FEATURE_KEYS], dtype=float)
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
                       outfolder: str, pdb_dir: str = None) -> pd.DataFrame:
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

        # Build residue name map from representative PDB if pdb_dir provided
        resname_map = {}
        if pdb_dir:
            file_name = str(rep.get('File name', ''))
            # Strip p2rank prediction suffix to get PDB filename
            # e.g. 'trajectory_test_data_30.pdb_predictions' -> 'trajectory_test_data_30.pdb'
            pdb_name = file_name.replace('_predictions', '')
            pdb_path = os.path.join(pdb_dir, pdb_name)
            if os.path.exists(pdb_path):
                resname_map = build_resname_map(pdb_path)
            else:
                _logger.warning("Representative PDB not found: %s", pdb_path)

        pocket_feats = get_pocket_features(residue_ids, resname_map=resname_map)
        scores_actives = [score_complementarity(pocket_feats, lf) for lf in active_features]
        scores_decoys = [score_complementarity(pocket_feats, decoy_feat) for decoy_feat in decoy_features]
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
