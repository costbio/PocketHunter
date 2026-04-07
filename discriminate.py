"""
discriminate.py — Pharmacophore-based active/decoy discrimination.

Given a pocket_clusters directory (output of cluster_pockets) and two SDF files
(actives and decoys), scores each cluster representative by its ability to
discriminate actives from decoys using residue-based pharmacophore vectors and
cosine similarity. Outputs ranked conformations with ROC-AUC, EF1%, and EF5%.

Scoring layers
--------------
Layer 0 (baseline):  Binary residue counting + cosine complementarity.
Layer 1 (+prob):     p2rank per-residue binding probabilities replace binary
                     weights — core residues contribute more than edge residues.
Layer 2 (+3D):       Rotation-invariant pairwise inter-feature distance matching
                     via a Gaussian kernel. Distances between complementary
                     pharmacophore features inside the pocket are compared to
                     distances between the matching feature types inside the ligand.
                     Combined score = 0.3 × cosine + 0.7 × 3D kernel.

Both layers activate automatically when the required data are present (p2rank
residues CSV for Layer 1; representative PDB with Cα coordinates for Layer 2).
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

# For each pocket feature type, the complementary ligand feature type
_COMPLEMENTARY = {
    'donor': 'acceptor',
    'acceptor': 'donor',
    'hydrophobic': 'hydrophobic',
    'aromatic': 'aromatic',
    'positive': 'negative',
    'negative': 'positive',
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


def build_ca_coord_map(pdb_path: str) -> dict:
    """
    Build a mapping from 'CHAIN_RESNUM' to Cα coordinates by parsing a PDB file.

    Returns a dict like {'A_1019': np.array([x, y, z]), ...}.
    Returns an empty dict if the file cannot be parsed.
    """
    try:
        structure = parsePDB(pdb_path, subset='ca', verbosity='none')
        if structure is None:
            _logger.warning("prody could not parse PDB file: %s", pdb_path)
            return {}
        coord_map = {}
        for res in structure.iterResidues():
            key = f"{res.getChid()}_{res.getResnum()}"
            coords = res.getCoords()
            if coords is not None and len(coords) > 0:
                coord_map[key] = coords[0]  # Cα is the only atom in 'ca' subset
        return coord_map
    except Exception as exc:
        _logger.warning("Failed to build Cα coord map from %s: %s", pdb_path, exc)
        return {}


def load_residue_probs(residues_csv_path: str) -> dict:
    """
    Load per-residue binding probabilities from a p2rank residues CSV file.

    Columns expected: chain, residue_label, residue_name, score, zscore, probability, pocket.
    Returns a dict mapping 'CHAIN_RESIDUELABEL' → probability (float, 0.0–1.0).
    Returns an empty dict if the file cannot be found or parsed.
    """
    if not residues_csv_path or not os.path.exists(residues_csv_path):
        return {}
    try:
        df = pd.read_csv(residues_csv_path)
        df.columns = [c.strip() for c in df.columns]
        result = {}
        for _, row in df.iterrows():
            chain = str(row.get('chain', '')).strip()
            label = str(row.get('residue_label', '')).strip()
            prob = float(row.get('probability', 0.0))
            key = f"{chain}_{label}"
            result[key] = prob
        return result
    except Exception as exc:
        _logger.warning("Failed to load residue probs from %s: %s", residues_csv_path, exc)
        return {}


def get_pocket_features(residue_ids: list, resname_map: dict = None,
                        residue_probs: dict = None) -> dict:
    """
    Compute pharmacophore feature counts (or probability-weighted sums) for a pocket.

    Args:
        residue_ids: List of p2rank residue ID strings (e.g. ['A_1019', 'A_1022']
                     or ['A_1_ALA', 'A_2_ARG']).
        resname_map: Optional dict mapping 'CHAIN_RESNUM' to 3-letter residue name,
                     built from the representative PDB via build_resname_map().
                     Required when residue IDs use 'CHAIN_RESNUM' format without
                     the residue name suffix.
        residue_probs: Optional dict mapping 'CHAIN_RESNUM' to per-residue p2rank
                       binding probability (0.0–1.0). When provided, each residue
                       contributes proportionally to its probability rather than
                       binary 0/1. Residues absent from the map fall back to 1.0.

    Returns:
        Dict with keys matching FEATURE_KEYS and float values.
    """
    features = {k: 0.0 for k in FEATURE_KEYS}
    for rid in residue_ids:
        resname = _extract_resname(rid, resname_map)
        if not resname:
            continue
        weight = residue_probs.get(rid, 1.0) if residue_probs is not None else 1.0
        if resname in _DONOR_RESNAMES:
            features['donor'] += weight
        if resname in _ACCEPTOR_RESNAMES:
            features['acceptor'] += weight
        if resname in _HYDROPHOBIC_RESNAMES:
            features['hydrophobic'] += weight
        if resname in _AROMATIC_RESNAMES:
            features['aromatic'] += weight
        if resname in _POSITIVE_RESNAMES:
            features['positive'] += weight
        if resname in _NEGATIVE_RESNAMES:
            features['negative'] += weight
    return features


def get_pocket_features_3d(residue_ids: list, resname_map: dict,
                            ca_coord_map: dict) -> list:
    """
    Return a list of (feature_type, xyz_coords) tuples for the pocket, using
    Cα positions as proxies for pharmacophore feature locations.

    A residue may contribute multiple entries (e.g. TYR contributes donor,
    acceptor, hydrophobic, and aromatic). Residues absent from ca_coord_map
    are skipped.

    Args:
        residue_ids: p2rank residue ID strings.
        resname_map: 'CHAIN_RESNUM' → 3-letter residue name.
        ca_coord_map: 'CHAIN_RESNUM' → np.array([x, y, z]) of Cα.

    Returns:
        List of (feat_type: str, coords: np.ndarray) tuples.
    """
    features_3d = []
    for rid in residue_ids:
        resname = _extract_resname(rid, resname_map)
        if not resname:
            continue
        coords = ca_coord_map.get(rid)
        if coords is None:
            continue
        if resname in _DONOR_RESNAMES:
            features_3d.append(('donor', coords))
        if resname in _ACCEPTOR_RESNAMES:
            features_3d.append(('acceptor', coords))
        if resname in _HYDROPHOBIC_RESNAMES:
            features_3d.append(('hydrophobic', coords))
        if resname in _AROMATIC_RESNAMES:
            features_3d.append(('aromatic', coords))
        if resname in _POSITIVE_RESNAMES:
            features_3d.append(('positive', coords))
        if resname in _NEGATIVE_RESNAMES:
            features_3d.append(('negative', coords))
    return features_3d


def get_ligand_features_3d(mol) -> list:
    """
    Return a list of (feature_type, xyz_coords) tuples for the ligand, using
    RDKit pharmacophore feature positions.

    Args:
        mol: RDKit Mol object. Must have a 3D conformer.

    Returns:
        List of (feat_type: str, coords: np.ndarray) tuples.
        Returns empty list if mol is None or has no 3D conformer.
    """
    if mol is None or mol.GetNumConformers() == 0:
        return []
    features_3d = []
    try:
        feats = _get_factory().GetFeaturesForMol(mol)
        for feat in feats:
            family = feat.GetFamily()
            key = _RDKIT_FEATURE_MAP.get(family)
            if key:
                pos = feat.GetPos()
                features_3d.append((key, np.array([pos.x, pos.y, pos.z])))
    except Exception as exc:
        _logger.warning("Failed to extract 3D features from molecule: %s", exc)
    return features_3d


def score_complementarity_3d(pocket_3d: list, ligand_3d: list,
                              sigma: float = 3.0) -> float:
    """
    Rotation-invariant 3D pharmacophore complementarity via pairwise distance matching.

    For each pair of features (i, j) in the pocket, we compute d_pocket = ||p_i - p_j||.
    We then find all ligand feature pairs (k, l) whose types are complementary to (i, j)
    and compute d_ligand = ||l_k - l_l||. The Gaussian kernel

        K(d_p, d_l) = exp(-(d_p - d_l)^2 / (2 * sigma^2))

    rewards similar inter-feature distances; sigma (Å) controls tolerance.
    The best-matching ligand pair is taken for each pocket pair.

    Returns a float in [0, 1] (normalised by evaluated pocket pairs).
    Returns 0.0 if pocket or ligand has fewer than 2 features.
    """
    if len(pocket_3d) < 2 or len(ligand_3d) < 2:
        return 0.0

    # Group ligand features by type for fast lookup
    ligand_by_type: dict = {}
    for feat_type, coords in ligand_3d:
        ligand_by_type.setdefault(feat_type, []).append(coords)

    total_score = 0.0
    n_pairs = 0

    for i in range(len(pocket_3d)):
        for j in range(i + 1, len(pocket_3d)):
            p_type_i, p_coords_i = pocket_3d[i]
            p_type_j, p_coords_j = pocket_3d[j]

            d_pocket = float(np.linalg.norm(p_coords_i - p_coords_j))

            l_type_i = _COMPLEMENTARY[p_type_i]
            l_type_j = _COMPLEMENTARY[p_type_j]

            l_feats_i = ligand_by_type.get(l_type_i, [])
            l_feats_j = ligand_by_type.get(l_type_j, [])

            if not l_feats_i or not l_feats_j:
                continue

            # Best (closest distance match) among all ligand pairs of the correct types
            best_kernel = 0.0
            for l_coords_i in l_feats_i:
                for l_coords_j in l_feats_j:
                    # Skip degenerate case when both types are the same and same point
                    if l_type_i == l_type_j and np.array_equal(l_coords_i, l_coords_j):
                        continue
                    d_ligand = float(np.linalg.norm(l_coords_i - l_coords_j))
                    kernel = float(np.exp(-((d_pocket - d_ligand) ** 2) / (2.0 * sigma ** 2)))
                    if kernel > best_kernel:
                        best_kernel = kernel

            total_score += best_kernel
            n_pairs += 1

    return total_score / n_pairs if n_pairs > 0 else 0.0


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


def score_combined(pocket_features: dict, ligand_features: dict,
                   pocket_3d: list, ligand_3d: list,
                   alpha: float = 0.3, sigma: float = 3.0) -> float:
    """
    Combined pharmacophore score: weighted sum of 1D cosine complementarity
    (Layer 0/1) and rotation-invariant 3D distance matching (Layer 2).

        score = alpha * cosine_score + (1 - alpha) * 3d_score

    Falls back to cosine-only if pocket_3d or ligand_3d is empty.
    """
    cosine = score_complementarity(pocket_features, ligand_features)
    if not pocket_3d or not ligand_3d:
        return cosine
    score_3d = score_complementarity_3d(pocket_3d, ligand_3d, sigma=sigma)
    return alpha * cosine + (1.0 - alpha) * score_3d


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

    # Precompute ligand feature vectors once for all molecules
    active_features = [get_ligand_features(mol) for mol in actives]
    decoy_features = [get_ligand_features(mol) for mol in decoys]

    # Auto-detect p2rank output directory for per-residue probability CSVs (Layer 1)
    p2rank_output_dir = None
    if pdb_dir:
        candidate = os.path.join(os.path.dirname(pdb_dir), 'pockets', 'p2rank_output')
        if os.path.isdir(candidate):
            p2rank_output_dir = candidate
            _logger.info("Layer 1 (probability weighting) enabled: %s", p2rank_output_dir)
        else:
            _logger.info(
                "p2rank output dir not found at %s; Layer 1 probability weighting disabled",
                candidate
            )

    rows = []
    for _, rep in df_reps.iterrows():
        residue_ids = [r for r in str(rep['residues']).split() if r]

        # Build residue name map from representative PDB
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

        # Layer 1: load per-residue p2rank probabilities for probability weighting
        residue_probs = None
        if p2rank_output_dir:
            file_name = str(rep.get('File name', ''))
            # residues CSV name: e.g. 'trajectory_test_data_30.pdb_residues.csv'
            residues_csv_name = file_name.replace('_predictions', '_residues.csv')
            residues_csv_path = os.path.join(p2rank_output_dir, residues_csv_name)
            residue_probs = load_residue_probs(residues_csv_path)
            if not residue_probs:
                _logger.warning("Residue probs not loaded from %s", residues_csv_path)

        # Layer 1: probability-weighted pharmacophore feature vector
        pocket_feats = get_pocket_features(
            residue_ids, resname_map=resname_map, residue_probs=residue_probs
        )

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
