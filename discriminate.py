"""
discriminate.py — Pharmacophore-based active/decoy discrimination.

Given a pocket_clusters directory (output of cluster_pockets) and two SDF files
(actives and decoys), scores each cluster representative by its ability to
discriminate actives from decoys using residue-based pharmacophore vectors and
complementarity scoring.

Scoring layers
--------------
Layer 0 (baseline):  Binary residue counting + cosine complementarity.
Layer 1 (+prob):     p2rank per-residue binding probabilities replace binary
                     weights — core residues contribute more than edge residues.
Layer 2 (+3D):       Rotation-invariant 3D pharmacophore complementarity via
                     Gaussian feature overlap. Sidechain functional group atom
                     coordinates serve as pharmacophore feature positions.
                     Combined score = alpha × cosine + (1−alpha) × 3D.

                     The Gaussian kernel formula derives from the RDKit FeatMaps
                     implementation (FeatMapParams.FeatProfile.Gaussian), which
                     itself is described in Landrum, Penzotti & Putta (2006).

Both layers activate automatically when the required data are present (p2rank
residues CSV for Layer 1; representative PDB for Layer 2).

References
----------
Gaussian kernel (feature-map scoring):
    Landrum, Penzotti & Putta, J. Comput. Aided Mol. Des. 20, 751–762 (2006).
    DOI: https://doi.org/10.1007/s10822-006-9085-8

Distance-pattern pharmacophore complementarity (Layer 2 methodology):
    Mahé, Ralaivola, Stoven & Vert, J. Chem. Inf. Model. 46, 2003–2014 (2006).
    DOI: https://doi.org/10.1021/ci060138m

Structure-based pharmacophore feature definitions:
    Wolber & Langer, J. Chem. Inf. Model. 45, 160–169 (2005).
    DOI: https://doi.org/10.1021/ci049885e

    IUPAC, Pure Appl. Chem. 70, 1129–1143 (1998).
    DOI: https://doi.org/10.1351/pac199870051129

Pharmacophore feature types (RDKit):
    RDKit BaseFeatures.fdef — Donor, Acceptor, Hydrophobe, Aromatic,
    PosIonizable, NegIonizable.

Pocket prediction:
    Krivák & Hoksza, J. Cheminform. 10, 39 (2018).
    DOI: https://doi.org/10.1186/s13321-018-0285-8

Evaluation metrics:
    Truchon & Bayly, J. Chem. Inf. Model. 47, 488–508 (2007).
    DOI: https://doi.org/10.1021/ci600426e
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

# ── Residue → pharmacophore feature mapping ─────────────────────────────────
# Reference: Wolber & Langer, JCIM 2005 (LigandScout feature definitions)
#           IUPAC pharmacophore definition, Pure Appl. Chem. 1998
#           RDKit BaseFeatures.fdef feature families

_DONOR_RESNAMES = {'SER', 'THR', 'TYR', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'TRP'}
_ACCEPTOR_RESNAMES = {'ASP', 'GLU', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'HIS'}
_HYDROPHOBIC_RESNAMES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}
_AROMATIC_RESNAMES = {'PHE', 'TRP', 'TYR', 'HIS'}
_POSITIVE_RESNAMES = {'LYS', 'ARG', 'HIS'}
_NEGATIVE_RESNAMES = {'ASP', 'GLU'}

# ── Sidechain atom name → pharmacophore feature type mapping ─────────────────
# Atom names follow PDB convention (ATOM records, HETATM excluded).
# Each list specifies the sidechain atoms whose centroid defines the
# pharmacophore feature position for that residue.
#
# Reference for atom naming: PDB Chemical Component Dictionary (wwPDB)
# Reference for which atoms define which features:
#   - H-bond donor/acceptor: Foldit Sidechain Bonding Table + Wolber 2005
#   - Hydrophobic carbons: sidechain C atoms excluding polar functional groups
#   - Aromatic: ring centroids per Wolber 2005 LigandScout definitions
#   - Charge: ionizable group centroids per IUPAC definitions

_SIDECHAIN_FEATURE_ATOMS = {
    # Hydroxyl sidechains — donor and acceptor share the same atom
    'SER':  {'donor': ['OG'],
             'acceptor': ['OG']},
    'THR':  {'donor': ['OG1'],
             'acceptor': ['OG1']},
    # Phenol — donor, acceptor, aromatic, hydrophobic
    'TYR':  {'donor': ['OH'],
             'acceptor': ['OH'],
             'aromatic': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
             'hydrophobic': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']},
    # Amide sidechains — each oxygen is acceptor, each nitrogen is donor
    'ASN':  {'donor': ['ND2'],
             'acceptor': ['OD1']},
    'GLN':  {'donor': ['NE2'],
             'acceptor': ['OE1']},
    # Positively ionizable donors
    'LYS':  {'donor': ['NZ'],
             'positive': ['NZ']},
    'ARG':  {'donor': ['NH1', 'NH2'],
             'positive': ['CZ', 'NH1', 'NH2']},
    # Imidazole — donor, acceptor, aromatic, positive (protonated at physiological pH)
    'HIS':  {'donor': ['ND1', 'NE2'],
             'acceptor': ['ND1', 'NE2'],
             'aromatic': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
             'positive': ['ND1', 'NE2']},
    # Indole — donor, aromatic, hydrophobic
    'TRP':  {'donor': ['NE1'],
             'aromatic': ['CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
             'hydrophobic': ['CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']},
    # Carboxylate sidechains — acceptor and negatively ionizable
    'ASP':  {'acceptor': ['OD1', 'OD2'],
             'negative': ['CG', 'OD1', 'OD2']},
    'GLU':  {'acceptor': ['OE1', 'OE2'],
             'negative': ['CD', 'OE1', 'OE2']},
    # Aliphatic hydrophobic
    'ALA':  {'hydrophobic': ['CB']},
    'VAL':  {'hydrophobic': ['CG1', 'CG2']},
    'ILE':  {'hydrophobic': ['CG2', 'CD1']},
    'LEU':  {'hydrophobic': ['CD1', 'CD2']},
    'MET':  {'hydrophobic': ['SD', 'CE']},
    'PRO':  {'hydrophobic': ['CG', 'CD']},
    # Aromatic + hydrophobic
    'PHE':  {'aromatic': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
             'hydrophobic': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']},
}

_RDKIT_FEATURE_MAP = {
    'Donor': 'donor',
    'Acceptor': 'acceptor',
    'Hydrophobe': 'hydrophobic',
    'LumpedHydrophobe': 'hydrophobic',
    'Aromatic': 'aromatic',
    'PosIonizable': 'positive',
    'NegIonizable': 'negative',
}

# Complementary feature mapping for protein-ligand interaction:
# pocket donor ↔ ligand acceptor, pocket positive ↔ ligand negative, etc.
# Self-complementary: hydrophobic, aromatic (pi-stacking).
# Reference: IUPAC pharmacophore definition + Wolber 2005
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

# Gaussian kernel width (sigma, in Å) for the distance-pattern scoring kernel.
# 2.5 Å gives a comfortable tolerance for intramolecular distance comparison;
# tighter values (< 1.5 Å) yield near-zero scores for most feature pairs.
# Overridden by the width parameter in run_discrimination (default 2.5 Å).
_GAUSSIAN_WIDTH = 2.5

_logger = logging.getLogger(__name__)


def _get_factory():
    """Return the singleton RDKit ChemicalFeatures factory (BaseFeatures.fdef)."""
    global _factory_instance
    if _factory_instance is None:
        if not os.path.exists(_FDEF_PATH):
            raise RuntimeError(
                f"RDKit feature definition file not found: {_FDEF_PATH}. "
                "Ensure rdkit data files are installed (e.g. conda install -c conda-forge rdkit)."
            )
        _factory_instance = ChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)
    return _factory_instance


# ═══════════════════════════════════════════════════════════════════════════════
# Residue identification
# ═══════════════════════════════════════════════════════════════════════════════

def _extract_resname(residue_id: str, resname_map: dict = None) -> str:
    """
    Extract the 3-letter residue name from a p2rank residue ID.

    Tries lookup in resname_map first (for 'CHAIN_RESNUM' format like 'A_1019'),
    then falls back to parsing the ID string (for 'CHAIN_NUM_RESNAME' format
    like 'A_1_ALA').
    """
    rid = residue_id.strip()

    if resname_map and rid in resname_map:
        return resname_map[rid]

    parts = rid.split('_')
    for part in reversed(parts):
        clean = part.strip().upper()
        if len(clean) == 3 and clean.isalpha():
            return clean

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

    Uses prody to read the structure (Cα subset, sufficient for residue identity).
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


# ═══════════════════════════════════════════════════════════════════════════════
# Sidechain coordinate extraction
# ═══════════════════════════════════════════════════════════════════════════════

def build_sidechain_coords_map(pdb_path: str) -> dict:
    """
    Build mapping from 'CHAIN_RESNUM' to sidechain pharmacophore feature
    coordinates by parsing a PDB file (all atoms).

    For each residue, extracts the 3D centroid of the sidechain atoms
    that define each pharmacophore feature type (donor, acceptor, hydrophobic,
    aromatic, positive, negative).

    Atom selection rules follow Wolber & Langer (2005) LigandScout feature
    definitions and the IUPAC pharmacophore definition (1998).
    Specific PDB atom names are derived from the wwPDB Chemical Component
    Dictionary.

    Returns:
        dict like {
            'A_1019': {
                'donor': np.array([x, y, z]),
                'acceptor': np.array([x, y, z]),
                'hydrophobic': np.array([x, y, z]),
                'aromatic': np.array([x, y, z]),
                'positive': np.array([x, y, z]),
                'negative': np.array([x, y, z]),
            },
            ...
        }
        Residues with no feature atoms map to empty dicts.
        Returns empty dict on parse failure.
    """
    try:
        structure = parsePDB(pdb_path, subset=None, verbosity='none')
        if structure is None:
            _logger.warning("prody could not parse PDB file: %s", pdb_path)
            return {}
    except Exception as exc:
        _logger.warning("Failed to parse %s for sidechain coords: %s", pdb_path, exc)
        return {}

    # Group atoms by residue, indexed by 'CHAIN_RESNUM'
    residue_atoms: dict = {}  # key → {atom_name: np.array([x,y,z])}
    residue_resnames: dict = {}  # key → 3-letter resname
    try:
        for atom in structure.iterAtoms():
            chain = atom.getChid()
            resnum = atom.getResnum()
            key = f"{chain}_{resnum}"
            name = atom.getName()
            coords = atom.getCoords()
            if coords is None:
                continue
            residue_atoms.setdefault(key, {})[name] = coords.copy()
            if key not in residue_resnames:
                residue_resnames[key] = atom.getResname()
    except Exception as exc:
        _logger.warning("Failed to iterate atoms in %s: %s", pdb_path, exc)
        return {}

    # For each residue, compute feature centroids from sidechain atoms
    coords_map: dict = {}
    for key, atoms in residue_atoms.items():
        resname = residue_resnames.get(key, '')
        feature_atom_map = _SIDECHAIN_FEATURE_ATOMS.get(resname, {})
        if not feature_atom_map:
            continue

        feat_coords: dict = {}
        for feat_type, atom_names in feature_atom_map.items():
            positions = []
            for aname in atom_names:
                if aname in atoms:
                    positions.append(atoms[aname])
            if positions:
                feat_coords[feat_type] = np.mean(positions, axis=0)

        if feat_coords:
            coords_map[key] = feat_coords

    return coords_map


# ═══════════════════════════════════════════════════════════════════════════════
# p2rank probability loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_residue_probs(residues_csv_path: str) -> dict:
    """
    Load per-residue binding probabilities from a p2rank residues CSV file.

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


# ═══════════════════════════════════════════════════════════════════════════════
# 1D pharmacophore feature counting (Layer 0 / Layer 1)
# ═══════════════════════════════════════════════════════════════════════════════

def get_pocket_features(residue_ids: list, resname_map: dict = None,
                        residue_probs: dict = None) -> dict:
    """
    Compute pharmacophore feature counts (or probability-weighted sums) for a pocket.

    Args:
        residue_ids: List of p2rank residue ID strings.
        resname_map: Optional dict mapping 'CHAIN_RESNUM' to 3-letter residue name.
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


# ═══════════════════════════════════════════════════════════════════════════════
# 3D pharmacophore feature extraction (Layer 2)
# ═══════════════════════════════════════════════════════════════════════════════

def get_pocket_features_3d(residue_ids: list, resname_map: dict,
                            sidechain_coords_map: dict) -> list:
    """
    Return a list of (feature_type, xyz_coords) tuples for the pocket,
    using sidechain functional group atom centroids as pharmacophore
    feature positions.

    Feature positions are sourced from the specific sidechain atoms that
    mediate each interaction type (e.g. OD1/OD2 for ASP acceptor, NZ for
    LYS donor). See _SIDECHAIN_FEATURE_ATOMS for the complete mapping.

    Reference: Wolber & Langer, JCIM 2005 (LigandScout feature positions)
               IUPAC pharmacophore definition, Pure Appl. Chem. 1998

    Args:
        residue_ids: p2rank residue ID strings.
        resname_map: 'CHAIN_RESNUM' → 3-letter residue name.
        sidechain_coords_map: 'CHAIN_RESNUM' → {feat_type: np.array([x,y,z]), ...}
                              from build_sidechain_coords_map().

    Returns:
        List of (feat_type: str, coords: np.ndarray) tuples.
    """
    features_3d = []
    for rid in residue_ids:
        resname = _extract_resname(rid, resname_map)
        if not resname or rid not in sidechain_coords_map:
            continue
        res_feats = sidechain_coords_map[rid]
        for feat_type in FEATURE_KEYS:
            if feat_type in res_feats:
                features_3d.append((feat_type, res_feats[feat_type]))
    return features_3d


def get_ligand_features_3d(mol) -> list:
    """
    Return a list of (feature_type, xyz_coords) tuples for the ligand,
    using RDKit pharmacophore feature positions.

    Feature positions are RDKit-calculated centroids of the atoms matching
    each pharmacophore feature definition in BaseFeatures.fdef.

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


# ═══════════════════════════════════════════════════════════════════════════════
# 3D Gaussian feature overlap scoring
# ═══════════════════════════════════════════════════════════════════════════════

def _gaussian_featfeat_score(distance: float, width: float = None) -> float:
    """
    Gaussian kernel for feature-feature overlap.

    Formula: score = exp(-d² / (2 * width²))

    This is the mathematical formula underlying RDKit's
    FeatMapParams.FeatProfile.Gaussian, as referenced in:
        Landrum, Penzotti & Putta, J. Comput. Aided Mol. Des. 20, 751–762 (2006).

    Args:
        distance: Euclidean distance between two feature positions (Å).
        width: Gaussian width parameter (sigma, in Å). Defaults to _GAUSSIAN_WIDTH.

    Returns:
        Float in [0, 1]. 1 = identical positions, →0 at large distances.
    """
    if width is None:
        width = _GAUSSIAN_WIDTH
    return float(np.exp(-(distance ** 2) / (2.0 * width ** 2)))


def score_complementarity_3d(pocket_3d: list, ligand_3d: list,
                              width: float = None) -> float:
    """
    Rotation-invariant 3D pharmacophore complementarity via pairwise distance
    matching with Gaussian feature overlap kernel.

    Approach:
        For each pair of features (i, j) in the pocket, compute the Euclidean
        distance d_pocket = ||pos_i - pos_j||. Find all ligand feature pairs
        (k, l) whose pharmacophore types are complementary to (type_i, type_j)
        and compute d_ligand = ||pos_k - pos_l||.

        The similarity of the two distance patterns is scored with a Gaussian
        kernel:
            K(d_p, d_l) = exp(-(d_p - d_l)² / (2 * width²))

        For each pocket pair, the best-matching ligand pair kernel value is
        taken. The final score is the mean over all pocket pairs.

    The distance-pattern comparison approach (scoring |d_pocket − d_ligand|
    rather than positional overlap) achieves rotation invariance without
    superposition, following the pharmacophore kernel framework of Mahé,
    Ralaivola, Stoven & Vert (J. Chem. Inf. Model. 46, 2003–2014, 2006;
    DOI: 10.1021/ci060138m). The Gaussian kernel function exp(−d²/2σ²) is
    from Landrum, Penzotti & Putta (J. Comput. Aided Mol. Des. 20, 751–762,
    2006; DOI: 10.1007/s10822-006-9085-8).

    Args:
        pocket_3d: List of (feat_type, coords) tuples for the pocket.
        ligand_3d: List of (feat_type, coords) tuples for the ligand.
        width: Gaussian width parameter (Å). Defaults to _GAUSSIAN_WIDTH.

    Returns:
        Float in [0, 1]. Returns 0.0 if fewer than 2 features on either side.
    """
    if len(pocket_3d) < 2 or len(ligand_3d) < 2:
        return 0.0

    # Group ligand features by complementary type for lookup
    ligand_by_type: dict = {}
    for feat_type, coords in ligand_3d:
        comp_type = _COMPLEMENTARY[feat_type]
        ligand_by_type.setdefault(comp_type, []).append(coords)

    total_score = 0.0
    n_pairs = 0

    for i in range(len(pocket_3d)):
        for j in range(i + 1, len(pocket_3d)):
            p_type_i, p_coords_i = pocket_3d[i]
            p_type_j, p_coords_j = pocket_3d[j]

            d_pocket = float(np.linalg.norm(p_coords_i - p_coords_j))

            # Find complementary ligand feature types
            l_feats_i = ligand_by_type.get(p_type_i, [])
            l_feats_j = ligand_by_type.get(p_type_j, [])

            if not l_feats_i or not l_feats_j:
                continue

            # Best (closest distance match) among all ligand pairs
            best_kernel = 0.0
            for l_coords_i in l_feats_i:
                for l_coords_j in l_feats_j:
                    if (p_type_i == p_type_j
                            and np.array_equal(l_coords_i, l_coords_j)):
                        continue
                    d_ligand = float(np.linalg.norm(l_coords_i - l_coords_j))
                    kernel = _gaussian_featfeat_score(
                        abs(d_pocket - d_ligand), width=width
                    )
                    if kernel > best_kernel:
                        best_kernel = kernel

            total_score += best_kernel
            n_pairs += 1

    return total_score / n_pairs if n_pairs > 0 else 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# 1D cosine complementarity
# ═══════════════════════════════════════════════════════════════════════════════

def get_ligand_features(mol) -> dict:
    """Extract pharmacophore feature counts from a ligand molecule using RDKit."""
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
    Compute pharmacophoric complementarity between a pocket and a ligand
    via cosine similarity with complementary feature mapping.

    H-bond and charge features are cross-mapped:
      pocket donor    ↔ ligand acceptor
      pocket acceptor ↔ ligand donor
      pocket positive ↔ ligand negative
      pocket negative ↔ ligand positive

    Hydrophobic and aromatic features are self-complementary.

    Returns float in [0, 1]. Returns 0.0 if either vector is all zeros.
    """
    # Remap pocket features to their complementary counterparts
    complementary_pocket = {
        'donor':       pocket_features['acceptor'],
        'acceptor':    pocket_features['donor'],
        'hydrophobic': pocket_features['hydrophobic'],
        'aromatic':    pocket_features['aromatic'],
        'positive':    pocket_features['negative'],
        'negative':    pocket_features['positive'],
    }
    p = np.array([complementary_pocket[k] for k in FEATURE_KEYS], dtype=float)
    l = np.array([ligand_features[k] for k in FEATURE_KEYS], dtype=float)
    norm_p = np.linalg.norm(p)
    norm_l = np.linalg.norm(l)
    if norm_p == 0 or norm_l == 0:
        return 0.0
    return float(np.dot(p, l) / (norm_p * norm_l))


# ═══════════════════════════════════════════════════════════════════════════════
# Combined scoring
# ═══════════════════════════════════════════════════════════════════════════════

def score_combined(pocket_features: dict, ligand_features: dict,
                   pocket_3d: list, ligand_3d: list,
                   alpha: float = 0.9, width: float = None) -> float:
    """
    Combined pharmacophore score: weighted sum of 1D cosine complementarity
    and rotation-invariant 3D Gaussian feature overlap.

        score = alpha × cosine + (1 − alpha) × 3D

    Falls back to cosine-only if pocket_3d or ligand_3d is empty.

    Args:
        pocket_features: Dict of 1D pharmacophore feature counts/weights.
        ligand_features: Dict of 1D pharmacophore feature counts.
        pocket_3d: List of (feat_type, coords) for pocket sidechain features.
        ligand_3d: List of (feat_type, coords) for ligand RDKit features.
        alpha: Weight of 1D cosine score. 1−alpha is weight of 3D score.
               Default 0.9 weights cosine heavily; use 0.5 for equal weighting.
        width: Gaussian width for 3D scoring (Å). Defaults to _GAUSSIAN_WIDTH.

    Returns:
        Float in [0, 1].
    """
    cosine = score_complementarity(pocket_features, ligand_features)
    if not pocket_3d or not ligand_3d:
        return cosine
    score_3d = score_complementarity_3d(pocket_3d, ligand_3d, width=width)
    return alpha * cosine + (1.0 - alpha) * score_3d


# ═══════════════════════════════════════════════════════════════════════════════
# Evaluation metrics
# ═══════════════════════════════════════════════════════════════════════════════

def _enrichment_factor(scores_actives: list, scores_decoys: list,
                       fraction: float) -> float:
    """
    Compute enrichment factor at a given top fraction.

    EF_x% = (actives in top x%) / (x% × total_actives)

    Reference: Truchon & Bayly, J. Chem. Inf. Model. 47, 488–508 (2007).
    """
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
    """
    Compute discrimination metrics: ROC-AUC, EF1%, EF5%.

    Reference: Truchon & Bayly, J. Chem. Inf. Model. 47, 488–508 (2007).
    """
    labels = [1] * len(scores_actives) + [0] * len(scores_decoys)
    scores = list(scores_actives) + list(scores_decoys)
    try:
        roc_auc = float(roc_auc_score(labels, scores))
    except ValueError:
        roc_auc = 0.5
    ef1 = _enrichment_factor(scores_actives, scores_decoys, 0.01)
    ef5 = _enrichment_factor(scores_actives, scores_decoys, 0.05)
    return {'roc_auc': roc_auc, 'ef1': ef1, 'ef5': ef5}


# ═══════════════════════════════════════════════════════════════════════════════
# Main pipeline
# ═══════════════════════════════════════════════════════════════════════════════

def _load_mols_from_sdf(sdf_path: str) -> list:
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    return [mol for mol in supplier if mol is not None]


def run_discrimination(cluster_dir: str, actives_sdf: str, decoys_sdf: str,
                       outfolder: str, pdb_dir: str = None,
                       alpha: float = 0.9, width: float = 2.5) -> pd.DataFrame:
    """
    Run pharmacophore-based active/decoy discrimination for all cluster
    representatives.

    Scoring pipeline:
      1. Load cluster representatives from cluster_representatives.csv.
      2. Load active and decoy molecules from SDF files.
      3. For each representative:
         a. Identify pocket residue IDs.
         b. Compute 1D pharmacophore feature vector (Layer 0/1).
         c. If PDB available, extract sidechain 3D feature coordinates (Layer 2).
         d. Score all actives and decoys against the pocket pharmacophore.
         e. Compute discrimination metrics (ROC-AUC, EF1%, EF5%).

    Args:
        cluster_dir: Path to pocket_clusters output directory.
        actives_sdf: Path to actives SDF file (max 2000 molecules).
        decoys_sdf: Path to decoys SDF file (max 2000 molecules).
        outfolder: Output directory for discrimination_results.csv.
        pdb_dir: Directory containing representative PDB files.
        alpha: Weight of 1D cosine in combined score (default 0.9).
        width: Gaussian width for 3D scoring (default 2.5 Å).

    Returns:
        DataFrame sorted by ROC-AUC (descending).

    References:
        - Feature extraction: RDKit BaseFeatures.fdef
        - Sidechain feature positions: Wolber & Langer (2005), IUPAC (1998)
        - Gaussian kernel: Landrum, Penzotti & Putta (2006)
        - Evaluation metrics: Truchon & Bayly (2007)
    """
    if width is None:
        width = _GAUSSIAN_WIDTH

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

    if len(actives) > 2000:
        raise ValueError(f"Too many actives: {len(actives)} (max 2000)")
    if len(decoys) > 5000:
        raise ValueError(f"Too many decoys: {len(decoys)} (max 5000)")

    # Precompute 1D and 3D ligand features once
    active_features_1d = [get_ligand_features(mol) for mol in actives]
    decoy_features_1d = [get_ligand_features(mol) for mol in decoys]
    active_features_3d = [get_ligand_features_3d(mol) for mol in actives]
    decoy_features_3d = [get_ligand_features_3d(mol) for mol in decoys]

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

        # Build residue name map and sidechain coordinate map from PDB
        resname_map = {}
        sidechain_coords_map = {}
        if pdb_dir:
            file_name = str(rep.get('File name', ''))
            pdb_name = file_name.replace('_predictions', '')
            pdb_path = os.path.join(pdb_dir, pdb_name)
            if os.path.exists(pdb_path):
                resname_map = build_resname_map(pdb_path)
                sidechain_coords_map = build_sidechain_coords_map(pdb_path)
            else:
                _logger.warning("Representative PDB not found: %s", pdb_path)

        # Layer 1: per-residue p2rank probability weighting
        residue_probs = None
        if p2rank_output_dir:
            file_name = str(rep.get('File name', ''))
            residues_csv_name = file_name.replace('_predictions', '_residues.csv')
            residues_csv_path = os.path.join(p2rank_output_dir, residues_csv_name)
            residue_probs = load_residue_probs(residues_csv_path)
            if not residue_probs:
                _logger.warning("Residue probs not loaded from %s", residues_csv_path)

        # Layer 1: probability-weighted pharmacophore feature vector
        pocket_feats = get_pocket_features(
            residue_ids, resname_map=resname_map, residue_probs=residue_probs
        )

        # Layer 2: 3D feature positions from sidechain functional group atoms
        pocket_3d = []
        if sidechain_coords_map:
            pocket_3d = get_pocket_features_3d(
                residue_ids, resname_map, sidechain_coords_map
            )

        scores_actives = [
            score_combined(
                pocket_feats, lf_1d, pocket_3d, lf_3d, alpha=alpha, width=width
            )
            for lf_1d, lf_3d in zip(active_features_1d, active_features_3d)
        ]
        scores_decoys = [
            score_combined(
                pocket_feats, df_1d, pocket_3d, df_3d, alpha=alpha, width=width
            )
            for df_1d, df_3d in zip(decoy_features_1d, decoy_features_3d)
        ]
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

    df_results = pd.DataFrame(rows).sort_values(
        'roc_auc', ascending=False
    ).reset_index(drop=True)
    df_results.to_csv(os.path.join(outfolder, 'discrimination_results.csv'), index=False)
    return df_results
