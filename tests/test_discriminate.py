import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from discriminate import (
    get_pocket_features,
    get_ligand_features,
    score_complementarity,
    compute_discrimination_metrics,
    run_discrimination,
    FEATURE_KEYS,
)


def test_get_pocket_features_returns_all_keys():
    features = get_pocket_features(['A_1_ALA'])
    assert set(features.keys()) == set(FEATURE_KEYS)

def test_get_pocket_features_hydrophobic():
    features = get_pocket_features(['A_1_ALA', 'A_2_VAL'])
    assert features['hydrophobic'] == 2

def test_get_pocket_features_donor():
    features = get_pocket_features(['A_1_ARG'])
    assert features['donor'] >= 1

def test_get_pocket_features_acceptor():
    features = get_pocket_features(['A_1_ASP'])
    assert features['acceptor'] >= 1

def test_get_pocket_features_aromatic():
    features = get_pocket_features(['A_1_PHE'])
    assert features['aromatic'] >= 1

def test_get_pocket_features_positive():
    features = get_pocket_features(['A_1_LYS'])
    assert features['positive'] >= 1

def test_get_pocket_features_negative():
    features = get_pocket_features(['A_1_GLU'])
    assert features['negative'] >= 1

def test_get_pocket_features_empty():
    features = get_pocket_features([])
    assert all(v == 0 for v in features.values())

def test_get_pocket_features_unknown_residue():
    features = get_pocket_features(['A_1_UNK'])
    assert all(v >= 0 for v in features.values())

def test_get_ligand_features_returns_all_keys():
    from rdkit import Chem
    mol = Chem.MolFromSmiles('c1ccccc1')
    mol = Chem.AddHs(mol)
    features = get_ligand_features(mol)
    assert set(features.keys()) == set(FEATURE_KEYS)

def test_get_ligand_features_aromatic_benzene():
    from rdkit import Chem
    mol = Chem.MolFromSmiles('c1ccccc1')
    mol = Chem.AddHs(mol)
    features = get_ligand_features(mol)
    assert features['aromatic'] >= 1

def test_get_ligand_features_donor_ethanol():
    from rdkit import Chem
    mol = Chem.MolFromSmiles('CCO')
    mol = Chem.AddHs(mol)
    features = get_ligand_features(mol)
    assert features['donor'] >= 1

def test_get_ligand_features_none_mol():
    features = get_ligand_features(None)
    assert all(v == 0 for v in features.values())

def test_score_complementarity_identical_vectors():
    pocket = {'donor': 2, 'acceptor': 0, 'hydrophobic': 0, 'aromatic': 0, 'positive': 0, 'negative': 0}
    ligand = {'donor': 2, 'acceptor': 0, 'hydrophobic': 0, 'aromatic': 0, 'positive': 0, 'negative': 0}
    score = score_complementarity(pocket, ligand)
    assert score == pytest.approx(1.0)

def test_score_complementarity_orthogonal_vectors():
    pocket = {'donor': 1, 'acceptor': 0, 'hydrophobic': 0, 'aromatic': 0, 'positive': 0, 'negative': 0}
    ligand = {'donor': 0, 'acceptor': 1, 'hydrophobic': 0, 'aromatic': 0, 'positive': 0, 'negative': 0}
    score = score_complementarity(pocket, ligand)
    assert score == pytest.approx(0.0)

def test_score_complementarity_zero_pocket():
    pocket = {k: 0 for k in FEATURE_KEYS}
    ligand = {'donor': 1, 'acceptor': 0, 'hydrophobic': 0, 'aromatic': 0, 'positive': 0, 'negative': 0}
    score = score_complementarity(pocket, ligand)
    assert score == pytest.approx(0.0)

def test_score_complementarity_range():
    pocket = {'donor': 2, 'acceptor': 1, 'hydrophobic': 3, 'aromatic': 1, 'positive': 0, 'negative': 0}
    ligand = {'donor': 1, 'acceptor': 2, 'hydrophobic': 1, 'aromatic': 0, 'positive': 1, 'negative': 0}
    score = score_complementarity(pocket, ligand)
    assert 0.0 <= score <= 1.0

def test_compute_metrics_perfect_discrimination():
    scores_actives = [1.0, 0.9, 0.8]
    scores_decoys = [0.3, 0.2, 0.1]
    metrics = compute_discrimination_metrics(scores_actives, scores_decoys)
    assert metrics['roc_auc'] == pytest.approx(1.0)

def test_compute_metrics_random_discrimination():
    import random
    random.seed(42)
    scores_actives = [random.random() for _ in range(100)]
    scores_decoys = [random.random() for _ in range(100)]
    metrics = compute_discrimination_metrics(scores_actives, scores_decoys)
    assert 0.0 <= metrics['roc_auc'] <= 1.0

def test_compute_metrics_keys_present():
    metrics = compute_discrimination_metrics([0.9, 0.8], [0.2, 0.1])
    assert 'roc_auc' in metrics
    assert 'ef1' in metrics
    assert 'ef5' in metrics

def test_compute_metrics_ef1_perfect():
    scores_actives = [1.0] * 10
    scores_decoys = [0.0] * 90
    metrics = compute_discrimination_metrics(scores_actives, scores_decoys)
    assert metrics['ef1'] > 1.0

def test_run_discrimination_missing_columns(tmp_path):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Write a valid minimal SDF so SDF loading doesn't fail first
    mol = Chem.MolFromSmiles('CCO')
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    writer = Chem.SDWriter(str(tmp_path / 'actives.sdf'))
    writer.write(mol)
    writer.close()
    writer2 = Chem.SDWriter(str(tmp_path / 'decoys.sdf'))
    writer2.write(mol)
    writer2.close()

    # CSV missing 'cluster' column — should raise before scoring
    bad_csv = tmp_path / 'cluster_representatives.csv'
    pd.DataFrame({'residues': ['A_1_ALA A_2_ARG'], 'Frame': [1]}).to_csv(bad_csv, index=False)

    with pytest.raises(ValueError, match='missing required columns'):
        from discriminate import run_discrimination
        run_discrimination(str(tmp_path), str(tmp_path / 'actives.sdf'),
                           str(tmp_path / 'decoys.sdf'), str(tmp_path / 'out'))


def test_run_discrimination_happy_path(tmp_path):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def _write_sdf(path, smiles_list):
        writer = Chem.SDWriter(str(path))
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            writer.write(mol)
        writer.close()

    _write_sdf(tmp_path / 'actives.sdf', ['c1ccccc1O', 'CCN'])   # 2 actives
    _write_sdf(tmp_path / 'decoys.sdf', ['CCCC', 'CCC', 'CC'])   # 3 decoys

    reps_csv = tmp_path / 'cluster_representatives.csv'
    pd.DataFrame({
        'cluster': [0, 1],
        'Frame': [10, 20],
        'residues': ['A_1_ARG A_2_ASP A_3_PHE', 'A_4_ALA A_5_VAL A_6_LEU'],
        'probability': [0.85, 0.72],
    }).to_csv(reps_csv, index=False)

    from discriminate import run_discrimination
    df = run_discrimination(
        cluster_dir=str(tmp_path),
        actives_sdf=str(tmp_path / 'actives.sdf'),
        decoys_sdf=str(tmp_path / 'decoys.sdf'),
        outfolder=str(tmp_path / 'out'),
    )

    assert len(df) == 2
    assert list(df.columns) >= ['cluster_id', 'frame', 'roc_auc', 'ef1', 'ef5']
    assert df['roc_auc'].between(0.0, 1.0).all()
    assert (tmp_path / 'out' / 'discrimination_results.csv').exists()
    # Results sorted by roc_auc descending
    assert df['roc_auc'].iloc[0] >= df['roc_auc'].iloc[1]
