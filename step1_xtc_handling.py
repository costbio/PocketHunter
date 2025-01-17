import os
import subprocess
import glob
import mdtraj as md


def xtc_to_pdb(xtc_file, topology, pdb_dir):
    """Convert XTC file to a series of PDB files."""
    traj = md.load_xtc(xtc_file, topology)
    for i, frame in enumerate(traj):
        pdb_file = os.path.join(pdb_dir, f"{os.path.splitext(os.path.basename(xtc_file))[0]}_{i}.pdb")
        frame.save_pdb(pdb_file)
        print(f"Frame {i} saved as {pdb_file}")
    return [os.path.join(pdb_dir, f"{os.path.splitext(os.path.basename(xtc_file))[0]}_{i}.pdb") for i in range(len(traj))]

def write_pdb_list(pdb_dir, output_file):
    """Write a list of PDB file paths relative to pdb_dir to a file."""
    pdb_list = glob.glob(os.path.join(pdb_dir, '*.pdb'))
    with open(output_file, 'w') as f:
        for pdb_file in pdb_list:
            rel_path = os.path.relpath(pdb_file, start=pdb_dir)
            f.write(rel_path + '\n')
    print(f"PDB list written to {output_file}")
    return output_file

def run_p2rank(pdb_list_file, output_dir, threads=4):
    """Run P2Rank with the specified list of PDB files."""
    command = f'prank predict {pdb_list_file} -o {output_dir} -threads {threads}'
    subprocess.call(command, shell=True)
    print(f"P2Rank output written to {output_dir}")
    return output_dir
