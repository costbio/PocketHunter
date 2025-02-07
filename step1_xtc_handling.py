import os
import subprocess
import glob
import mdtraj as md

def xtc_to_pdb(xtc_file, topology, stride, outfolder, config):
    """Extract frames from XTC file to a series of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting XTC to PDB extraction for {xtc_file}")
        traj = md.load_xtc(xtc_file, topology, stride=stride)
        for i, frame in enumerate(traj):
            real_frame_number = (i+1)*stride
            pdb_file = os.path.join(outfolder, f"{os.path.splitext(os.path.basename(xtc_file))[0]}_{real_frame_number}.pdb")
            frame.save_pdb(pdb_file)
        return [os.path.join(outfolder, f"{os.path.splitext(os.path.basename(xtc_file))[0]}_{(i+1)*stride}.pdb") for i in range(len(traj))]
        logger.info(f"Extraction completed successfully. Saved {len(pdb_files)} PDB files.")
    except Exception as e:
        logger.error(f"Error during XTC to PDB extraction: {str(e)}")
        raise

def write_pdb_list(pdb_dir, output_file, config):
    try:
        """Write a list of PDB file paths relative to pdb_dir to a file."""
        logger = config['logger']
        logger.info(f"Starting PDB list creation for directory {pdb_dir}")
        
        pdb_list = glob.glob(os.path.join(pdb_dir, '*.pdb'))
        if not pdb_list:
            logger.warning(f"No PDB files found in {pdb_dir}")
            return None
            
        with open(output_file, 'w') as f:
            for pdb_file in pdb_list:
                rel_path = os.path.basename(pdb_file)
                f.write(rel_path + '\n')
                logger.debug(f"Added {rel_path} to PDB list")
                
        logger.info(f"PDB list written to {output_file} containing {len(pdb_list)} files")
        return output_file
        
    except Exception as e:
        logger.error(f"Error during PDB list creation: {str(e)}")
        raise

def run_p2rank(pdb_list_file, output_dir, numthreads, config):
    """Run P2Rank with the specified list of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting P2Rank with {numthreads} threads")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        command = f'{dir_path}/tools/p2rank/prank predict {pdb_list_file} -o {output_dir} -threads {numthreads}'
        logger.debug(f"Executing command: {command}")
        
        # Redirect stdout and stderr to DEVNULL to suppress output
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.info(f"P2Rank completed successfully. Output written to {output_dir}")
        return output_dir
        
    except Exception as e:
        logger.error(f"Error during P2Rank execution: {str(e)}")
        raise
