import os
import subprocess
import glob
import mdtraj as md
import tarfile
import shutil

def xtc_to_pdb(xtc_file, topology, stride, outfolder, overwrite, config):
    """Extract frames from XTC file to a series of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting XTC to PDB extraction for {xtc_file}")
        traj = md.load_xtc(xtc_file, topology, stride=stride)

        os.makedirs(outfolder, exist_ok=overwrite)
            
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

def compress_folder(infolder, output_filename, config):
    logger = config['logger']
    logger.info(f"Compressing folder: {infolder} to {output_filename}")
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(infolder, arcname=os.path.basename(infolder))

    # Delete infolder
    logger.info(f"Deleting folder: {infolder}")
    shutil.rmtree(infolder)
    logger.info(f"Compression completed successfully.")

def run_p2rank(pdb_list_file, output_dir, numthreads, novis, compress, overwrite, config):
    """Run P2Rank with the specified list of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting P2Rank with {numthreads} threads")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        p2rank_output_dir = os.path.join(output_dir, 'p2rank_output')
        os.makedirs(p2rank_output_dir, exist_ok=overwrite)

        command = f'{dir_path}/tools/p2rank/prank predict {pdb_list_file} -o {p2rank_output_dir} -threads {numthreads}'
        
        if novis == True:
            command += ' -visualizations 0'
        logger.info(f"Executing command: {command}")
        
        # Redirect stdout and stderr to DEVNULL to suppress output
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.info(f"P2Rank completed successfully. Output written to {output_dir}")

        if compress:
            compress_folder(p2rank_output_dir, f"{p2rank_output_dir}.tar.gz", config)
            return f"{p2rank_output_dir}.tar.gz"
        else:
            return output_dir
        
    except Exception as e:
        logger.error(f"Error during P2Rank execution: {str(e)}")
        raise
