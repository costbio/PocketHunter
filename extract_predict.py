import os
import subprocess
import glob
import mdtraj as md
import tarfile
import shutil
import process_predicted

def xtc_to_pdb(xtc_file, topology, stride, outfolder, overwrite, config):
    """Extract frames from XTC file to a series of PDB files."""
    try:
        logger = config['logger']
        logger.info(f"Starting XTC to PDB extraction for {xtc_file}")
        traj = md.load_xtc(xtc_file, topology, stride=stride)
        
        total_frames = len(traj)
        logger.info(f"Total frames to extract: {total_frames}")
        
        # Calculate progress reporting intervals (every 10%)
        progress_interval = max(1, total_frames // 10)
        
        # Ensure output directory exists
        os.makedirs(outfolder, exist_ok=True)
        pdb_files = []
        for i, frame in enumerate(traj):
            real_frame_number = (i+1)*stride
            pdb_file = os.path.join(outfolder, f"{os.path.splitext(os.path.basename(xtc_file))[0]}_{real_frame_number}.pdb")
            frame.save_pdb(pdb_file)
            pdb_files.append(pdb_file)
            # Log progress every 10%
            if (i + 1) % progress_interval == 0 or (i + 1) == total_frames:
                progress_percent = ((i + 1) / total_frames) * 100
                logger.info(f"Extraction progress: {i + 1}/{total_frames} frames ({progress_percent:.1f}%)")
        
        logger.info(f"Extraction completed successfully. Saved {len(pdb_files)} PDB files.")
        return pdb_files
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

def extract_folder(archive_path, destination_path, config):
    """
    Extracts a folder from a tar.gz archive to a specified destination.

    Args:
        archive_path (str): The path to the tar.gz archive.
        destination_path (str): The path to the destination directory.
        config (dict): A configuration dictionary containing a logger.
    """
    logger = config['logger']
    logger.info(f"Extracting folder from: {archive_path} to {destination_path}")

    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            # Get the top-level folder name from the archive.
            top_level_folder = tar.getnames()[0].split('/')[0]

            # Extract the entire archive to the destination.
            tar.extractall(path=destination_path)

            # Move the extracted top-level folder to the destination.
            extracted_folder_path = os.path.join(destination_path, top_level_folder)
            final_destination = os.path.join(destination_path, os.path.basename(top_level_folder))

            if extracted_folder_path != final_destination: # only move if they are different.
              shutil.move(extracted_folder_path, final_destination)

        logger.info(f"Extraction completed successfully to: {final_destination}")

    except FileNotFoundError:
        logger.error(f"Error: Archive not found at {archive_path}")
    except tarfile.ReadError:
        logger.error(f"Error: Invalid tar.gz archive at {archive_path}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

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
        
        # Compiling results into a CSV
        logger.info("Compiling pocket data into a CSV")
        process_predicted.merge_to_csv(output_dir, config)

        if compress:
            compress_folder(p2rank_output_dir, f"{p2rank_output_dir}.tar.gz", config)
            return f"{p2rank_output_dir}.tar.gz"
        else:
            return output_dir
        
    except Exception as e:
        logger.error(f"Error during P2Rank execution: {str(e)}")
        raise
