import os
import csv
import re


def extract_frame_number(filename):
    """Extract frame number from filenames like '2XIR_traj_1.pdb_predictions.csv'."""
    match = re.search(r'(\d+)(?=\.[^.]+$)', filename)

    if match:
        return match.group(1)
    else:
        return None

def load_pdb_list(pdb_list_file):
    with open(pdb_list_file, 'r') as f:
        # Assuming each line in pdb_list.ds is a PDB file path
        pdb_files = [line.strip() for line in f.readlines()]
    return pdb_files

def merge_to_csv(outfolder, pdb_file_list, config):
    """
    Writes the pocket predictions data into a CSV file.
    
    Args:
    - processed_dir (str): Directory where processed files are stored.
    - pdb_list (list): List of PDB file paths.
    
    Returns:
    - None
    """

    logger = config['logger']
    # Load the pdb list from the .ds file
    pdb_list = load_pdb_list(pdb_file_list)
    # Define the output CSV filename
    csv_filename = os.path.join(outfolder, "pockets.csv")
    
    # Ensure the directory exists
    os.makedirs(os.path.dirname(csv_filename), exist_ok=True)

    logger.info('Compiling pocket predictions into summary csv files.')
    # Open CSV file for writing
    with open(csv_filename, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)

        # Write header
        csvwriter.writerow(["File name", 'Frame', "pocket_index", "probability", "residues"])

        # Iterate through PDB files
        for pdb_file in pdb_list:
            #print(f"pdb file:{pdb_file}")
            # Define paths for prediction and residues files
            predictions_file = os.path.join(outfolder, os.path.basename(pdb_file).replace('.pdb', '.pdb_predictions.csv'))
            residues_file = os.path.join(outfolder, os.path.basename(pdb_file).replace('.pdb', '.pdb_residues.csv'))
            #print(f"predictions:{predictions_file} and residues:{residues_file}")

            # Skip if predictions file doesn't exist
            if not os.path.exists(predictions_file):
                continue

            # Extract the frame number from the PDB file name
            frame_number = extract_frame_number(os.path.basename(pdb_file))
            if not frame_number:
                logger.info(f"Frame number not found for {pdb_file}. Skipping.")
                continue

            # Open predictions CSV file
            with open(predictions_file, 'r') as pred_file:
                reader = csv.DictReader(pred_file)
                for row in reader:
                    # Check if probability is higher than 0.5
                    # TO-DO: THIS CAN BE MADE ADJUSTABLE BY THE USER!

                    if float(row[' probability']) > 0.5:
                        pocket_index = row['  rank']  # Pocket index
                        full_name = os.path.basename(pdb_file)[:-4]  # Get the PDB file name as the File name
                        probability = row[' probability']
                        residues = row[' residue_ids']

                        # Write data to the CSV file
                        csvwriter.writerow([full_name, frame_number, pocket_index, probability, residues])
    return csv_filename
