import os
import csv
import re
import glob

def extract_frame_number(filename):
    """Extracts the number before '.pdb' and after '_' from a filename."""
    match = re.search(r'_(\d+)\.pdb', filename)

    if match:
        return match.group(1)
    else:
        return None

def load_pdb_list(pdb_list_file):
    with open(pdb_list_file, 'r') as f:
        # Assuming each line in pdb_list.ds is a PDB file path
        pdb_files = [line.strip() for line in f.readlines()]
    return pdb_files

def merge_to_csv(outfolder, config):

    logger = config['logger']

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

        # Get a list of all prediction files in the output folder.
        p2rank_output_folder = os.path.join(outfolder, 'p2rank_output')
        prediction_list = glob.glob(p2rank_output_folder+'/*.pdb_predictions.csv')

        # Iterate through PDB files
        for predictions_file in prediction_list:
            
            # Skip if predictions file doesn't exist
            if not os.path.exists(predictions_file):
                logger.info(f"Predictions file {predictions_file} not found. Skipping.")
                continue

            # Extract the frame number from the PDB file name
            frame_number = extract_frame_number(os.path.basename(predictions_file))
            if not frame_number:
                logger.info(f"Frame number not found for {predictions_file}. Skipping.")
                continue

            # Open predictions CSV file
            with open(predictions_file, 'r') as pred_file:
                reader = csv.DictReader(pred_file)
                for row in reader:
                    # Check if probability is higher than 0.5
                    # TO-DO: THIS CAN BE MADE ADJUSTABLE BY THE USER!

                    if float(row[' probability']) > 0.5:
                        pocket_index = row['  rank']  # Pocket index
                        full_name = os.path.basename(predictions_file)[:-4]  # Get the PDB file name as the File name
                        probability = row[' probability']
                        residues = row[' residue_ids']

                        # Write data to the CSV file
                        csvwriter.writerow([full_name, frame_number, pocket_index, probability, residues])
    return csv_filename
