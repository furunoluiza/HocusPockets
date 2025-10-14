import os
import csv
import pandas as pd

def merge_cavity_scores(struct_dir, reference_csv, output_csv):
    """
    Combines all *_cavity_scores.csv files found in struct_dir into a single table,
    adding the 'Category' column from a reference CSV based on Protein_Code matching.
    
    Parameters
    ----------
    struct_dir : str
        Path to the directory containing protein subfolders with *_cavity_scores.csv files.
    reference_csv : str
        Path to the CSV file containing columns ['Protein_Code', 'Category'].
    output_csv : str
        Path to save the merged output CSV file.
    """
    
    # Load the reference file (must contain 'Protein_Code' and 'Category')
    ref_df = pd.read_csv(reference_csv)
    if 'Protein_Code' not in ref_df.columns or 'Category' not in ref_df.columns:
        raise ValueError("The reference CSV must contain 'Protein_Code' and 'Category' columns.")
    
    all_rows = []

    # Iterate through all protein directories in struct_dir
    for name in os.listdir(struct_dir):
        protein_dir = os.path.join(struct_dir, name)
        if not os.path.isdir(protein_dir):
            continue
        
        cav_dir = os.path.join(protein_dir, "cavities")
        score_file = os.path.join(cav_dir, f"{name}_cavity_scores.csv")

        if not os.path.isfile(score_file):
            print(f"[WARNING] Score file not found for {name}: {score_file}")
            continue

        # Read the score file
        try:
            df = pd.read_csv(score_file)
        except Exception as e:
            print(f"[ERROR] Could not read {score_file}: {e}")
            continue

        # Normalize possible column names
        if "Protein_Code" not in df.columns and "Protein" in df.columns:
            df.rename(columns={"Protein": "Protein_Code"}, inplace=True)
        
        # Find the category for this protein
        protein_entry = ref_df.loc[ref_df["Protein_Code"] == name]
        if protein_entry.empty:
            print(f"[INFO] Protein {name} not found in reference CSV — skipping.")
            continue

        category_value = protein_entry["Category"].values[0]
        df["Category"] = category_value

        all_rows.append(df)

    # Concatenate all results
    if all_rows:
        merged_df = pd.concat(all_rows, ignore_index=True)
        merged_df.to_csv(output_csv, index=False)
        print(f"[INFO] Combined file saved successfully → {output_csv}")
        print(f"Total entries: {len(merged_df)}")
    else:
        print("[WARNING] No valid _cavity_scores.csv files were found or matched.")

# Example usage:
merge_cavity_scores(
    struct_dir="/path/to/results"
    reference_csv="/path/to/reference.csv",
    output_csv="all_cavities_with_category.csv"
)
