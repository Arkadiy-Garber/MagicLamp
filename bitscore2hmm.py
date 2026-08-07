#!/usr/bin/env python3
import os
import argparse
import pandas as pd


def update_hmm_files(bitscore_file, hmm_folder):
    """Creates new HMM files with TC values inside the header, below the CKSUM line."""
    # Load bitscore data
    bitscore_df = pd.read_csv(bitscore_file)

    # Iterate over each row in the bitscore file
    for _, row in bitscore_df.iterrows():
        hmm_name = str(row[0])  # First column: HMM filename (without .hmm extension)
        hmm_name_short = hmm_name.split(".")[0]
        bit_score = float(row[1])  # Second column: Bit score

        hmm_path = os.path.join(hmm_folder, f"{hmm_name}")
        new_hmm_path = os.path.join(hmm_folder, f"{hmm_name_short}.TC.hmm")

        if not os.path.exists(hmm_path):
            print(f"Warning: HMM file {hmm_path} not found.")
            continue

        # Read the HMM file
        with open(hmm_path, "r") as f:
            lines = f.readlines()

        updated_lines = []
        cksum_found = False

        for i, line in enumerate(lines):
            updated_lines.append(line)

            # Insert TC line **right after** the CKSUM line
            if line.startswith("CKSUM"):
                updated_lines.append(f"TC    {bit_score} {bit_score - 0.9};\n")
                cksum_found = True

        # If CKSUM line was not found, print a warning
        if not cksum_found:
            print(f"Warning: No CKSUM line found in {hmm_name}.hmm, TC line was not inserted.")

        # Write the new HMM file
        with open(new_hmm_path, "w") as f:
            f.writelines(updated_lines)

        print(f"Created {hmm_name_short}.TC.hmm with TC    {bit_score} {bit_score - 0.9};")

    print("Processing complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create new HMM files with TC values inserted after the CKSUM line.")
    parser.add_argument("--bitscore", required=True, help="Path to the bitscores.csv file.")
    parser.add_argument("--hmm_folder", required=True, help="Path to the folder containing HMM files.")

    args = parser.parse_args()

    update_hmm_files(args.bitscore, args.hmm_folder)