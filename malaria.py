"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

malaria.py

Description: This program merges information from a BLASTX and a FASTA file into a new file and prints the
            updated content. The user provides the three file paths and the program parses the FASTA and
            BLASTX files to match IDs and retrieve protein description information from the latter. The
            output format (printed and saved as a file) follows the given FASTA file and the protein descri-
            ption information is added to the header of each entry. Sequences for which protein description
            is "null" are excluded from the final output.
            sys.exit() is used instead of raise <Error> so as not to be considered crushing


User-defined functions: None
Non-standard modules:
    os: python built-in module for interacting with the operating system - here used to check paths

Procedure:
    1. Parse command-line arguments (input FASTA file, input BLASTX file, output file)
    2. Read BLASTX file and create {ID: hitDescription} dictionary
    3. Read FASTA file and use extracted sequence IDs to retrieve protein description from BLASTX dictionary
    4. Update header lines to include protein description
    5. Exclude header + sequence lines if protein description is "null"
    6. Print output and save to user-defined output file

Input 1: FASTA file (tab-separated header)
Input 2: BLASTX file (tab-separated fields)
Output: .txt file in FASTA format with protein name added at the end of the tab-separated header

Usage: python malaria.py malaria.fna malaria.blastx.tab output.txt

Version: 1.00
Date: 2026-10-12
Name: Alexandra Bekou
--------------------------------

Limitations:
    • Accepts FASTA input file only as .fna
    • Accepts BLASTX input file only as .blastx.tab
    • Input files are only checked for extension and not for content compatibility
    • Duplicate BLASTX entries for the same ID are not checked
    • Sequence length is not checked to confirm correct sequence entry in both files

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

# Import modules
import os
import sys

# Check that correct number of arguments exist
if len(sys.argv) != 4:
    sys.exit(f"Expected 3 arguments, got {len(sys.argv) - 1}.")

# Assign user-defined input and output file paths to variables for downstream use
input_fasta_path = sys.argv[1]
input_blastx_path = sys.argv[2]
output_file_path = sys.argv[3]

# Check that FASTA and BLASTX input file paths exist
if not os.path.isfile(input_fasta_path):
    sys.exit("FASTA file not found.")

if not os.path.isfile(input_blastx_path):
    sys.exit("BLASTX file not found.")

# Check for correct input file formats
if not input_fasta_path.endswith(".fna"):
    sys.exit("Provide valid .fna file.")

if not input_blastx_path.endswith(".blastx.tab"):
    sys.exit("Provide valid .blastx.tab file.")


"""
A. Parse BLASTX file to extract protein description (hitDescription column) and query IDs 
to populate initialized dictionary.
"""

# Initialize BLASTX look-up dictionary
blast_dict = {}

# Open BLASTX input file defined by user
with open(input_blastx_path, "r") as blast_in_file:
    # Exclude header line
    header = blast_in_file.readline()

    # Read rest of BLASTX file
    for line in blast_in_file:
        # Convert tab-separated line to list
        words = line.split("\t")
        # Extract query ID
        id_key = words[0]
        # Extract protein description from hitDescription column
        description = words[9].rstrip()
        # Check for empty string of protein description
        if not description:
            sys.exit(f"Protein description does not exist in BLASTX file for {id_key}.")
        # Populate initialized dictionary with {queryID: hitDescription}
        blast_dict[id_key] = description

# Check that BLASTX dict is not empty
if not blast_dict:
    sys.exit("Parsing BLASTX file resulted in an empty dictionary.")


"""
 B.Parse FASTA file, extract sequence IDs and look for matches in BLASTX dictionary.
 Append protein description to headers of non-null descriptions. Save output to file
 and print updated content.
"""

# Initialize variable as False for excluding null protein description entries from output
write_seq = False

# Open input FASTA file and output file defined by user
with (
    open(input_fasta_path, "r") as fasta_in_file,
    open(output_file_path, "w") as out_file,
):
    # Read FASTA file
    for line in fasta_in_file:
        # Check header lines
        if line.startswith(">"):
            # Initialize new header list (will be joined later and provide the final header)
            new_words = []
            # Convert tab-separated line into list
            words = line.split("\t")
            # Extract sequence ID
            clean_id = words[0].lstrip(">")
            # Retrieve protein description from BLASTX look-up dictionary
            hit_description = blast_dict[clean_id]

            # Check if hit description is null to decide if to exclude sequences from output
            # I. NOT NULL
            if hit_description != "null" and hit_description is not None:
                # Update boolean variable to write sequence with non-null protein description
                write_seq = True
                # Append line and protein description to new header list
                new_words.append(line.rstrip())
                new_words.append("protein=" + hit_description)
                # Finalize updated header of sequence entry
                new_line = "\t".join(new_words)
                # Print updated header to output
                print(new_line)
                # Write updated header to user-defined output file
                out_file.write(new_line + "\n")
            # II. NULL PROTEIN DESCRIPTION
            elif hit_description == "null":
                # Update boolean variable to False in case it was changed to True during previous iterations
                write_seq = False

        # Check sequence lines
        else:
            # Check boolean variable to decide writing sequence to output
            # I. Protein description exists and is not null
            if write_seq:
                # Print sequence and write to output file
                print(line.rstrip())
                out_file.write(line)
            # II. Protein description exists and is null
            else:
                # Skip sequence
                continue
