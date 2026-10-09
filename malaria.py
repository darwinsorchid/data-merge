"""
┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
File_name: malaria.py
Author: Alexandra Bekou
Date:

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

Description: Adds protein description from a blastx file to the corresponding sequence header in a fasta file.
             The output is printed and saved to a user-specified .txt file.

Usage: python malaria.py malaria.fna malaria.blastx.tab output.txt

Procedure:
    1. Parse command-line arguments (input FASTA file, input BLASTX file, output file)
    2. Read FASTA file sequence by sequence and get ids
    3. Parse BLASTX file and find matching ids
    4. Generate (ids, protein name) tuple
    5. Go through FASTA file and append protein name at the end of the header with the matching id
    6. Print output and save to user-defined output file

Input 1: FASTA (tab-separated header)
Input 2: BLASTX (tab-separated fields)
Output: .txt file with FASTA header with protein name + sequence (tab-separated header)

┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
"""

import sys


input_fasta_path = sys.argv[1]
input_blastx_path = sys.argv[2]
output_file_path = sys.argv[3]


print(f"Fasta file: {input_fasta_path}")
print(f"BlastX file: {input_blastx_path}")
print(f"Output_file: {output_file_path}")


# Create BLASTX look-up ids

blast_dict = {}

with open(input_blastx_path, "r") as blast_in_file:
    # Exclude header
    header = blast_in_file.readline()
    for line in blast_in_file:
        words = line.split("\t")
        # print(words)
        id_key = words[0]
        description = words[9].rstrip()
        blast_dict[id_key] = description

# print(blast_dict)

# Check if BLASTX dict is not empty
if not blast_dict:
    raise ValueError("Parsing BLASTX file resulted in an empty dictionary.")

write_seq = False

# Parse FASTA and look for matches in BLASTX
with (
    open(input_fasta_path, "r") as fasta_in_file,
    open(output_file_path, "w") as out_file,
):
    for line in fasta_in_file:
        # print(line)
        if line.startswith(">"):
            new_words = []
            words = line.split("\t")
            clean_id = words[0].lstrip(">")
            # print(clean_id)
            hit_description = blast_dict[clean_id]
            # print(repr(hit_description))
            if hit_description != "null":
                write_seq = True
                new_words.append(line.rstrip())
                new_words.append("protein=" + hit_description)
                new_line = "\t".join(new_words)
                print(new_line)
                out_file.write(new_line + "\n")
            else:
                write_seq = False

        else:
            if write_seq:
                print(line)
                out_file.write(line)
            else:
                continue
