#!/usr/bin/env python3
import sys
from Bio import SeqIO

def detect_file_format(input_file, output_file):
    """Detect if the file is GenBank (gbk) or FASTA (fa)."""
    out = open(output_file, "w")
    try:
        # Try parsing as GenBank
        with open(input_file, "r") as file:
            for record in SeqIO.parse(file, "genbank"):
                print("gbk")
                out.write("gbk")
                return

        # Try parsing as FASTA
        with open(input_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                print("fa")
                out.write("fa")
                return

    except Exception as e:
        pass

    out.close()
    print("Unknown format")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python detect_format.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    detect_file_format(input_file, output_file)