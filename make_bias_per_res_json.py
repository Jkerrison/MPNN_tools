import json
from Bio import SeqIO
import pandas as pd
import sys

#AA ordering from MPNN code
mpnn_alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
mpnn_alphabet_dict = {'A': 0,'C': 1,'D': 2,'E': 3,'F': 4,'G': 5,'H': 6,'I': 7,'K': 8,'L': 9,'M': 10,'N': 11,'P': 12,'Q': 13,'R': 14,'S': 15,'T': 16,'V': 17,'W': 18,'Y': 19,'X': 20}

def generate_bias_matrix(sequence):
    #Make a matrix of default bias (1) NxM where M is the length of the sequence and N is the number of residues
    bias_matrix = []
    for i, residue in enumerate(sequence):
        bias_row = [1.0] * 21
        bias_matrix.append(bias_row)
    bias_df = pd.DataFrame(bias_matrix,columns= list(mpnn_alphabet))
    bias_df.index = bias_df.index + 1

    #simply set a cell position to a bias value
    bias_df.at[2, 'A']+=0.39

    #Set all positions where the residue is 'T' in starting seq to be bias towards 'T'
    for i, residue in enumerate(sequence):
        if residue == 'T' or residue == 'S':
            bias_df.at[i + 1, 'T'] +=0.39
            bias_df.at[i + 1, 'S'] +=0.39

    # Add a polar bias to all positions 
    columns_to_update = ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y']
    for column in columns_to_update:
        bias_df[column] += 0.39

    bias_df['C'] -= 1

    #Here we have the sequences of 50% conservation, a bias is added to these positions
    fifty_percent = '-LO--R-UM-IPS---UM--' #sequence of 50% conservation
    for i, char in enumerate(fifty_percent):
        if char != '-'or char != 'X' or char != '_':
            bias_df.at[i + 1, char] += 0.39

    #Convert back to matrix for addition to json
    bias_matrix = bias_df.values.tolist()
    print(bias_df)
    return bias_matrix


def fasta_to_jsonl(fasta_file, output_file):
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id
            sequence = str(record.seq)
            chain = "A"  # Assuming single chain 'A', modify as needed
            bias_matrix = generate_bias_matrix(sequence)
            entry = {
                name: {
                    chain: bias_matrix
                }
            }
            json.dump(entry, outfile)
            outfile.write('\n')

def main():
    fasta_file = sys.argv[1]
    json_output_file = sys.argv[2]
    fasta_to_jsonl(fasta_file, json_output_file)

if __name__ == "__main__":
    main()