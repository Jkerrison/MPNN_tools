from Bio import SeqIO
from Bio.AlignIO import read
from Bio.Align import AlignInfo
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import json
import os

def parse_fasta(fasta_file):
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return record.seq

def parse_alignment(alignment_file):
    with open(alignment_file, "r") as file:
        alignment = read(file, "fasta")
    return alignment

def create_pssm(alignment):
    summary_align = AlignInfo.SummaryInfo(alignment)
    pssm = summary_align.pos_specific_score_matrix()
    return pssm

def pssm_to_dataframe(pssm):
    pssm_dict = {str(pos): {aa: score for aa, score in scores.items()} for pos, scores in enumerate(pssm, start=1)}
    df = pd.DataFrame(pssm_dict).T
    return df

def normalize_dataframe(df):
    scaler = MinMaxScaler(feature_range=(0, 2))
    normalized_df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns, index=df.index)
    return normalized_df

def main(fasta_file, alignment_file, output_file):
    protein_sequence = parse_fasta(fasta_file)
    alignment = parse_alignment(alignment_file)
    pssm = create_pssm(alignment)
    df = pssm_to_dataframe(pssm)
    normalized_df = normalize_dataframe(df)

    print(normalized_df)
    bias_matrix = normalized_df.values.tolist()
    with open(f'{output_file}.jsonl', 'w') as outfile:
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = os.path.basename(fasta_file)
            sequence = str(record.seq)
            chain = "A"  # Assuming single chain 'A', modify as needed
            bias_matrix = bias_matrix
            entry = {
                name: {
                    chain: bias_matrix
                }
            }
            json.dump(entry, outfile)
            outfile.write('\n')

    normalized_df.to_excel(f'{output_file}.xlsx', index_label="Position")

if __name__ == "__main__":
    fasta_file = "path/to/reference/fasta.fa"
    alignment_file = "path/to/alignment-aln.fa"
    output_file = "path/to/output"
    main(fasta_file, alignment_file, output_file)
