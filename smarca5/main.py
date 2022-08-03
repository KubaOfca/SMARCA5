"""Main module."""
import pandas as pd  # type: ignore


def read_ref_seq_fasta() -> None:
    ref_seq_fasta = ""
    with open(r"..\pliki-krok\SMARCA5_fasta.txt", "r") as file:
        for line in file:
            if line[0] != ">":
                ref_seq_fasta += line.strip()


def create_lps(peptide: str, peptide_len: int) -> list[int]:
    lps = [0] * peptide_len
    prefix_len = 0
    i = 1

    while i < peptide_len:
        if peptide[prefix_len] == peptide[i]:
            prefix_len += 1
            lps[i] = prefix_len
            i += 1
        elif prefix_len != 0:
            prefix_len = lps[prefix_len - 1]
        else:
            lps[i] = 0
            i += 1

    return lps


def search_peptide_in_protein_seq(peptide: str, protein_seq: str) -> None:
    n = len(protein_seq)
    m = len(peptide)

    lps = create_lps(peptide, m)
    i = 0
    j = 0

    while (n - i) >= (m - j):
        if j == m:
            print(f"found at {i - j}")
            j = lps[j - 1]
        elif protein_seq[i] == peptide[j]:
            i += 1
            j += 1
        elif j > 0:
            j = lps[j - 1]
        else:
            i += 1


result_of_experiment = pd.read_excel(
    r"..\pliki-krok\SMARCA5.xlsx", header=None, sheet_name="SMARCA5"
)
