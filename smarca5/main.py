"""Main module."""
import pandas as pd  # type: ignore
from typing import Dict, Tuple
import ttkbootstrap as tk  # type: ignore
from tkinter import filedialog as fd


def read_ref_seq_fasta() -> str:
    """
    Read reference protein sequence in fasta format.

    :return: str
    """
    ref_seq_fasta = ""
    with open(r"..\pliki-krok\SMARCA5_fasta.txt", "r") as file:
        for line in file:
            if line[0] != ">":
                ref_seq_fasta += line.strip()
    return ref_seq_fasta


def create_lps(peptide: str, peptide_len: int) -> list[int]:
    """
    Create lps table.

    :param peptide: peptide sequence
    :param peptide_len: length of peptide sequence
    :return: list containing integers which store information
    about length of prefix
    """
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


def search_peptide_in_protein_seq(
    peptide: pd.Series,
    protein_seq: str,
    result: Dict[Tuple[int, int], pd.Series],
) -> None:
    """
    KMP (Knuth–Morris–Pratt algorithm) for pattern searching.

    :param result:
    :param peptide: pattern
    :param protein_seq: text where we want to find pattern
    :return: None
    """
    n = len(protein_seq)
    m = len(peptide.Sequence)

    lps = create_lps(peptide.Sequence, m)
    i = 0
    j = 0

    while (n - i) >= (m - j):
        if j == m:
            coords = (i - j, i)
            if coords in result:
                result[coords]["Amount"] += 1
            else:
                result[coords] = peptide
                result[coords]["Amount"] += 1
            j = lps[j - 1]
        elif protein_seq[i] == peptide.Sequence[j]:
            i += 1
            j += 1
        elif j > 0:
            j = lps[j - 1]
        else:
            i += 1


def find_peptide_in_protein_seq() -> None:
    """
    Find position of peptide in protein sequence.

    :return:
    """
    ref_seq_fasta = read_ref_seq_fasta()
    result_of_experiment = pd.read_excel(
        r"..\pliki-krok\SMARCA5.xlsx", sheet_name="SMARCA5"
    ).loc[:, ["Sequence", "Proteins", "Experiment"]]
    result: Dict[
        Tuple[int, int], pd.Series
    ] = {}  # coord 0 is first letter in string
    for index, row in result_of_experiment.iterrows():
        row["Amount"] = 0
        search_peptide_in_protein_seq(
            peptide=row, protein_seq=ref_seq_fasta, result=result
        )
    print(result)


# GUI


def browse(entry: tk.Entry) -> None:
    """
    Browse button behavior.

    :param entry: tk.Entry object where selected file path will be displayed
    :return:
    """
    entry.delete(0, tk.END)
    entry.insert(
        0,
        fd.askopenfilename(
            filetypes=[
                ("Excel files", ".xlsx .xls"),
                ("CSV files", ".csv"),
                ("Fasta files", ".fasta"),
                ("TXT files", ".txt"),
                ("All files", "*.*"),
            ]
        ),
    )


root = tk.Window(themename="darkly")

title_label = tk.Label(root, text="Peptide position finder in protein")
protein_ref_frame = tk.LabelFrame(
    root, text="Select path of protein sequence fasta file"
)
protein_entry = tk.Entry(protein_ref_frame, width=80)
protein_browse_button = tk.Button(
    protein_ref_frame, text="Browse", command=lambda: browse(protein_entry)
)
peptide_frame = tk.LabelFrame(
    root, text="Select path of peptide sequence file"
)
peptide_entry = tk.Entry(peptide_frame, width=80)
peptide_browse_button = tk.Button(
    peptide_frame, text="Browse", command=lambda: browse(peptide_entry)
)
analyze_button = tk.Button(root, text="Analyze", bootstyle="danger")

title_label.pack(padx=20, pady=20)
protein_ref_frame.pack(padx=20, pady=20)
protein_entry.grid(row=0, column=0, padx=20, pady=20)
protein_browse_button.grid(row=0, column=1, padx=20, pady=20)
peptide_frame.pack(padx=20, pady=20)
peptide_entry.grid(row=0, column=0, padx=20, pady=20)
peptide_browse_button.grid(row=0, column=1, padx=20, pady=20)
analyze_button.pack(padx=20, pady=20)

root.mainloop()
