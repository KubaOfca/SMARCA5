"""Main module."""
import bs4  # type: ignore
import colour  # type: ignore
import pandas as pd  # type: ignore
import typing
from typing import List, Dict, Tuple, Set
import ttkbootstrap as tk  # type: ignore
from tkinter import filedialog as fd
import webbrowser
from bs4 import BeautifulSoup
from colour import Color
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl  # type: ignore
import re
import os

import alignment
import html_report
import html_page

# flake8: noqa E203
LINE_LENGTH_DISPLAYED = 80
SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES = 3
LEFT_TEXT_OFFSET = 3
UNIQUE_PATTERN = re.compile(r".*Unique: (.*)\n.*")


# TODO: Better name convention of .html files
# TODO: function to handle operations from 86line
# TODO: GUI as class


def load_file(file_name: str) -> str:
    """
    Load file in correct format for creating exe one-file.

    :param: file_name: filename to load
    :return: str of correct path to file
    """
    return os.path.join(os.path.dirname(__file__), file_name)


def read_ref_seq_fasta() -> str:
    """
    Read reference protein sequence in fasta format.

    :return: reference sequence string
    """
    ref_seq_fasta = ""
    with open(rf"{protein_entry.get()}", "r") as file:
        for line in file:
            if line[0] != ">":
                ref_seq_fasta += line.strip()
    return ref_seq_fasta


# Main function - program start
def generate_alignment_report(dataset, type_of_experiment, alignment_obj, html_report_config):
    for index, experiment in enumerate(dataset):
        if experiment != "All":
            peptides_metadata_filtered_by_experiment = alignment_obj.peptides_metadata.loc[
                alignment_obj.peptides_metadata.loc[:, "Experiment"].str.contains(rf'{experiment}')]
        else:
            peptides_metadata_filtered_by_experiment = alignment_obj.peptides_metadata

        page = html_page.ReportPage(html_report_config, type_of_experiment, experiment,
                                    peptides_metadata_filtered_by_experiment, alignment_obj)
        page.fill_alignment_window()
        page.color_bar_fig.savefig(
            load_file(rf"html_files\Amount_of_peptide-type_{type_of_experiment}-{experiment}.png"))
        img_tag = page.soup.find("img")
        img_tag["src"] = f"Amount_of_peptide-type_{type_of_experiment}-{experiment}.png"
        with open(load_file(rf"html_files\type_{type_of_experiment}-{experiment}.html"), "w") as page_to_save:
            page_to_save.write(page.soup.prettify(formatter="html"))  # soup.encode(formatter="html")


@typing.no_type_check
def find_peptide_in_protein_seq() -> None:
    """
    Find position of peptide in protein sequence.

    :return: None
    """
    global LINE_LENGTH_DISPLAYED
    LINE_LENGTH_DISPLAYED = int(spin_box.get())
    protein_fasta = read_ref_seq_fasta()
    peptides_dataframe = pd.read_excel(rf"{peptide_entry.get()}").loc[
                         :, ["Sequence", "Proteins", "Experiment"]
                         ]
    peptides_dataframe["Start"] = ""
    peptides_dataframe["End"] = ""

    alignment_obj = alignment.Alignment(protein_seq=protein_fasta, peptides_metadata=peptides_dataframe)
    alignment_obj.search_peptide_in_protein_seq()

    with open(load_file(r"html_files\index.html"), "r") as html:
        soup = BeautifulSoup(html, "html.parser")

    # config global options of report like interline, offset ect.
    html_report_config = html_report.HtmlReport(soup,
                                                alignment_obj.group_names,
                                                alignment_obj.sample_names,
                                                LEFT_TEXT_OFFSET,
                                                SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES,
                                                LINE_LENGTH_DISPLAYED
                                                )

    # generate individual reports divided by groups or samples of experiment
    generate_alignment_report(alignment_obj.group_names, "group", alignment_obj, html_report_config)
    generate_alignment_report(alignment_obj.sample_names, "sample", alignment_obj, html_report_config)
    # open one tab (type All)
    webbrowser.open_new_tab(load_file(rf"html_files\type_group-All.html"))

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
                ("All files", "*.*"),
                ("Excel files", ".xlsx .xls"),
                ("CSV files", ".csv"),
                ("Fasta files", ".fasta"),
                ("TXT files", ".txt"),
            ]
        ),
    )


root = tk.Window(themename="darkly")
root.title("SMARCA5")
root.wm_iconbitmap(load_file(r"icons\protein_app.ico"))

title_label = tk.Label(root, text="Peptide position finder in protein")
protein_ref_frame = tk.LabelFrame(
    root, text="Select path of protein sequence fasta file"
)
protein_entry = tk.Entry(protein_ref_frame, width=80)
protein_browse_button = tk.Button(
    protein_ref_frame, text="Browse", command=lambda: browse(protein_entry)
)
peptide_frame = tk.LabelFrame(root, text="Select path of peptide sequence file")
peptide_entry = tk.Entry(peptide_frame, width=80)
peptide_browse_button = tk.Button(
    peptide_frame, text="Browse", command=lambda: browse(peptide_entry)
)
button_frame = tk.Frame(root)
analyze_button = tk.Button(
    button_frame,
    text="Analyze",
    bootstyle="danger",
    command=find_peptide_in_protein_seq,
)

spin_box = tk.Spinbox(button_frame, from_=1, to=1000)
spin_box.insert(0, "80")

sping_box_label = tk.Label(button_frame, text="Line length:")

title_label.pack(padx=20, pady=20)
protein_ref_frame.pack(padx=20, pady=20)
protein_entry.grid(row=0, column=0, padx=20, pady=20)
protein_browse_button.grid(row=0, column=1, padx=20, pady=20)
peptide_frame.pack(padx=20, pady=20)
peptide_entry.grid(row=0, column=0, padx=20, pady=20)
peptide_browse_button.grid(row=0, column=1, padx=20, pady=20)
button_frame.pack(padx=20, pady=20)
sping_box_label.grid(row=0, column=0, padx=10, pady=20)
spin_box.grid(row=0, column=1, padx=10, pady=20)
analyze_button.grid(row=1, column=0, columnspan=2, padx=20, pady=20)
root.mainloop()
