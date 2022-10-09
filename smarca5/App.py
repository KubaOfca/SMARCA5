import ttkbootstrap as tk  # type: ignore
import bs4  # type: ignore
import colour  # type: ignore
import pandas as pd  # type: ignore
import ttkbootstrap as tk  # type: ignore
from tkinter import filedialog as fd
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl  # type: ignore
import os
import webbrowser
from bs4 import BeautifulSoup
import atexit
import glob

import alignment
import html_report
import html_page

LEFT_TEXT_OFFSET = 3


class App:
    def __init__(self, master):
        atexit.register(delete_html_files)
        master.title("SMARCA5")
        master.wm_iconbitmap(load_file(r"icons\protein_app.ico"))
        self.title_label = tk.Label(
            master, text="Peptide position finder in protein"
        )
        self.protein_ref_frame = tk.LabelFrame(
            master, text="Select path of protein sequence fasta file"
        )
        self.protein_entry = tk.Entry(self.protein_ref_frame, width=80)
        self.protein_browse_button = tk.Button(
            self.protein_ref_frame,
            text="Browse",
            command=lambda: self.browse(self.protein_entry),
        )
        self.peptide_frame = tk.LabelFrame(
            master, text="Select path of peptide sequence file"
        )
        self.peptide_entry = tk.Entry(self.peptide_frame, width=80)
        self.peptide_browse_button = tk.Button(
            self.peptide_frame,
            text="Browse",
            command=lambda: self.browse(self.peptide_entry),
        )
        self.button_frame = tk.Frame(master)
        self.analyze_button = tk.Button(
            self.button_frame,
            text="Analyze",
            bootstyle="danger",
            command=self.find_peptide_in_protein_seq,
        )

        self.line_length_spinbox = tk.Spinbox(
            self.button_frame, from_=1, to=10000
        )
        self.line_length_spinbox.insert(0, "80")
        self.line_length_spinbox_label = tk.Label(
            self.button_frame, text="Line length (character):"
        )

        self.interline_spinbox = tk.Spinbox(
            self.button_frame, from_=2, to=10000
        )
        self.interline_spinbox.insert(0, "20")
        self.interline_spinbox_label = tk.Label(
            self.button_frame, text="Interline (pixel):"
        )

        self.title_label.pack(padx=20, pady=20)
        self.protein_ref_frame.pack(padx=20, pady=20)
        self.protein_entry.grid(row=0, column=0, padx=20, pady=20)
        self.protein_browse_button.grid(row=0, column=1, padx=20, pady=20)
        self.peptide_frame.pack(padx=20, pady=20)
        self.peptide_entry.grid(row=0, column=0, padx=20, pady=20)
        self.peptide_browse_button.grid(row=0, column=1, padx=20, pady=20)
        self.button_frame.pack(padx=20, pady=20)
        self.line_length_spinbox_label.grid(row=0, column=0, padx=10, pady=20)
        self.line_length_spinbox.grid(row=0, column=1, padx=10, pady=20)
        self.interline_spinbox_label.grid(row=1, column=0, padx=10, pady=20)
        self.interline_spinbox.grid(row=1, column=1, padx=10, pady=20)
        self.analyze_button.grid(
            row=2, column=0, columnspan=2, padx=20, pady=20
        )

    def read_ref_seq_fasta(self) -> str:
        """
        Read reference protein sequence in fasta format.

        :return: reference sequence string
        """
        ref_seq_fasta = ""
        with open(rf"{self.protein_entry.get()}", "r") as file:
            for line in file:
                if line[0] != ">":
                    ref_seq_fasta += line.strip()
        return ref_seq_fasta

    def browse(self, entry) -> None:
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

    def generate_alignment_report(
        self, dataset, type_of_experiment, alignment_obj, html_report_config
    ):
        for index, experiment in enumerate(dataset):
            if experiment != "All":
                peptides_metadata_filtered_by_experiment = (
                    alignment_obj.peptides_metadata.loc[
                        alignment_obj.peptides_metadata.loc[
                            :, "Experiment"
                        ].str.contains(rf"{experiment}")
                    ]
                )
            else:
                peptides_metadata_filtered_by_experiment = (
                    alignment_obj.peptides_metadata
                )

            page = html_page.ReportPage(
                html_report_config,
                type_of_experiment,
                experiment,
                peptides_metadata_filtered_by_experiment,
                alignment_obj,
            )
            page.fill_alignment_window()
            page.color_bar_fig.savefig(
                load_file(
                    rf"html_files\Amount_of_peptide-type_{type_of_experiment}-{experiment}.png"
                )
            )
            plt.cla()
            plt.close(page.color_bar_fig)
            img_tag = page.soup.find("img")
            img_tag[
                "src"
            ] = f"Amount_of_peptide-type_{type_of_experiment}-{experiment}.png"
            with open(
                load_file(
                    rf"html_files\type_{type_of_experiment}-{experiment}.html"
                ),
                "w",
            ) as page_to_save:
                page_to_save.write(
                    page.soup.prettify(formatter="html")
                )  # soup.encode(formatter="html")

    def find_peptide_in_protein_seq(self) -> None:
        """
        Find position of peptide in protein sequence.

        :return: None
        """
        protein_fasta = self.read_ref_seq_fasta()
        peptides_dataframe = pd.read_excel(rf"{self.peptide_entry.get()}").loc[
            :, ["Sequence", "Proteins", "Experiment"]
        ]
        peptides_dataframe["Start"] = ""
        peptides_dataframe["End"] = ""

        alignment_obj = alignment.Alignment(
            protein_seq=protein_fasta, peptides_metadata=peptides_dataframe
        )
        alignment_obj.search_peptide_in_protein_seq()

        with open(load_file(r"html_files\index.html"), "r") as html:
            soup = BeautifulSoup(html, "html.parser")

            # config global options of report like interline, offset ect.
            html_report_config = html_report.HtmlReport(
                soup,
                alignment_obj.group_names,
                alignment_obj.sample_names,
                LEFT_TEXT_OFFSET,
                int(self.interline_spinbox.get()),
                int(self.line_length_spinbox.get()),
            )

            # generate individual reports divided by groups or samples of experiment
            self.generate_alignment_report(
                alignment_obj.group_names,
                "group",
                alignment_obj,
                html_report_config,
            )
            self.generate_alignment_report(
                alignment_obj.sample_names,
                "sample",
                alignment_obj,
                html_report_config,
            )
            # open one tab (type All)
            webbrowser.open_new_tab(
                load_file(r"html_files\type_group-All.html")
            )


def load_file(file_name: str) -> str:
    """
    Load file in correct format for creating exe one-file.

    :param: file_name: filename to load
    :return: str of correct path to file
    """
    return os.path.join(os.path.dirname(__file__), file_name)


def delete_html_files():
    print("Delete html files")
    files = glob.glob(load_file(r"html_files\*"))
    for f in files:
        if f != load_file(r"html_files\index.html") and f != load_file(
            r"html_files\style.css"
        ):
            os.remove(f)
