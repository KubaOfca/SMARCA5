"""Main module."""
import bs4  # type: ignore
import colour  # type: ignore
import pandas as pd  # type: ignore
from typing import List, Dict, Tuple, Set
import ttkbootstrap as tk  # type: ignore
from tkinter import filedialog as fd
import webbrowser
from bs4 import BeautifulSoup
from colour import Color
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl  # type: ignore


def read_ref_seq_fasta() -> str:
    """
    Read reference protein sequence in fasta format.

    :return: reference sequence string
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
    result_of_experiment: pd.DataFrame,
    protein_seq: str,
) -> Tuple[List[pd.Series], Dict[str, int]]:
    """
    KMP (Knuth–Morris–Pratt algorithm) for pattern searching.

    :param result_of_experiment: dataframe with result of experiment
    :param protein_seq: text where we want to find pattern
    :return: List[List[pd.Series], Dict[str, int]]
    """
    result: List[pd.Series] = []  # coord 0 is first letter in string
    amount: Dict[str, int] = {}
    for index, peptide in result_of_experiment.iterrows():
        n = len(protein_seq)
        m = len(peptide.Sequence)

        lps = create_lps(peptide.Sequence, m)
        i = 0
        j = 0

        while (n - i) >= (m - j):
            if j == m:
                coords = (i - j, i)
                peptide["Coords"] = coords
                result.append(peptide)
                if peptide["Sequence"] in amount:
                    amount[peptide["Sequence"]] += 1
                else:
                    amount[peptide["Sequence"]] = 1
                j = lps[j - 1]
            elif protein_seq[i] == peptide.Sequence[j]:
                i += 1
                j += 1
            elif j > 0:
                j = lps[j - 1]
            else:
                i += 1

    result.sort(key=lambda x: x["Coords"])
    return result, amount


def create_table_of_information_about_peptide(
    result: List[pd.Series],
    current_peptide_seq: str,
    amount: Dict[str, int],
    experiment_types: Set[str],
) -> str:
    """
    Create tooltip with extra information about peptide.

    :param result: peptide list
    :param current_peptide_seq: peptide sequence
    about which we want to find additional information
    :param amount: dictionary where key is sequence and value is amount
    :param experiment_types: list of all found experiment types
    :return: tooltip formated text about peptide extra information
    """
    peptide_proteins = []
    peptide_experiments = []
    not_existing_experiments = []
    for peptide in result:
        if peptide["Sequence"] == current_peptide_seq:
            if peptide["Proteins"] not in peptide_proteins:
                peptide_proteins.append(peptide["Proteins"])
            if peptide["Experiment"] not in peptide_experiments:
                peptide_experiments.append(peptide["Experiment"])
    for experiment in experiment_types:
        if experiment not in peptide_experiments:
            not_existing_experiments.append(experiment)
    proteins_str = ", ".join(peptide_proteins).replace(";", ", ")
    experiment_str = ", ".join(peptide_experiments)
    not_existing_experiment_str = ", ".join(not_existing_experiments)
    return (
        f"Proteins: {proteins_str}\n"
        f"Experiment: {experiment_str}\n"
        f"Experiments not found: {not_existing_experiment_str}\n"
        f"Amount: {amount[current_peptide_seq]}"
    )


def create_color_bar_image(
    amount: Dict[str, int], amount_list: List[str]
) -> List[colour.Color]:
    """
    Create color bar image.

    :param amount: dictionary where key is sequence and value is amount
    :param amount_list: peptide sequence list in descending order
    :return: list of colors used in color bar
    """
    red = Color("red")
    colors = list(red.range_to(Color("pink"), len(amount)))

    fig, ax = plt.subplots(figsize=(2, 4))
    fig.subplots_adjust(right=0.5)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "", [x.hex for x in reversed(colors)]
    )
    norm = mpl.colors.Normalize(
        vmin=amount[amount_list[-1]], vmax=amount[amount_list[0]]
    )

    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)
    cb1.set_label("Amount of peptide")
    fig.savefig("Amount_of_peptide.png")

    return colors


def add_peptide_on_html_page(
    soup: bs4.BeautifulSoup,
    peptide: pd.Series,
    colors: List[colour.Color],
    amount_list: List[str],
    tooltip_text: str,
    alignment: str,
) -> None:
    """
    Add peptide to div in order to display it on page.

    :param soup: html page handler
    :param peptide: peptide information (pd.Series)
    :param colors: list of color bar colors
    :param amount_list: peptide sequence list in descending order
    :param tooltip_text: string containing peptide extra information
    :param alignment: peptide sequence alignment to ref sequence
    :return: None
    """
    div_center = soup.find("div", {"class": "center"})
    p_tag = soup.new_tag("p")
    peptide_seq = peptide["Sequence"]
    p_tag[
        "style"
    ] = f"color: {colors[amount_list.index(f'{peptide_seq}')].hex}"
    p_tag["title"] = tooltip_text
    div_center.append(p_tag)
    p_tag.append(soup.new_string(alignment))


def create_html_report(
    result: List[pd.Series],
    ref_seq_fasta: str,
    amount: Dict[str, int],
    amount_list: List[str],
    all_possible_experiment_type: Set[str],
    colors: List[colour.Color],
) -> None:
    """
    Create html report.

    :param result: peptide list
    :param ref_seq_fasta: reference sequence of protein
    :param amount: dictionary where key is sequence and value is amount
    :param amount_list: peptide sequence list in descending order
    :param all_possible_experiment_type: list of all found experiment types
    :param colors: list of color bar colors
    :return: None
    """
    with open("index.html", "r") as html:
        soup = BeautifulSoup(html, "html.parser")
    soup.div.p.string = ref_seq_fasta
    current_seq = ""

    for peptide in result:
        if peptide["Sequence"] == current_seq:
            continue
        else:
            current_seq = peptide["Sequence"]
            tooltip_text = create_table_of_information_about_peptide(
                result, current_seq, amount, all_possible_experiment_type
            )
        alignment = "-" * len(ref_seq_fasta)
        alignment = (
            alignment[: peptide["Coords"][0]]
            + current_seq
            + alignment[peptide["Coords"][1] :]  # noqa: E203
        )
        alignment = alignment.replace("-", "\xa0")
        add_peptide_on_html_page(
            soup, peptide, colors, amount_list, tooltip_text, alignment
        )

    # save changes to page.html and display page with result
    with open("page.html", "w") as page:
        page.write(
            soup.prettify(formatter="html")
        )  # soup.encode(formatter="html")
        webbrowser.open_new_tab("page.html")


def find_peptide_in_protein_seq() -> None:
    """
    Find position of peptide in protein sequence.

    :return:
    """
    ref_seq_fasta = read_ref_seq_fasta()
    result_of_experiment = pd.read_excel(
        r"..\pliki-krok\SMARCA5.xlsx", sheet_name="SMARCA5"
    ).loc[:, ["Sequence", "Proteins", "Experiment"]]

    result, amount = search_peptide_in_protein_seq(
        result_of_experiment, ref_seq_fasta
    )
    amount_list = sorted(amount, key=amount.get, reverse=True)  # type: ignore
    all_possible_experiment_type = {x["Experiment"] for x in result}
    colors = create_color_bar_image(amount, amount_list)
    create_html_report(
        result,
        ref_seq_fasta,
        amount,
        amount_list,
        all_possible_experiment_type,
        colors,
    )


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
analyze_button = tk.Button(
    root,
    text="Analyze",
    bootstyle="danger",
    command=find_peptide_in_protein_seq,
)

title_label.pack(padx=20, pady=20)
protein_ref_frame.pack(padx=20, pady=20)
protein_entry.grid(row=0, column=0, padx=20, pady=20)
protein_browse_button.grid(row=0, column=1, padx=20, pady=20)
peptide_frame.pack(padx=20, pady=20)
peptide_entry.grid(row=0, column=0, padx=20, pady=20)
peptide_browse_button.grid(row=0, column=1, padx=20, pady=20)
analyze_button.pack(padx=20, pady=20)

root.mainloop()
