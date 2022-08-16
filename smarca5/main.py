"""Main module."""
import bs4  # type: ignore
import colour  # type: ignore
import pandas as pd  # type: ignore
import typing
from typing import List, Dict, Tuple, Set, Union
import ttkbootstrap as tk  # type: ignore
from tkinter import filedialog as fd
import webbrowser
from bs4 import BeautifulSoup
from colour import Color
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl  # type: ignore
import re
import numpy
# flake8: noqa E203
LINE_LENGTH_DISPLAYED = 80
SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES = 4
LEFT_TEXT_OFFSET = 33
UNIQUE_PATTERN = re.compile(r".*Unique: (.*)\n.*")

def read_ref_seq_fasta() -> str:
    """
    Read reference protein sequence in fasta format.

    :return: reference sequence string
    """
    ref_seq_fasta = ""
    with open(r"..\pliki-SMARCA5\SMARCA5_fasta.txt", "r") as file:
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
            tmp_peptide_protein_list = peptide["Proteins"].split(";")
            for protein in tmp_peptide_protein_list:
                if protein not in peptide_proteins:
                    peptide_proteins.append(protein)
            if peptide["Experiment"] not in peptide_experiments:
                peptide_experiments.append(peptide["Experiment"])
    for experiment in experiment_types:
        if experiment not in peptide_experiments:
            not_existing_experiments.append(experiment)
    proteins_str = ", ".join(peptide_proteins)
    experiment_str = ", ".join(peptide_experiments)
    not_existing_experiment_str = ", ".join(not_existing_experiments)
    if len(peptide_proteins) == 1:
        unique_str = "True"
    else:
        unique_str = "False"
    return (
        f"Proteins: {proteins_str}\n"
        f"Unique: {unique_str}\n"
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

    fig, ax = plt.subplots(figsize=(1, 4))
    fig.subplots_adjust(right=0.5)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "", [x.hex for x in reversed(colors)]
    )
    # norm = mpl.colors.Normalize(
    #     vmin=amount[amount_list[-1]], vmax=amount[amount_list[0]]
    # )
    bounds = [amount[x] for x in reversed(amount_list)]
    norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)

    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)
    cb1.set_label("Amount of peptide")
    fig.savefig("Amount_of_peptide.png")

    return colors


@typing.no_type_check
def search_peptide_height(
    peptide: pd.Series,
    config,
) -> int:
    """
    Find peptide height in order to display it correctly.

    :param peptide: peptide information (pd.Series)
    :param config: configuration of html properties
    :return: height of peptide (top parameter in html)
    """
    max_height = (
        config["protein_line_height"] + SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES
    )
    x, y = peptide["Coords"]
    for key, value in config["peptides_already_displayed"].items():
        if key[0] >= y or key[1] > x:
            tmp = value + SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES
            if tmp > max_height:
                max_height = tmp
    config["peptides_already_displayed"][(x, y)] = max_height
    return max_height


@typing.no_type_check
def add_peptide_on_html_page(
    soup: bs4.BeautifulSoup,
    peptide: pd.Series,
    peptide_color,
    tooltip_text: str,
    is_transfer_element: bool,
    config,
) -> None:
    """
    Add peptide to div in order to display it on page.

    :param soup: html page handler
    :param peptide: peptide information (pd.Series)
    :param peptide_color: color of peptide
    :param tooltip_text: string containing peptide extra information
    :param is_transfer_element: if True left margin is set to 0
    :param config: configuration of html properties
    :return: None
    """
    peptide_height = search_peptide_height(peptide, config)
    if peptide_height == config["protein_line_height"]:
        peptide_height += SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES

    span_tag = soup.new_tag("span")
    peptide_sequence_str = peptide["Sequence"]
    if is_transfer_element:
        left_pos = 0
    else:
        left_pos = (peptide["Coords"][0]) % LINE_LENGTH_DISPLAYED
    unique_bool = UNIQUE_PATTERN.search(tooltip_text)
    if unique_bool.group(1) == "True":
        underline = f"{peptide_color} underline"
    else:
        underline = "none"
    span_tag["style"] = (
        f"color: {peptide_color};"
        f" left: {LEFT_TEXT_OFFSET + left_pos}ch;"
        f" top: {peptide_height}ch;"
        f" text-decoration: {underline};"
    )
    span_tag["title"] = tooltip_text
    config["div_center"].append(span_tag)

    if (
        peptide["Coords"][0]
        < config["current_first_index_of_protein_sequence"]
        and config["current_first_index_of_protein_sequence"]
        < peptide["Coords"][1]
    ):
        span_tag.append(
            soup.new_string(
                peptide["Sequence"][
                    : config["current_first_index_of_protein_sequence"]
                    - peptide["Coords"][1]
                ]
            )
        )
        new_coords_x = config["current_first_index_of_protein_sequence"]
        new_coords_y = peptide["Coords"][1]
        config["list_of_element_to_transfer"].append(
            (
                (new_coords_x, new_coords_y),
                peptide["Sequence"][
                    config["current_first_index_of_protein_sequence"]
                    - peptide["Coords"][1] :
                ],
                peptide_color,
                tooltip_text,
            )
        )
    else:
        span_tag.append(soup.new_string(peptide_sequence_str))


@typing.no_type_check
def add_next_line_of_protein_sequence(
    soup: bs4.BeautifulSoup,
    ref_seq_fasta: str,
    config,
) -> None:
    """
    Add next line contain protein sequence.

    :param soup: html page handler
    :param ref_seq_fasta: protein sequence
    :param config: configuration of html properties
    :return: None
    """
    # TODO: possible bug
    current_max_height = max(
        config["peptides_already_displayed"].values(), default=0
    )
    if current_max_height != 0:
        current_max_height += SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES

    span_tag = soup.new_tag("span")
    span_tag["style"] = (
        f"color: white;"
        f" left: {LEFT_TEXT_OFFSET}ch;"
        f" top: {current_max_height}ch"
    )
    span_tag.append(
        soup.new_string(
            ref_seq_fasta[
                config["current_first_index_of_protein_sequence"] : config[
                    "current_first_index_of_protein_sequence"
                ]
                + LINE_LENGTH_DISPLAYED
            ]
        )
    )
    config["div_center"].append(span_tag)
    config["current_first_index_of_protein_sequence"] += LINE_LENGTH_DISPLAYED
    config["protein_line_height"] = current_max_height


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
    display_config = {
        "current_first_index_of_protein_sequence": 0,
        "protein_line_height": 0,
        "peptides_already_displayed": {},
        "list_of_element_to_transfer": [],
        "div_center": soup.find("div", {"class": "center"}),
    }
    add_next_line_of_protein_sequence(soup, ref_seq_fasta, display_config)
    current_seq = ""
    for peptide in result:
        # TODO: nazwa do poprawy current first...
        if (
            peptide["Coords"][0]
            >= display_config["current_first_index_of_protein_sequence"]
        ):
            add_next_line_of_protein_sequence(
                soup, ref_seq_fasta, display_config
            )

            for element_to_transfer in display_config[
                "list_of_element_to_transfer"
            ]:
                transfer_peptide_info = {
                    "Sequence": element_to_transfer[1],
                    "Coords": element_to_transfer[0],
                }
                transfer_peptide_info_series = pd.Series(transfer_peptide_info)
                add_peptide_on_html_page(
                    soup,
                    transfer_peptide_info_series,
                    element_to_transfer[2],
                    element_to_transfer[3],
                    True,
                    display_config,
                )

            display_config["list_of_element_to_transfer"] = []

        if peptide["Sequence"] == current_seq:
            continue

        # sprawdzic gdy nie ma w linijce zadnego peptydu
        while (
            peptide["Coords"][0]
            > display_config["current_first_index_of_protein_sequence"]
        ):
            add_next_line_of_protein_sequence(
                soup, ref_seq_fasta, display_config
            )

        current_seq = peptide["Sequence"]
        tooltip_text = create_table_of_information_about_peptide(
            result, current_seq, amount, all_possible_experiment_type
        )
        peptide_color = colors[amount_list.index(f"{current_seq}")].hex
        add_peptide_on_html_page(
            soup,
            peptide,
            peptide_color,
            tooltip_text,
            False,
            display_config,
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

    :return: None
    """
    ref_seq_fasta = read_ref_seq_fasta()
    result_of_experiment = pd.read_excel(
        r"..\pliki-SMARCA5\SMARCA5.xlsx", sheet_name="SMARCA5"
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
