from copy import deepcopy
from math import ceil
import pandas as pd
from colour import Color
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl  # type: ignore
import re

UNIQUE_PATTERN = re.compile(r".*Unique: (.*)\n.*")
UNDERLINE_MENU_STYLE = "text-decoration: underline; text-underline-offset: 5px; text-decoration-color:red"


class ReportPage:
    def __init__(
        self,
        html_report_config,
        type_of_experiment,
        experiment,
        peptides_from_given_experiment,
        alignment_obj,
    ):
        self.soup = deepcopy(html_report_config.base_soup)
        self.number_of_amino_acids_already_display = 0
        self.first_empty_line_height = (
            0  # first empty line under already created alignment
        )
        self.peptides_already_displayed_info = {}
        self.df_of_peptides_to_transfer = pd.DataFrame(
            columns=[
                "Sequence",
                "FullSequence",
                "Proteins",
                "Experiment",
                "Start",
                "End",
            ]
        )
        self.div_center = self.soup.find("div", {"class": "center"})
        self.soup.title.string = experiment
        self.type_of_experiment = type_of_experiment
        self.peptides_from_given_experiment = peptides_from_given_experiment
        self.alignment_obj = alignment_obj
        self.color_bar_fig = None
        self.html_report_config = html_report_config
        self.amount_of_peptides = self.count_amount_of_the_same_peptides()
        self.color_mapped_with_amount = self.create_color_bar_image()
        self.add_underline_to_menu_element_of_current_experiment()

    def add_underline_to_menu_element_of_current_experiment(self):
        a_tag = self.soup.find(
            "a",
            {
                "href": {
                    f"type_{self.type_of_experiment}-{self.soup.title.string}.html"
                }
            },
        )
        a_tag["style"] = UNDERLINE_MENU_STYLE

    def fill_alignment_window(self):
        number_of_protein_lines = ceil(
            len(self.alignment_obj.protein_seq)
            / self.html_report_config.line_length
        )
        for _ in range(number_of_protein_lines):
            self.add_next_segment_of_alignment()

    def add_next_segment_of_alignment(self):
        self.add_next_protein_line()
        self.add_next_peptides_line()

    def add_next_protein_line(self):
        span_tag = self.soup.new_tag("span")
        span_tag["style"] = (
            f"color: #EEEEEE;"
            f" left: {self.html_report_config.text_left_offset}ch;"
            f" top: {self.first_empty_line_height}px"
        )
        span_tag.append(
            self.soup.new_string(
                self.alignment_obj.protein_seq[
                    self.number_of_amino_acids_already_display : self.number_of_amino_acids_already_display
                    + self.html_report_config.line_length
                ]
            )
        )
        self.number_of_amino_acids_already_display += (
            self.html_report_config.line_length
        )
        self.first_empty_line_height += self.html_report_config.interline
        self.div_center.append(span_tag)

    def add_next_peptides_line(self):
        # there should be a function
        # TODO: list of element to transfer should look the same as peptide
        for (
            index,
            peptide_to_transfer,
        ) in self.df_of_peptides_to_transfer.iterrows():
            start = peptide_to_transfer["Start"]
            end = peptide_to_transfer["End"]
            current_peptide_height = self.search_peptide_height(
                start, end, peptide_to_transfer["FullSequence"]
            )
            # prepare info ect. about peptide
            peptide_span_tag = self.create_span_tag_of_peptide(
                peptide_to_transfer, current_peptide_height, True
            )
            self.div_center.append(peptide_span_tag)

        self.df_of_peptides_to_transfer = pd.DataFrame(
            columns=[
                "Sequence",
                "FullSequence",
                "Proteins",
                "Experiment",
                "Start",
                "End",
            ]
        )

        peptides_from_a_given_segment = (
            self.peptides_from_given_experiment.loc[
                (
                    self.peptides_from_given_experiment["Start"]
                    >= self.number_of_amino_acids_already_display
                    - self.html_report_config.line_length
                )
                & (
                    self.peptides_from_given_experiment["Start"]
                    < self.number_of_amino_acids_already_display
                )
            ]
        )

        for index, peptide in peptides_from_a_given_segment.iterrows():
            if (
                peptide["Sequence"] not in self.peptides_already_displayed_info
                and peptide["Start"]
                < self.number_of_amino_acids_already_display
            ):
                start = peptide["Start"]
                end = peptide["End"]
                current_peptide_height = self.search_peptide_height(
                    start, end, peptide["Sequence"]
                )
                # prepare info ect. about peptide
                peptide_span_tag = self.create_span_tag_of_peptide(
                    peptide, current_peptide_height
                )
                self.div_center.append(peptide_span_tag)
        self.first_empty_line_height = (
            max(
                [
                    height
                    for height, coords in self.peptides_already_displayed_info.values()
                ],
                default=self.first_empty_line_height,
            )
            + self.html_report_config.interline
        )

        self.peptides_already_displayed_info = {}

    def search_peptide_height(self, start, end, peptide_seq):
        overlay_peptides_height = []
        for height, coords in self.peptides_already_displayed_info.values():
            start_displayed, end_displayed = coords
            if (start_displayed <= start < end_displayed) or (start_displayed <= end < end_displayed):
                overlay_peptides_height.append(height)

        if overlay_peptides_height:
            current_peptide_height = (
                max(overlay_peptides_height)
                + self.html_report_config.interline
            )
        else:
            current_peptide_height = self.first_empty_line_height

        if end >= self.number_of_amino_acids_already_display:
            end = self.number_of_amino_acids_already_display
        self.peptides_already_displayed_info[peptide_seq] = [
            current_peptide_height,
            (start, end),
        ]
        return current_peptide_height

    def create_span_tag_of_peptide(
        self, peptide, current_peptide_height, is_transfer=False
    ):
        span_tag = self.soup.new_tag("span")
        start_position_in_line_of_peptide = (
            peptide["Start"] % self.html_report_config.line_length
        )
        if not is_transfer:
            peptide_color = self.color_mapped_with_amount[
                self.amount_of_peptides[peptide["Sequence"]]
            ]
        else:
            peptide_color = self.color_mapped_with_amount[
                self.amount_of_peptides[peptide["FullSequence"]]
            ]
        tooltip_text = self.create_table_of_information_about_peptide(
            peptide, is_transfer
        )
        span_tag["title"] = tooltip_text
        unique_bool = UNIQUE_PATTERN.search(tooltip_text)
        if unique_bool.group(1) == "True":
            underline = f"{peptide_color} underline"
        else:
            underline = "none"

        span_tag["style"] = (
            f"color: {peptide_color};"
            f" left: {start_position_in_line_of_peptide + self.html_report_config.text_left_offset}ch;"
            f" top: {current_peptide_height}px;"
            f" text-decoration: {underline};"
        )
        # seq of peptide displayed if is longer than already displayed protein line
        if (
            peptide["Start"]
            < self.number_of_amino_acids_already_display
            < peptide["End"]
        ):
            peptide_sequence_before_transfer = peptide["Sequence"][
                : self.number_of_amino_acids_already_display - peptide["Start"]
            ]
            span_tag.append(
                self.soup.new_string(peptide_sequence_before_transfer)
            )

            # transfer
            peptide_sequence_after_transfer = peptide["Sequence"][
                self.number_of_amino_acids_already_display - peptide["Start"] :
            ]
            start_coord_after_transfer = (
                self.number_of_amino_acids_already_display
            )

            peptide_to_transfer_data = pd.DataFrame(
                {
                    "Sequence": [peptide_sequence_after_transfer],
                    "FullSequence": [peptide["Sequence"]],
                    "Proteins": [peptide["Proteins"]],
                    "Experiment": [peptide["Experiment"]],
                    "Start": [start_coord_after_transfer],
                    "End": [peptide["End"]],
                }
            )
            self.df_of_peptides_to_transfer = pd.concat(
                [self.df_of_peptides_to_transfer, peptide_to_transfer_data]
            )
        else:
            span_tag.append(self.soup.new_string(peptide["Sequence"]))

        return span_tag

    def create_table_of_information_about_peptide(
        self, current_peptide, is_transfer
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
        if is_transfer:
            current_peptide_seq = current_peptide["FullSequence"]
        else:
            current_peptide_seq = current_peptide["Sequence"]
        proteins = current_peptide["Proteins"].split(";")
        if len(proteins) == 1:
            unique_str = "True"
        else:
            unique_str = "False"

        peptide_experiments = []
        for index, peptide in self.peptides_from_given_experiment.iterrows():
            if current_peptide_seq == peptide["Sequence"]:
                if peptide["Experiment"] not in peptide_experiments:
                    peptide_experiments.append(peptide["Experiment"])
        # calc how many the same peptides is in this dataset
        amount = self.amount_of_peptides[current_peptide_seq]
        proteins_str = ", ".join(proteins)
        peptide_experiments_str = ", ".join(peptide_experiments)

        return (
            f"Sequence: {current_peptide_seq}\n"
            f"Proteins: {proteins_str}\n"
            f"Unique: {unique_str}\n"
            f"Experiment: {peptide_experiments_str}\n"
            f"Amount: {amount}"
        )

    def create_color_bar_image(self):
        """
        Create color bar image.

        :param amount: dictionary where key is sequence and value is amount
        :param amount_list: peptide sequence list in descending order
        :return: list of colors used in color bar
        """
        unique_amount_of_peptides = sorted(
            list({value for value in self.amount_of_peptides.values()})
        )  # [1, 5, 6, 10...]
        color_map_with_amount = {}
        blue = Color("LightBlue")
        colors = list(
            blue.range_to(
                value=Color("MediumBlue"), steps=len(unique_amount_of_peptides)
            )
        )
        if len(colors) == 1:
            colors.append(Color("MediumBlue"))
            unique_amount_of_peptides.append(2)
        # blue - min, green - max
        for index, amount in enumerate(unique_amount_of_peptides):
            color_map_with_amount[amount] = colors[index].hex

        fig, ax = plt.subplots(figsize=(1, 4))
        fig.subplots_adjust(right=0.5)

        c_map = mpl.colors.LinearSegmentedColormap.from_list(
            "", [x.hex for x in colors]
        )

        bounds = [x for x in unique_amount_of_peptides]
        bounds.insert(0, 0)
        norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)

        cb1 = mpl.colorbar.ColorbarBase(ax=ax, cmap=c_map, norm=norm)
        cb1.set_label("Amount of peptide")
        self.color_bar_fig = fig
        return color_map_with_amount

    def count_amount_of_the_same_peptides(self):
        amount_of_peptides = dict(
            self.peptides_from_given_experiment["Sequence"].value_counts()
        )
        return amount_of_peptides
