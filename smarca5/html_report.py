from bs4 import BeautifulSoup
from main import load_file
from copy import deepcopy
from math import ceil

class HtmlReport:
    def __int__(self, groups, samples, interline, text_left_offset, line_length):
        self.base_soup = soup_obj_generator()
        self.interline = interline
        self.text_left_offset = text_left_offset
        self.line_length = line_length
        self.protein_seq_color = ""
        self.groups = groups
        self.samples = samples
        self.fill_right_menus(self.groups, "Groups")
        self.fill_right_menus(self.samples, "Samples")

    def fill_right_menus(self, menu_elements, menu_type
                         ) -> None:
        """
        Create right panel with links to other page result.

        :param soup: html page handler
        :param experiment: list of experiments
        :param current: actual experiment name page
        :return:
        """
        ul_tag = self.base_soup.find("ul", {"id": f"{menu_type}"})
        for element_name in menu_elements:
            li_tag = self.base_soup.new_tag("li")
            a_tag = self.base_soup.new_tag("a", href=f"{element_name}.html")
            a_tag.string = element_name
            li_tag.append(a_tag)
            ul_tag.append(li_tag)


class ReportPage(HtmlReport):
    def __init__(self, experiment, protein_seq, peptides):
        super().__init__()  # may cause problems
        self.soup = deepcopy(self.base_soup)
        self.number_of_amino_acids_already_display = 0
        self.first_empty_line_height = 0 # first empty line under already created alignment (some kind of temporary EOF)
        self.peptides_already_displayed = {}
        self.list_of_element_to_transfer = [] # ?
        self.div_center = self.soup.find("div", {"class": "center"})
        self.soup.title.string = experiment
        self.protein_seq = protein_seq
        self.peptides = peptides
    #def underline correct page in menu
    #def generate correct colorbar
    #def generate contex menu
    #def generate aligment

    def fill_alignment_window(self):
        number_of_protein_lines = ceil(len(self.protein_seq) / self.line_length)
        for _ in range(number_of_protein_lines):
            self.add_next_line_of_alignment()
            pass

    def add_next_protein_line(self):
        # add protein
        span_tag = self.soup.new_tag("span")
        span_tag["style"] = (
            f"color: #EEEEEE;"
            f" left: {self.text_left_offset}ch;"
            f" top: {self.first_empty_line_height}ch"
        )
        span_tag.append(self.soup.new_string(
            self.protein_seq[self.number_of_amino_acids_already_display:
                             self.number_of_amino_acids_already_display + self.line_length]))
        self.number_of_amino_acids_already_display += self.line_length
        self.div_center.append(span_tag)
        self.first_empty_line_height += self.interline

    def add_next_peptides_line(self):
        pass

    def add_next_line_of_alignment(self):
        self.add_next_protein_line()
        # add peptides

        # below add peptides function test
        #TODO: search peptide heigh logic
        for peptide in self.peptides:
            if peptide["Coords"][0] >= self.number_of_amino_acids_already_display:
                for element_to_transfer in display_config["list_of_element_to_transfer"]:
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
                add_next_line_of_protein_sequence(soup, ref_seq_fasta, display_config)

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

        while (
                len(ref_seq_fasta) > display_config["current_first_index_of_protein_sequence"]
        ):
            add_next_line_of_protein_sequence(soup, ref_seq_fasta, display_config)
        # save changes to page.html and display page with result
        with open(load_file(rf"html_files\{exp}.html"), "w") as page:
            page.write(soup.prettify(formatter="html"))  # soup.encode(formatter="html")
            if is_first_report:
                webbrowser.open_new_tab(load_file(rf"html_files\{exp}.html"))

def soup_obj_generator():
    with open(load_file(r"html_files\index.html"), "r") as html:
        soup = BeautifulSoup(html, "html.parser")
    return soup
