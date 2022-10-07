from bs4 import BeautifulSoup
from main import load_file

SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES = 3
class HtmlRaport:
    def __int__(self):
        self.soup = soup_obj_generator()
        title = ""
        display_config = {
            "current_first_index_of_protein_sequence": 0,
            "protein_line_height": -SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES,
            "peptides_already_displayed": {},
            "list_of_element_to_transfer": [],
            "div_center": self.soup.find("div", {"class": "center"}),
        }

class RaportPage
def soup_obj_generator():
    with open(load_file(r"html_files\index.html"), "r") as html:
        soup = BeautifulSoup(html, "html.parser")
    return soup
