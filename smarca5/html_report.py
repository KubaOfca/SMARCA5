from bs4 import BeautifulSoup
from main import load_file

SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES = 3


class HtmlReport:
    def __int__(self, interline, text_left_offset, line_length):
        self.soup = soup_obj_generator()
        self.interline = interline
        self.text_left_offset = text_left_offset
        self.line_length = line_length
        self.protein_seq_color = ""


class ReportPage(HtmlReport):
    pass


def soup_obj_generator():
    with open(load_file(r"html_files\index.html"), "r") as html:
        soup = BeautifulSoup(html, "html.parser")
    return soup
