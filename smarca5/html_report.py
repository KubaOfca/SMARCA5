class HtmlReport:
    def __init__(
        self, soup, groups, samples, text_left_offset, interline, line_length
    ):
        self.base_soup = soup
        self.interline = interline
        self.text_left_offset = text_left_offset
        self.line_length = line_length
        self.protein_seq_color = ""
        self.groups = groups
        self.samples = samples
        self.fill_right_menus(self.groups, "Groups")
        self.fill_right_menus(self.samples, "Samples")

    def fill_right_menus(self, menu_elements, menu_type) -> None:
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
            if menu_type == "Groups":
                a_tag = self.base_soup.new_tag(
                    "a", href=f"type_group-{element_name}.html"
                )
            else:
                a_tag = self.base_soup.new_tag(
                    "a", href=f"type_sample-{element_name}.html"
                )
            a_tag.string = element_name
            li_tag.append(a_tag)
            ul_tag.append(li_tag)
