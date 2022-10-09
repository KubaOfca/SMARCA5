"""Main module."""
import ttkbootstrap as tk  # type: ignore
import re
from App import App


LINE_LENGTH_DISPLAYED = 80
SIZE_OF_THE_TOP_SPACE_BETWEEN_THE_LINES = 3
LEFT_TEXT_OFFSET = 3
UNIQUE_PATTERN = re.compile(r".*Unique: (.*)\n.*")


# flake8: noqa E203


# TODO: Better name convention of .html files
# TODO: function to handle operations from 86line
# TODO: GUI as class
# TODO: too many figure open problem


if __name__ == '__main__':
    root = tk.Window(themename="darkly")
    app = App(root)
    root.mainloop()
