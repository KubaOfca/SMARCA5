"""Main module."""
import ttkbootstrap as tk  # type: ignore
from App import App


if __name__ == "__main__":
    root = tk.Window(themename="darkly")
    app = App(root)
    root.mainloop()
