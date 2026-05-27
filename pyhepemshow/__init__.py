import importlib.resources
import os

def get_hepemshow_data_file():
    pkg_data = os.path.join(os.path.dirname(__file__), "data", "hepem_data.json")
    return str(pkg_data)
