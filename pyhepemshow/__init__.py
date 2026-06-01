import importlib.resources

from . import forward, forward_ad, reverse


def get_hepemshow_data_file():
    resource = importlib.resources.files("pyhepemshow").joinpath("data/hepem_data.json")
    with importlib.resources.as_file(resource) as path:
        return str(path)

