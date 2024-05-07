import numpy as np
from yee3d.Grid import Grid

"""
Source class: all sources inherit from this class
Grid
update
"""
class Source:
    def __init__(self, name:str = None) -> None:
        self.Grid = None
        self.name = name
        pass

    def bindToGrid(self, g: Grid):
        pass

    def update(self):
        pass
