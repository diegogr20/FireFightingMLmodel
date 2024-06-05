import numpy as np
from random import random
from cellSimN import Cell

class Grid:

    def __init__(self, DX=10, XMIN=0, l1=0, l2=0, l3=0, T1=0):
        self.DX = DX
        self.XMIN = XMIN
        self.L1 = int(l1/DX)
        self.L2 = int(l2/DX)
        self.L3 = int(l3/DX)
        self.T1 = T1
        self.LB = 5
        self.NX = self.L1 + self.L2 + self.L3 + self.LB + 1
        # Cell centers
        self.XCELL = np.empty(self.NX-1)
        # All Cells
        self.CELLS = np.empty((self.NX-1), dtype=object)
        for i in range(self.NX-1):
            self.XCELL[i] = XMIN + DX*(i + 0.5)
            self.CELLS[i] = Cell()
            if ((i < self.L1 or (i >= (self.L1 + self.L2) and i < (self.L1 + self.L2 + self.L3))) and self.T1 == 1):
                self.CELLS[i].tree_cell()
            elif i >= (self.L1 + self.L2 + self.L3) or (i >= self.L1 and i < (self.L1 + self.L2)):
                self.CELLS[i].cannot_burn()
                

    def get_cells_prop(self, prop_name):
        DATA = np.empty(self.NX - 1)
        for i in range(self.NX - 1):
            DATA[i] = self.CELLS[i][prop_name]
        return DATA

    def get_heat_released(self):
        return self.get_cells_prop("QMAXTR")

    def get_burn_state(self):
        return self.get_cells_prop("BURNSTAT")

    def get_npart(self):
        return Cell.NPARTMAX*(self.NX - 1)

    def __repr__(self):
        return f"Grid(NX={self.NX})"
