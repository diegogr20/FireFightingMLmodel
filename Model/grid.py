import numpy as np
from random import random
from cell import Cell

class Grid:
    """Object to model the domain grid.

    Params:
    -------

    IN REAL APPLICATION, this comes from outside UNITS: m
    NX     (int, optional): # of grid lines (NX-1) is the number of cells. Defaults to 101
    NY     (int, optional): # of grid columns (NY-1) is the number of cells. Defaults to 101
    XMIN (float, optional): lower bound of grid in x-dir. Defaults to 0
    XMAX (float, optional): upper bound of grid in x-dir. Defaults to 2000
    YMIN (float, optional): lower bound of grid in y-dir. Defaults to 0
    YMAX (float, optional): upper bound of grid in y-dir. Defaults to 2000
    """

    def __init__(self, DX=10, DY=10, XMIN=0, XMAX=1000, YMIN=0, YMAX=1000):
        self.DX = DX
        self.DY = DY
        self.XMIN = XMIN
        self.YMIN = YMIN
        self.XMAX = XMAX
        self.YMAX = YMAX
        self.NX = int(((XMAX - XMIN) / (DX)) + 1)
        self.NY = int(((YMAX - YMIN) / (DY)) + 1)
        self.HIx = np.arange(76,94,1)
        self.HIy = np.arange(6,42,1)
        self.HITx = np.arange(74,96,1)
        self.HITy = np.arange(4,44,1)
        q, t = np.meshgrid(self.HITx,self.HITy)
        self.HITT = np.array([q, t]).reshape([2,len(self.HITx)*len(self.HITy)]).T
        # Cell centers
        self.XCELL = np.empty((self.NX-1, self.NY-1))
        self.YCELL = np.empty((self.NX-1, self.NY-1))
        # All Cells
        self.CELLS = np.empty((self.NX-1, self.NY-1), dtype=object)
        for i in range(self.NX-1):
            for j in range(self.NY-1):
                self.XCELL[i, j] = XMIN + self.DX*(i + 0.5)
                self.YCELL[i, j] = YMIN + self.DY*(j + 0.5)
                self.CELLS[i, j] = Cell()
        
        self.HIT = []
        for k in range(len(self.HITT)):
            self.HIT = self.HIT + [self.XCELL[self.HITT[k,0],self.HITT[k,1]],self.YCELL[self.HITT[k,0],self.HITT[k,1]]]
        self.HIT = np.array(self.HIT).reshape(len(self.HITT),2)

    def set_layout(self):
        for i in self.HIx:
            for j in self.HIy:
                self.CELLS[i,j].cannot_burn()

        #Tree area 1
        x1 = np.arange(0,8,1)
        y1 = np.arange(0,52,1)
        for i in x1:
            for j in y1:
                self.CELLS[i,j].tree_cell()
        
        #Tree area 2
        x2 = np.arange(32,100,1)
        y2 = np.arange(62,100,1)
        for i in x2:
            for j in y2:
                self.CELLS[i,j].tree_cell()
        
        #Tree area 3
        x3 = np.arange(8,32,1)
        y13 = []
        y23 = []
        for i in x3:
            y13 = y13 + [int(1.260869565*i + 21.91304348)]
            y23 = y23 + [int(2.086956522*i + 34.30437483)]
        
        for i in range(len(x3)):
            yt = np.arange(y13[i], y23[i] + 1, 1)
            for j in yt:
                self.CELLS[x3[i],j].tree_cell()


    def get_cells_prop(self, prop_name):
        """Get an array of `prop_name` values for each cell of the grid."""
        DATA = np.empty((self.NX - 1, self.NY - 1))
        for i in range(self.NX - 1):
            for j in range(self.NY - 1):
                DATA[i, j] = self.CELLS[i, j][prop_name]
        return DATA

    def get_heat_released(self):
        return self.get_cells_prop("QMAXTR")

    def get_burn_state(self):
        return self.get_cells_prop("BURNSTAT")

    def get_npart(self):
        """Return the number of particle."""
        return Cell.NPARTMAX*(self.NX - 1)*(self.NY - 1)

    def __repr__(self):
        return f"Grid(NX={self.NX}, NY={self.NY})"
