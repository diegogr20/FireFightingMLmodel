# Cell constant parameters
IGNTIME    = 30.0
TM         = 10.0

class Cell:
    # Maximum number of particle
    NPARTMAX   = 1

    def __init__(self):
        # Cell's parameters
        self.TREEAREA  = 1
        self.CLOCK     = 0.0
        self.QMAXTR    = 10.0
        self.TIGNTR    = IGNTIME
        self.TMEM      = TM
        self.BURNSTAT  = 0.0
        self.BURNSTAT2 = 1.0
        self.NPARTTR   = Cell.NPARTMAX


    def cannot_burn(self):
        """ Define a cell that cannot burn."""
        self.QMAXTR  = 0.0
        self.QMAXBLD = 0.0

    def tree_cell(self):
        """Define the cells with tree areas"""
        self.TIGNTR = IGNTIME*6
        self.TMEM = TM*2

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        raise KeyError(f"Cell doesn't have key '{key}'")

    def __repr__(self):
        return f"Cell(burning_state={self.BURNSTAT:.2f})"
