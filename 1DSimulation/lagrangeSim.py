import numpy as np
from particleSim import Particle
from cellSim import Cell

STTHR      = 0.2
RELT       = 0.2
LANGFACTOR = 0.15

A         = 0.3
TLENGTH   = 50.0


def normal(mu=0, sigma=1):
    return np.random.normal(mu, sigma)

class Lagrange:
    def __init__(self, grid, U=0):
        # self.grid_copy = deepcopy(grid)
        self.grid       = grid
        self.UXM        = U
        self.UPRIME     = U*A
        self.TFREQ      = self.UPRIME/TLENGTH
        self.CLANG      = np.sqrt(2.0*self.UPRIME**3 / TLENGTH)
        self.ITYPESPARK = None
        self.DT = 2 * (self.grid.DX / self.UXM)
        self.SQRTDT = np.sqrt(self.DT)

        self.NPART = 0
        self.particles = None # particles array
        self._init_particles()

    def _init_particles(self):

        self.NPART = self.grid.get_npart()
        self.TERRAIN = self.grid.get_cells_prop('TMEM')
        self.particles = np.empty(self.NPART, dtype=object)
        for i in range(self.NPART):
            tmem = self.TERRAIN[i]
            ux = self.UXM + self.UPRIME*normal()*LANGFACTOR
            self.particles[i] = Particle(tmem=tmem, ux=ux)
            
        NX = self.grid.NX
        icounter = 0
        for i in range(NX - 1):
            x    = self.grid.XCELL[i]
            self.particles[icounter].update(x=x)
            icounter += 1


    def start_fire(self):
        
        self.grid.CELLS[0].BURNSTAT = 1
        for ip in range(self.NPART):
            particle = self.particles[ip]
            X = self.grid.XCELL[0]
            x = particle.x
            cell = self.grid.CELLS[0]
            cell.BURNSTAT2 = 0
            cell.BURNSTAT  = 1
            if (abs(X - x) < self.grid.DX):
                particle.update(state=1, factor=LANGFACTOR)

    def _move_particle(self, ip, xlarge):
        
        particle = self.particles[ip]
        props = ["state", "x", "ux", "factor", "tmem"]
        state, x, ux, factor, tmem = particle.get_from_keys(props)
        xlarge = xlarge
        if state > STTHR:
            DU  = 0 -(ux - self.UXM)*2.0*self.TFREQ*self.DT + self.CLANG*self.SQRTDT*normal()
            UXP = ux + DU
            XP  = x + UXP*self.DT*factor
            STP = state*np.exp(-self.DT/tmem)
            particle.update(ux=UXP, x=XP, state=STP)
            if XP > xlarge:
                xlarge = XP
        return xlarge

            
    def _ignite_cells(self, istep, ip):

        particle = self.particles[ip]
        state, x = particle.get_from_keys(["state", "x"])
        if state > STTHR:
            for i in range(self.grid.NX-1):
                if abs(x - self.grid.XCELL[i]) <= self.grid.DX/2:
                    INDX = i
            cell = self.grid.CELLS[INDX]
            if (cell.QMAXTR > 0) and cell.BURNSTAT == 0 and cell.BURNSTAT2 == 1:
                cell.BURNSTAT = 1
                cell.CLOCK    = istep

    def _launch_particles(self, istep):

        for i in range(self.grid.NX-1):
            cell = self.grid.CELLS[i]
            TLOCAL = istep - cell.CLOCK
            TCRIT  = cell.TIGNTR * (1 + RELT*normal())
            if cell.BURNSTAT == 1 and TLOCAL > TCRIT and cell.BURNSTAT2 == 1:
                LOCALF = LANGFACTOR
                indp = (i)*Cell.NPARTMAX
                cell.CLOCK = istep
                for k in range(cell.NPARTTR):
                    self.particles[k + indp].update(state=1.0, factor=LOCALF)
                cell.BURNSTAT2 = 0

    def propagate(self, istep):
        xlarge = 0
        for ip in range(self.NPART):
            xlarge = self._move_particle(ip,xlarge)
            self._ignite_cells(istep, ip)
        self._launch_particles(istep)
        return xlarge

    def get_particles_props(self, *props, array=None):

        n = len(props)
        if not isinstance(array, np.ndarray):
            array = self.particles
        if n == 0:
            return []
        elif n == 1:
            key = props[0]
            properties = np.empty(len(array))
        else:
            properties = np.empty((n, len(array)))

        for ip, particle in enumerate(array):
            if n == 1:
                properties[ip] = particle[key]
            else:
                for i, key in enumerate(props):
                    properties[i, ip] = particle[key]
        return properties


    def launch(self):

        istep = self.DT
        loop = 0
        while loop == 0:
            xbig = self.propagate(istep)
            istep += self.DT
            if (abs(xbig-self.grid.XCELL[len(self.grid.XCELL) - 5])<= self.grid.DX/2):
                loop = 1
            ST = self.get_particles_props('state')
            BRN = self.grid.get_cells_prop('BURNSTAT')
            BRN2 = self.grid.get_cells_prop('BURNSTAT2')
            if not((BRN*BRN2 == 1).any()) and (ST < STTHR).all():
                loop = 1
                istep = 'inf'
        return istep

