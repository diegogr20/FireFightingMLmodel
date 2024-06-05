import numpy as np
import matplotlib.pyplot as plt
from particle import Particle
from cell import Cell
from Extinguish import Airplane
from mlmodel import MLmodel

# Message to show when the simulation is launched

MESSAGE = """
==========================================================================================
Lagrangian stochastic model for fire propagation in forests and the forest/urban interface
==========================================================================================

Model parameters :
{grid}, NSTEPS = {nsteps}
"""

# INITIALISATIONS
#
# STTHR =  Threshold to decide if particle is "on". THIS COULD BE USED TO
# MODEL QUENCHING. OR, if made f(space) to model ignitability of the region.
#
# LANGFACTOR = Maximum factor (0-1) of the air Langevin walk step that the flame
# particles experience.
# RELT = randomize a bit the ignition time criterion so the propagation is
# not too jumpy
STTHR      = 0.2
RELT       = 0.2
LANGFACTOR = 0.15

# Creating array for extinguishing actions to happen
STEPS     = 600
exting = np.zeros(9)
for i in range(9):
    exting[i] = int((0.1 + 0.1*i)*STEPS)


# WIND SPEED AND LANGEVIN PARAMETERS
# All in SI units
# UX = mean wind speed in X-dir (later, it can be f(Z))
# UY = mean wind speed in Y-dir
# A = relative turbulence intensity (u'/<U>) - will be adjusted for proper
#     Langevin walk later
# FIRECONST = factor to model slower fire speed, kept 1 for now
# The rest are factors used for modelling the fire.

UX        = 12.0
VY        = -13.0
FIRECONST = 1.0
UXM       = UX*FIRECONST
VYM       = VY*FIRECONST
USPEED    = np.sqrt(UXM**2 + VYM**2)
A         = 0.3
UPRIME    = USPEED*A
TLENGTH   = 50.0
TFREQ     = UPRIME / TLENGTH
CLANG     = np.sqrt(2.0*UPRIME**3 / TLENGTH)


def normal(mu=0, sigma=1):
    return np.random.normal(mu, sigma)

class Lagrange:
    """Implement the lagragian model for fire propagation in forests and the forest/urban interface

    Params:
    -------
    grid       (grid.Grid): grid representing the surface area
    NSTEPS (int, optional): number of time to propagate the model to/
                            number of timesteps of the simulation


    Attributes:
    -----------
    DT       (float): timestep in (s)
    SQRDT    (float): square root of DT
    NPART      (int): number of particles in the grid
    ITYPESPARK (int): number to specify the type of the initial spark
                      1 for fire line, 0 for a point

    TIME  (np.ndarray[float]): the time axis for the whole simulation
    BURNX (np.ndarray[float]): burning progress averaged in y-dir as a function of time & x-dir


    Methods:
    --------
    start_fire(ITYPESPARK): start the fire by setting either a front line
                            or a point fire.
    propagate(istep): propagate the model for the given timestep.


    """
    def __init__(self, grid, NSTEPS=STEPS):
        # self.grid_copy = deepcopy(grid)
        self.grid       = grid
        self.NSTEPS     = NSTEPS
        self.ITYPESPARK = None
        self.DT = 2 * (self.grid.DX / USPEED)
        self.SQRTDT = np.sqrt(self.DT)
        self.TIME = np.empty(self.NSTEPS)
        for i in range(self.NSTEPS):
            self.TIME[i] = (i + 1) * self.DT

        self.BURNX = np.empty((NSTEPS, self.grid.NX-1))
        # SAVE BURNSTAT AT EVERY TIMESTEP
        self.BURNEVOL = np.empty((NSTEPS, self.grid.NX-1, self.grid.NY-1))

        self.NPART = 0
        self.particles = None # particles array
        self._init_particles()

    def _init_particles(self):
        """Initialisation of particle's global state.

            Use large memory for the time being (all particles stay lit for a long
            time)
            Give all particles some random velocity component, so the drift in the
            later Langevin walk makes sense.
        """
        self.NPART = self.grid.get_npart()
        self.TERR = self.grid.get_cells_prop('TMEM')
        self.TERRAIN = self.TERR.reshape(self.NPART)
        self.particles = np.empty(self.NPART, dtype=object)
        for i in range(self.NPART):
            tmem = self.TERRAIN[i]
            ux = UXM + UPRIME*normal()*LANGFACTOR
            vy = VYM + UPRIME*normal()*LANGFACTOR
            self.particles[i] = Particle(tmem=tmem, ux=ux, vy=vy)
        #
        # PUT THE PARTICLES IN THE CELLS.
        # LOOP OVER CELLS AND DEFINE THEIR PARTICLES.
        # FOR NOW, ONLY POSITION DEPENDS ON SPACE HEIGHT & MEMORY DO NOT.
        # FIRST THE TREE PARTICLES, THEN THE BUILDING PARTICLES.
        #
        NX = self.grid.NX
        NY = self.grid.NY
        icounter = 0
        for i in range(NX - 1):
            for j in range(NY - 1):
                cell = self.grid.CELLS[i, j]
                x    = self.grid.XCELL[i, j]
                y    = self.grid.YCELL[i, j]
                for k in range(cell.NPARTMAX):
                    self.particles[k + icounter].update(x=x, y=y)
                icounter += cell.NPARTMAX


    def start_fire(self, ITYPESPARK=1):
        """Start the fire

            ITYPESPARK = 0 (point)  1 (line)
        """
        self.ITYPESPARK = ITYPESPARK
        if ITYPESPARK == 0:
            IF = 0
            IT = []
            JRANGE = [49]
        elif ITYPESPARK == 1:
            IF = 1
            IT = np.arange(1,30,1)
            JRANGE = np.arange(48,99,1)
        for JF in JRANGE:
            self.grid.CELLS[IF, JF].BURNSTAT = 1
            for ip in range(self.NPART):
                particle = self.particles[ip]
                X = self.grid.XCELL[IF, JF]
                Y = self.grid.YCELL[IF, JF]
                x, y = particle.x, particle.y
                if (abs(X - x) < self.grid.DX) and (abs(Y - y) < self.grid.DY):
                    particle.update(state=1, factor=LANGFACTOR)
        for k in IT:
            self.grid.CELLS[k, 98].BURNSTAT = 1
            for ip in range(self.NPART):
                particle = self.particles[ip]
                X = self.grid.XCELL[k, 98]
                Y = self.grid.YCELL[k, 98]
                x, y = particle.x, particle.y
                if (abs(X - x) < self.grid.DX) and (abs(Y - y) < self.grid.DY):
                    particle.update(state=1, factor=LANGFACTOR)

    def _move_particle(self, ip):
        """
            RANDOM WALK & DECAY STATUS TO ACCOUNT FOR MEMORY.
            IF PARTICLES LEAVE DOMAIN, IGNORE THEM.
            MOVE ONLY PARTICLES THAT ARE "ON".
            RANDOM WALK FOLLOWS CNF PAPER.

            NEW FEATURE: Multiply step by a FACTOR, which reflects the "age" of the
            particle relative to the clock of the fire in each cell. Hence, 0 causes
            no propagation, 1 full movement by the wind turbulence (note: in this
            implementation, FIRECONST = 1).

        """
        # RANDOM WALK
        # ADVANCE ONLY THE PARTICLES THAT ARE "ON" (i.e. ABOVE STTHR).
        #
        particle = self.particles[ip] # get particle
        props = ["state", "x", "y", "ux", "vy", "factor", "tmem"]
        state, x, y, ux, vy, factor, tmem = particle.get_from_keys(props)
        if state > STTHR:
            DU  = 0 -(ux - UXM)*2.0*TFREQ*self.DT + CLANG*self.SQRTDT*normal()
            DV  = 0 -(vy - VYM)*2.0*TFREQ*self.DT + CLANG*self.SQRTDT*normal()
            UXP = ux + DU
            VYP = vy + DV
            XP  = x + UXP*self.DT*factor
            YP  = y + VYP*self.DT*factor
            STP = state*np.exp(-self.DT/tmem)
            particle.update(ux=UXP, vy=VYP, x=XP, y=YP, state=STP)
        if x > self.grid.XMAX - self.grid.DX:
            particle.update(x=self.grid.XMAX - self.grid.DX, state=0.)
        elif x < self.grid.XMIN + self.grid.DX:
            particle.update(x=self.grid.XMIN + self.grid.DX, state=0.)
        if y > self.grid.YMAX - self.grid.DY:
            particle.update(y=self.grid.YMAX - self.grid.DY, state=0.)
        elif y < self.grid.YMIN + self.grid.DY:
            particle.update(y=self.grid.YMIN + self.grid.DY, state=0.)

    def _ignite_cells(self, istep, ip):
        """
            IGNITE CELLS VISITED BY BURNING PARTICLES, BUT ONLY IF VISITED
            FOR THE FIRST TIME.
            NOTE: WORKS ONLY FOR CONSTANT NUMBER OF PARTICLES (NPARTMAX), OR MUST
            RE-CALCULATE THE INDEX indp DIFFERENTLY.
            NOTE: STUPID WAY TO GET INDY & INDX, BUT USE OF "ROUND" DID NOT WORK. WORKS
            FOR STRUCTURED GRIDS ONLY.

            To see random walk alone, put state > STTHR*10 below
            to avoid doing any ignitions.

        """
        particle = self.particles[ip] # get particle
        state, x, y = particle.get_from_keys(["state", "x", "y"])
        if state > STTHR:
            for i in range(self.grid.NX-1):
                if abs(x - self.grid.XCELL[i, 0]) <= self.grid.DX/2:
                    INDX = i
            for j in range(self.grid.NY-1):
                if abs(y - self.grid.YCELL[0, j]) <= self.grid.DY/2:
                    INDY = j
            cell = self.grid.CELLS[INDX, INDY]
            if (cell.QMAXTR > 0) and cell.BURNSTAT == 0 and cell.BURNSTAT2 == 1:
                cell.BURNSTAT = 1
                cell.CLOCK    = self.TIME[istep]

    def _launch_particles(self, istep):
        """
            Launch new particles from the ignited cells. Randomise a bit the
            ignition delay time.
        """
        for i in range(self.grid.NX-1):
            for j in range(self.grid.NY-1):
                cell = self.grid.CELLS[i, j]
                TLOCAL = self.TIME[istep] - cell.CLOCK
                TCRIT  = cell.TIGNTR * (1 + RELT*normal())
                if cell.BURNSTAT == 1 and TLOCAL > TCRIT and cell.BURNSTAT2 == 1:
                    LOCALF = LANGFACTOR
                    indp = (i*(self.grid.NY - 1) + j)*Cell.NPARTMAX
                    cell.CLOCK = self.TIME[istep]
                    for k in range(cell.NPARTTR):
                        self.particles[k + indp].update(state=1.0, factor=LOCALF)
                    cell.BURNSTAT2 = 0

    def save_burning_state(self, istep):
        if self.ITYPESPARK == 1:
            for i in range(self.grid.NX-1):
                for j in range(self.grid.NY-1):
                    self.BURNX[istep, i] += self.grid.CELLS[i, j].BURNPROG
                self.BURNX[istep, i] /= self.grid.NY - 1
        self.BURNEVOL[istep, :, :] = self.grid.get_burn_state()[:, :]

    def propagate(self, istep):
        for ip in range(self.NPART):
            self._move_particle(ip)
            self._ignite_cells(istep, ip)
        self._launch_particles(istep)

    def extinguish(self, x, y):
        airplane = Airplane()
        Area, angle = airplane.calculate_params(x, y, USPEED)
        allx, ally, states = self.get_particles_props('x','y','state')
        allxrel = allx - x
        allyrel = ally - y
        A = np.array([[np.cos(angle), -np.sin(angle)],[np.sin(angle),np.cos(angle)]])

        for i in range(0,len(allx)):
            B = np.array([[allxrel[i]],[allyrel[i]]])
            rot_coord = np.linalg.solve(A,B)
            if states[i]>STTHR and (abs(rot_coord[0])<Area[1]*0.5) and (abs(rot_coord[1])<Area[0]*0.5):
                particle = self.particles[i] # get particle
                t = particle.get_from_keys(["tmem"])
                tm = t[0]*0.5
                self.particles[i].update(tmem=tm)
                
    def get_particles_props(self, *props, array=None):
        """Get all particles props given the props name.

        Args:
            array (np.ndarray or None, optional): the array of particles from
                                                  which we want to fetch the props.
                                                  Defaults to None.

        Returns:
            properties (np.ndarray): the properties of the particles in array
        """
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



    def scatter_plot_part(self, speed=0.0001):
        """Plot the current grid state."""
        FPARTX, FPARTY, STA = self.get_particles_props('x', 'y', 'state')
        plt.clf()
        b = []
        lg = []
        g = []
        grey = []
        r = []
        for i in range(self.grid.NX-1):
            for j in range(self.grid.NY-1):
                cell = self.grid.CELLS[i, j]
                if cell.QMAXTR == 0 and cell.BURNSTAT == 0:
                    b = b + [[self.grid.XCELL[i,j],self.grid.YCELL[i,j]]]                   
                elif cell.TIGNTR == 30 and cell.BURNSTAT == 0:
                    lg = lg + [[self.grid.XCELL[i,j],self.grid.YCELL[i,j]]] 
                elif cell.BURNSTAT == 0:
                    g = g + [[self.grid.XCELL[i,j],self.grid.YCELL[i,j]]]
                elif cell.BURNSTAT2 == 0: 
                    grey = grey + [[self.grid.XCELL[i,j],self.grid.YCELL[i,j]]]
                else:
                    r = r + [[self.grid.XCELL[i,j],self.grid.YCELL[i,j]]]
        b = np.array(b)
        lg = np.array(lg)
        g = np.array(g)
        grey = np.array(grey)
        r = np.array(r)
        plt.xlabel("Distance / m")
        plt.ylabel("Distance / m")
        plt.scatter(b[:,0],b[:,1],color = "tab:brown")
        plt.scatter(lg[:,0],lg[:,1],color = "lightgreen")
        plt.scatter(g[:,0],g[:,1],color = "green")
        if len(grey)>0:
            plt.scatter(grey[:,0],grey[:,1],color = "grey")
        plt.scatter(r[:,0],r[:,1],color = "orange", s = 10)
        plt.scatter(FPARTX[STA >= 0.2], FPARTY[STA >= 0.2], c = "darkred", s = 6)
        #plt.scatter(FPARTX, FPARTY, STA*5 + 0.01, c=STA, cmap="jet")
        #plt.clim(0, 1)
        #plt.colorbar()
        plt.pause(speed)

    def plot_igntcell(self, istep):
        Burning = np.zeros((self.grid.NX-1, self.grid.NY-1))
        for i in range(self.grid.NX-1):
            for j in range(self.grid.NY-1):
                cell = self.grid.CELLS[i, j]
                if cell.BURNSTAT2 == 0 and cell.BURNSTAT == 1:
                    Burning[j, i] = 100*(1 - (self.TIME[istep] - cell.CLOCK)/(cell.TBURNDUR/self.DT))
        FCELLX = self.grid.XCELL[:,0]
        FCELLY = self.grid.YCELL[0,:]
        PAX, PAY, PST = self.get_particles_props('x', 'y', 'state')
        plt.clf()
        plt.contourf(FCELLX, FCELLY, Burning, levels = 300, vmin=0, vmax=100, cmap="YlOrRd")
        plt.clim(0, 100)
        plt.colorbar(ticks=(0,25,50,75,100), extend="both", cmap="YlOrRd")
        for i in range(len(PAX)):
            if PST[i]>=0.2:
                plt.scatter(PAX[i],PAY[i], color="blue", s=2)
        plt.show()
        plt.pause(2)

    def launch(self):
        """Launch the simulation."""
        print(MESSAGE.format(grid=self.grid, nsteps=self.NSTEPS))

        for istep in range(self.NSTEPS):
            self.propagate(istep)
            if istep in(exting):
                Ml = MLmodel(self.grid,self.particles, UXM, VYM, self.grid.HIT)
                x,y = Ml.find_shortest()
                print(x,y)
                self.extinguish(x,y)
            #if istep == 0:
            self.scatter_plot_part()
        print("[END of simulation]")
        plt.show()

    def get_fire_line(self, epsilon=0.01, line_only=False):
        """Extract the fire line."""
        # Get active particles
        ACTIVE_PARTICLES = self.particles[self.particles > STTHR]
        # get maximum burning state
        MAX = max(self.particles, key=lambda p: p.state)
        Max = MAX.state

        # fetch  fire  line
        FIRE_LINE = ACTIVE_PARTICLES[Max - epsilon <= ACTIVE_PARTICLES]
        return FIRE_LINE if line_only else ACTIVE_PARTICLES, FIRE_LINE, MAX

    def interpolate_fire_line(self):
        """Interpolate fire line to get a fire particle state at each cell center."""
        fire_line = self.get_fire_line(line_only=True)
        LINEX, LINEY = self.get_particles_props('x', 'y', array=fire_line)
        coef = np.polyfit(LINEY, LINEX, deg=min(10, len(LINEX)))
        fit  = np.poly1d(coef)
        return fit(self.grid.YCELL[0])


    def plot_fire_line(self, epsilon=0.01, speed=10, stop=200):
        """ Plot the fire line at specific time :
            Launch a simulation and plot the fire line at each step.
        """
        stop = min(stop, self.NPART)
        # PLOTS setup
        _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 5))

        for istep in range(stop):
            self.propagate(istep)
            if istep % speed == 0:
                print(f"{istep = }")
                active, line, _ = self.get_fire_line(epsilon)
                XP, YP, ST = self.get_particles_props('x', 'y', 'state', array=active)
                LINEX, LINEY, LINEST = self.get_particles_props('x', 'y', 'state', array=line)

                coef = np.polyfit(LINEY, LINEX, deg=min(10, len(LINEX)))
                fit = np.poly1d(coef)
                data = fit(self.grid.YCELL[0])

                STinterp = np.interp(self.grid.YCELL[0], LINEY, LINEST)
                plt.cla()
                z1 = ax1.scatter(XP, YP, c=ST, cmap="jet")
                z2 = ax2.scatter(LINEX, LINEY, c=LINEST, cmap="jet")
                z3 = ax3.scatter(data, self.grid.YCELL[0], c=STinterp, cmap="jet")

                if istep == 0: plt.colorbar(z1)
                for ax, z in zip((ax1, ax2, ax3), (z1, z2, z3)):
                    ax.set_xlim(0, self.grid.XMAX)
                    z.set_clim(0, 1)
                plt.pause(1)
            if istep == stop - 1:
                print("[END of SIMULATION]")
                plt.pause(3)
                plt.close()
        plt.show()
