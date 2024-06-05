from grid import Grid
from lagrange import Lagrange


# INITIALIZATION OF THE GRID
grid = Grid(XMAX=1000, YMAX=1000)

# SET NON FLAMMABLE CELLS
# using ratio = 0.5 (ratio is the fraction of non flammable cells);
# this also put few roads.
grid.set_layout()

# INITIALIZE THE LAGRAGIAN MODEL ON THE GRID
# we can specify parameter NSTEPS if we want but by defaut it
# is set to 400
model = Lagrange(grid, NSTEPS=700)

# START FIRE
# ITYPESPARK = 0 (point)  1 (line)
model.start_fire(ITYPESPARK=1)

# LAUNCH SIMULATION
model.launch()

#model.plot_fire_line()
