import os
import joblib
import numpy as np
from particle import Particle


class MLmodel:
    def __init__(self, grid, particles, Ux, Uy, importance):
        self.grid       = grid
        self.particles  = particles
        self.area       = importance
        self.Wx         = Ux
        self.Wy         = Uy
        self.m          = self.Wy/self.Wx

        self.ppos       = []
        mindist         = np.sqrt(2*((self.grid.DX/2)**2))
        for i in range(len(self.particles)):
            p = self.particles[i]
            st, xp, yp = p.get_from_keys(["state", "x", "y"])
            if st > 0.2:
                c = yp - self.m*xp
                k = 0
                while k < len(self.area):
                    xa, ya = self.area[k,0], self.area[k,1]
                    k += 1
                    dist = (abs((-self.m*xa) + ya - c))/np.sqrt((self.m**2) + 1)
                    if dist <= mindist:
                        length = np.sqrt((xp-xa)**2 + (yp - ya)**2)
                        self.ppos = self.ppos + [[xp, yp, length]]
                        k += len(self.area)
        
        self.ppos = np.array(self.ppos)
        self.params = np.zeros((len(self.ppos),6))
        self.params[:,5] = np.sqrt((self.Wx)**2 + (self.Wy)**2)
        for j in range(len(self.ppos)):
            for i in range(self.grid.NX-1):
                if abs(self.ppos[j,0] - self.grid.XCELL[i, 0]) <= self.grid.DX/2:
                    INDX = i
            for k in range(self.grid.NY-1):
                if abs(self.ppos[j,1] - self.grid.YCELL[0, k]) <= self.grid.DY/2:
                    INDY = k

            cell = self.grid.CELLS[INDX, INDY]
            tau = cell.TIGNTR
            if tau > 30:
                self.params[j,3] = 1

            dx = np.sqrt(2*(10**2))
            it = int(self.ppos[j,2]/dx)
            x  = self.ppos[j,0]
            y  = self.ppos[j,1]
            change = 0
            l = 1
            Loop = True
            while Loop == True and l<it:
                x = x + dx/(np.sqrt(1 + (self.m)**2))
                y = y + dx/(self.m*np.sqrt(1 + (self.m)**2))
                for i in range(self.grid.NX-1):
                    if abs(x - self.grid.XCELL[i, 0]) <= self.grid.DX/2:
                        INDX = i
                for k in range(self.grid.NY-1):
                    if abs(y - self.grid.YCELL[0, k]) <= self.grid.DY/2:
                        INDY = k

                cell = self.grid.CELLS[INDX, INDY]
                tau2 = cell.TIGNTR
                if (tau2 != tau) and (change == 0):
                    L1 = l*dx
                    change = 1
                    self.params[j,0] = L1
                    tau = tau2
                    if tau > 30:
                        self.params[j,4] = 1
                elif (tau2 != tau) and (change == 1):
                    L2 = l*dx
                    change = 2
                    self.params[j,1] = L2
                    L3 = self.ppos[j,2] - L2 - L1
                    self.params[j,2] = L3
                    Loop = False
                l += 1
            if change == 0:
                self.params[j,0] = self.ppos[j,2]
            elif change == 1:
                self.params[j,1] = self.ppos[j,2] - L1


    def find_shortest(self):

        RFRModel = joblib.load("./RFRmodel.joblib")
        scale = np.array([[1401,1439.8,1193.8,0.577,0.723,49.3],[1587.92,1446.76,1382.59,0.4939,0.7661,30.02]])
        paramsc = []

        for i in range(len(self.params)):
            paramsc.append((self.params[i]-scale[0])/scale[1])

        paramsc = np.array(paramsc)
        t_predic = RFRModel.predict(paramsc)
        min = 100000000000

        for i in range(len(t_predic)):
            if t_predic[i] < min:
                min = t_predic[i]
                ind = i      
        x, y = self.ppos[ind,0], self.ppos[ind,1]

        return x, y 