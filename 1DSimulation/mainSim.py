from gridSim import Grid
from lagrangeSim import Lagrange
from ArraySim import Parameters
import xlsxwriter
import numpy as np

# CREATE ARRAY OF PARAMETERS
X = Parameters()
Xf = np.array([X[0:3000]]).reshape(3000,6)
for i in range(3000,len(X)):
    if X[i,4] != 2:
        Xf = np.vstack([Xf,X[i]])

#Length = 11212

# CREATE THE EXCEL WORKBOOK WITH THE REQUIRED COLUMNS
workbook = xlsxwriter.Workbook('10251-11212.xlsx')
worksheet = workbook.add_worksheet()
worksheet.write('A1', 'L1')
worksheet.write('B1', 'L2')
worksheet.write('C1', 'L3')
worksheet.write('D1', 'T1')
worksheet.write('E1', 'T2')
worksheet.write('F1', 'Ux')                  
worksheet.write('G1', 'time')

for i in range(10250,len(Xf)):
    # INITIALIZATION OF THE GRID
    grid = Grid(l1=Xf[i,0], l2=Xf[i,1], l3=Xf[i,2], T1=Xf[i,3], T2=Xf[i,4])

    # INITIALIZE THE LAGRAGIAN MODEL ON THE GRID
    model = Lagrange(grid, U=Xf[i,5])

    # START FIRE
    model.start_fire()

    # LAUNCH SIMULATION AND GET TIME
    time = model.launch()

    # WRITE PARAMETERS AND TIME INTO WORKBOOK
    worksheet.write(i+1,0, Xf[i,0])
    worksheet.write(i+1,1, Xf[i,1])
    worksheet.write(i+1,2, Xf[i,2])
    worksheet.write(i+1,3, Xf[i,3])
    worksheet.write(i+1,4, Xf[i,4])
    worksheet.write(i+1,5, Xf[i,5])
    worksheet.write(i+1,6, time)
    print(i)
    
workbook.close()
    

