import numpy as np

class Airplane:

    def __init__(self, speed=66, DoorArea=3.45, volume=4.55, flowrate=2.1, waterx=0, watery=0):
        self.speed     = speed
        self.DoorArea  = DoorArea
        self.volume    = volume
        self.flowrate  = flowrate
        self.waterx    = waterx
        self.watery    = watery
        self.occupied  = 0
        self.U1        = self.flowrate/self.DoorArea
        self.reload    = self.volume/self.flowrate
        self.ext       = 0.5
        
        self.f1        = 2     # legendre paper
        self.f2        = 27    # 27 for gravitiy, 58 for pressurised jet
        self.rhoa      = 1.225
        self.rhol      = 1000

    def calculate_params(self, target_x, target_y, windspeed):

        flight_angle = np.tan((self.watery - target_y)/(self.waterx - target_x))
        Urel = self.speed - windspeed*np.cos(flight_angle)
        q = (self.rhol*(self.U1**2))/(self.rhoa*Urel)
        Lambda = 5*(self.DoorArea**0.5)*self.f2*(q**0.2)
        droplength = 5*(self.speed*self.reload + self.f1*Lambda)
       
        Area = np.array([Lambda, droplength])
        return Area, flight_angle