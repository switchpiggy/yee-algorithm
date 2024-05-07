from typing import List
import numpy as np
import os
import matplotlib.pyplot as plt

"""
Auxiliary function for calculating the dipole source
"""

def dipoleSourceCalc(time, location) -> float:
    SC = 1/np.sqrt(3)
    # print(np.pi * ((SC * time - location)/1000 - 1.0))
    arg = (np.pi * ((SC * time - location)/0.1 - 1.0)) ** 2
    return (1 - 2 * arg) * np.exp(-arg)

"""
Auxiliary function for calculating sinusoidal source
"""

def harmonicSourceCalc(time, location, k: float = 1, A: float = 1) -> float:
    SC = 1/np.sqrt(3)
    return A * np.sin(2.0 * np.pi / k * (SC * time - location));


"""
Grid class: most of the heavy lifting is done here. The Ampere's Law and Faraday's Law update equations are handled by the Grid class
    as well as putting together the boundary conditions, objects, and sources.
"""

class Grid:
    def _init_3d_array(self, val, sizeX, sizeY, sizeZ):
        return [[[val] * sizeZ for _ in range(sizeY)] for _ in range(sizeX)]
        
    #initialize homogeneous grid of size (sx, sy, sx) 
    def __init__(self, sx: int, sy: int, sz: int, maxTime: int, permittivity: float = 1, permeability: float = 1, boundaryConditions: List = []):
        self.dims = (sx, sy, sz)
        self.maxTime = maxTime
        self.time = 0
        self.permittivity = permittivity #relative permittivity of free space
        self.permeability = permeability #relative permeability of free space
        self.imp0 = 377.0 * np.sqrt(self.permeability/self.permittivity) #impedance of free space, equal to sqrt(mu/epsilon)
        self.SC = 1/np.sqrt(3)

        #initialize E and H fields
        self.Ex = self._init_3d_array(0, sx - 1, sy, sz)
        self.Ey = self._init_3d_array(0, sx, sy - 1, sz)
        self.Ez = self._init_3d_array(0, sx, sy, sz - 1)

        self.Hx = self._init_3d_array(0, sx, sy - 1, sz - 1)
        self.Hy = self._init_3d_array(0, sx - 1, sy, sz - 1)
        self.Hz = self._init_3d_array(0, sx - 1, sy - 1, sz)

        #initialize update coefficients for homo grid
        self.ExE = self._init_3d_array(1, sx - 1, sy, sz)
        self.EyE = self._init_3d_array(1, sx, sy - 1, sz)
        self.EzE = self._init_3d_array(1, sx, sy, sz - 1)

        self.HxH = self._init_3d_array(1, sx, sy - 1, sz - 1)
        self.HyH = self._init_3d_array(1, sx - 1, sy, sz - 1)
        self.HzH = self._init_3d_array(1, sx - 1, sy - 1, sz)

        self.ExH = self._init_3d_array(self.SC * self.imp0, sx - 1, sy, sz)
        self.EyH = self._init_3d_array(self.SC * self.imp0, sx, sy - 1, sz)
        self.EzH = self._init_3d_array(self.SC * self.imp0, sx, sy, sz - 1)

        self.HxE = self._init_3d_array(self.SC / self.imp0, sx, sy - 1, sz - 1)
        self.HyE = self._init_3d_array(self.SC / self.imp0, sx - 1, sy, sz - 1)
        self.HzE = self._init_3d_array(self.SC / self.imp0, sx - 1, sy - 1, sz)

        #boundary conditions (none yet)
        self.boundaryConditions = boundaryConditions

    def getTime(self):
        return self.time
    
    def getMaxTime(self):
        return self.maxTime
    
    def shape(self):
        return self.dims

    def updateE(self):
        #take curl of H and update E
        for i in range(0, self.dims[0] - 1):
            for j in range(1, self.dims[1] - 1):
                for k in range(1, self.dims[2] - 1):
                    self.Ex[i][j][k] = self.ExE[i][j][k] * self.Ex[i][j][k] + self.ExH[i][j][k] * \
                    ((self.Hz[i][j][k] - self.Hz[i][j - 1][k]) - (self.Hy[i][j][k] - self.Hy[i][j][k - 1]))
        
        for i in range(1, self.dims[0] - 1):
            for j in range(0, self.dims[1] - 1):
                for k in range(1, self.dims[2] - 1):
                    self.Ey[i][j][k] = self.EyE[i][j][k] * self.Ey[i][j][k] - self.EyH[i][j][k] * ((self.Hz[i][j][k] - self.Hz[i - 1][j][k]) - (self.Hx[i][j][k] - self.Hx[i][j][k - 1]))
        
        for i in range(1, self.dims[0] - 1):
            for j in range(1, self.dims[1] - 1):
                for k in range(0, self.dims[2] - 1):
                    self.Ex[i][j][k] = self.EzE[i][j][k] * self.Ez[i][j][k] + self.EzH[i][j][k] * ((self.Hy[i][j][k] - self.Hy[i - 1][j][k]) - (self.Hx[i][j][k] - self.Hx[i][j - 1][k]))
        #apply boundary conditions
    
    def updateH(self):
        for i in range(0, self.dims[0]):
            for j in range(0, self.dims[1] - 1):
                for k in range(0, self.dims[2] - 1):
                    # print(i)
                    self.Hx[i][j][k] = self.HxH[i][j][k] * self.Hx[i][j][k] + self.HxE[i][j][k] * ((self.Ez[i][j + 1][k] - self.Ez[i][j][k]) - (self.Ey[i][j][k + 1] - self.Ey[i][j][k]))
        
        for i in range(0, self.dims[0] - 1):
            for j in range(0, self.dims[1]):
                for k in range(0, self.dims[2] - 1):
                    self.Hy[i][j][k] = self.HyH[i][j][k] * self.Hy[i][j][k] - self.HyE[i][j][k] * ((self.Ez[i + 1][j][k] - self.Ez[i][j][k]) - (self.Ex[i][j][k + 1] - self.Ex[i][j][k]))
        
        for i in range(0, self.dims[0] - 1):
            for j in range(0, self.dims[1] - 1):
                for k in range(0, self.dims[2]):
                    self.Hz[i][j][k] = self.HzH[i][j][k] * self.Hz[i][j][k] + self.HzE[i][j][k] * ((self.Ey[i + 1][j][k] - self.Ey[i][j][k]) - (self.Ex[i][j + 1][k] - self.Ex[i][j][k]))

    def applyBoundaryConditions(self):
        for bc in self.boundaryConditions:
            bc.update()
    
    def step(self):
        self.time += 1
        self.updateH()
        self.updateE()
        self.Ex[(self.dims[0] - 1)//2][(self.dims[1] - 1)//2][(self.dims[2] - 1)//2] += dipoleSourceCalc(self.time, 0)
        self.applyBoundaryConditions()

    def reset(self):
        self.__init__(self.dims[0], self.dims[1], self.dims[2], self.maxTime, self.permittivity, self.permeability, self.boundaryConditions)

    def run(self, endTime: int = None):
        if endTime == None:
            endTime = self.maxTime

        if self.time >= self.maxTime:
            print("Simulation already over. Please reset in order to restart the grid.")
            return
        
        if endTime > self.maxTime:
            print("End time is greater than the maximum allowed time. Please choose a different end time.")
            return
        
        while self.time < endTime:
            self.step()