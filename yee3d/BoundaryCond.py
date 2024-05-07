from yee3d.Grid import Grid

class BoundaryCond:
    def __init__(self, g: Grid = None, name: str = None) -> None:
        self.grid = g
        self.name = name

    def update(self):
        return
    
    def bindToGrid(self, g: Grid):
        self.grid = g
        self.grid.boundaryConditions.append(self)

    def unbindGrid(self):
        self.grid = None

class FirstOrderABC(BoundaryCond):
    def __init__(self, g: Grid = None, name: str = 'first_order_ABC') -> None:
        super().__init__(name)
        if self.grid != None:
            self.bindToGrid(g)
    
    def bindToGrid(self, g: Grid):
        super().bindToGrid(g)
        self.EzXM = self.EzX0 = self.EyX0 = self.EyXM = [[0] * self.grid.dims[2] for _ in range(self.grid.dims[1])]
        self.ExY0 = self.EzY0 = self.ExYN = self.EzYN = [[0] * self.grid.dims[2] for _ in range(self.grid.dims[0])]
        self.EyZ0 = self.ExZ0 = self.EyZP = self.ExZP = [[0] * self.grid.dims[1] for _ in range(self.grid.dims[0])]
        self.coef = (self.grid.SC - 1)/(self.grid.SC + 1)

    def unbindGrid(self):
        return super().unbindGrid()
    
    def update(self):
        #update Ez, Ey for X0, XM
        for y in range(self.grid.dims[1]):
            for z in range(self.grid.dims[2]):
                if y < self.grid.dims[1] - 1:
                    #update Ey
                    self.grid.Ey[0][y][z] = self.EyX0[y][z] + self.coef * (self.grid.Ey[1][y][z] - self.grid.Ey[0][y][z])
                    self.EyX0[y][z] = self.grid.Ey[0][y][z]
                    self.grid.Ey[self.grid.dims[0] - 1][y][z] = self.EyXM[y][z] + self.coef * \
                        (self.grid.Ey[self.grid.dims[0] - 2][y][z] - self.grid.Ey[self.grid.dims[0] - 1][y][z])
                    self.EyXM[y][z] = self.grid.Ey[self.grid.dims[0] - 1][y][z]
                if z < self.grid.dims[2] - 1:
                    #update Ez
                    self.grid.Ez[0][y][z] = self.EzX0[y][z] + self.coef * (self.grid.Ez[1][y][z] - self.grid.Ez[0][y][z])
                    self.EzX0[y][z] = self.grid.Ez[0][y][z]
                    self.grid.Ez[self.grid.dims[0] - 1][y][z] = self.EzXM[y][z] + self.coef * \
                        (self.grid.Ez[self.grid.dims[0] - 2][y][z] - self.grid.Ez[self.grid.dims[0] - 1][y][z])
                    self.EzXM[y][z] = self.grid.Ez[self.grid.dims[0] - 1][y][z]

        #update Ex, Ez for Y0, YM
        for x in range(self.grid.dims[0]):
            for z in range(self.grid.dims[2]):
                if x < self.grid.dims[0] - 1:
                    #update Ey
                    self.grid.Ex[x][0][z] = self.ExY0[x][z] + self.coef * (self.grid.Ex[x][1][z] - self.grid.Ex[x][0][z])
                    self.ExY0[x][z] = self.grid.Ex[x][0][z]
                    self.grid.Ex[x][self.grid.dims[1] - 1][z] = self.ExYN[x][z] + self.coef * (self.grid.Ex[x][self.grid.dims[1] - 2][z] - self.grid.Ex[x][self.grid.dims[1] - 1][z])
                    self.ExYN[y][z] = self.grid.Ex[x][self.grid.dims[1] - 1][z]
                if z < self.grid.dims[2] - 1:
                    #update Ez
                    self.grid.Ez[x][0][z] = self.EzY0[x][z] + self.coef * (self.grid.Ez[x][1][z] - self.grid.Ez[x][0][z])
                    self.EzY0[y][z] = self.grid.Ez[x][0][z]
                    self.grid.Ez[x][self.grid.dims[1] - 1][z] = self.EzYN[x][z] + self.coef * (self.grid.Ez[x][self.grid.dims[1] - 2][z] - self.grid.Ez[x][self.grid.dims[1] - 1][z])
                    self.EzYN[y][z] = self.grid.Ez[x][self.grid.dims[1] - 1][z]
                    
        #update Ey, Ex for Z0, Z
        for x in range(self.grid.dims[0]):
            for y in range(self.grid.dims[1]):
                if x < self.grid.dims[0] - 1:
                    #update Ey
                    self.grid.Ex[x][y][0] = self.ExY0[x][y] + self.coef * (self.grid.Ex[x][y][1] - self.grid.Ex[x][y][0])
                    self.ExZ0[x][y] = self.grid.Ex[x][y][0]
                    self.grid.Ex[x][y][self.grid.dims[2] - 1] = self.ExZP[x][y] + self.coef * (self.grid.Ex[x][y][self.grid.dims[2] - 2] - self.grid.Ex[x][y][self.grid.dims[2] - 1])
                    self.ExZP[y][z] = self.grid.Ex[x][y][self.grid.dims[2] - 1]
                if y < self.grid.dims[1] - 1:
                    #update Ez
                    self.grid.Ey[x][y][0] = self.EyZ0[x][y] + self.coef * (self.grid.Ey[x][y][1] - self.grid.Ey[x][y][0])
                    self.EyZ0[x][y] = self.grid.Ey[x][y][0]
                    self.grid.Ey[x][y][self.grid.dims[2] - 1] = self.EyZP[x][y] + self.coef * (self.grid.Ey[x][y][self.grid.dims[2] - 2] - self.grid.Ey[x][y][self.grid.dims[2] - 1])
                    self.EyZP[y][z] = self.grid.Ey[x][y][self.grid.dims[2] - 1]

class TFSF(BoundaryCond):
    def __init__(self, g: Grid = None, name: str = None) -> None:
        super().__init__(g, name)