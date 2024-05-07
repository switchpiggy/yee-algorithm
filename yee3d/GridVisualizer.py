from yee3d.Grid import Grid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import imageio
import uuid

"""
GridVisualizer class: handles the visualization tasks. There are a wide range of ways to visualize a 3D FDTD simulation, such as 
using GIFs and 2D snapshots of the field.
"""

class GridVisualizer:
    def __init__(self, g: Grid = None, name:str = 'default') -> None:
        self.grid = g
        self.name = name
        self.savepaths = []

    def _return_field(self, field: str):
        if field == 'Ex':
            return self.grid.Ex
        elif field == 'Ey':
            return self.grid.Ey
        elif field == 'Ez':
            return self.grid.Ez
        elif field == 'Hx':
            return self.grid.Hx
        elif field == 'Hy':
            return self.grid.Hy
        elif field == 'Hz':
            return self.grid.Hz
        else:
            return None
        
    def cleanup(self):
        for path in self.savepaths:
            print(path)
            os.system(f'rm -rf {path}')

    def bindToGrid(self, g: Grid):
        self.grid = g

    def unbindGrid(self):
        self.grid = None

    def snapshot(self, dim: int = 0, slice_index: int = 0, field:str = 'Ex', path: str = '/', save: bool = False):
        matplotlib.use('agg')
        fig, ax = plt.subplots(1)
        F = self._return_field(field)

        if F == None:
            print("Field name invalid.")
            return
        
        if dim < 0 or dim > 2:
            print("Dimension invalid. Please enter an integer from 0 to 2. 0 => x, 1 => y, 2 => z")
            return
        
        if slice_index < 0 or slice_index >= self.grid.dims[dim] - 1:
            print("Slice index invalid.")
            return
        
        match dim:
            case 0:
                F = F[slice_index][:][:]
            case 1:
                F = F[:][slice_index][:]
            case 2:
                F = F[:][:][slice_index]
        
        im = ax.imshow(F, interpolation="nearest", cmap="plasma")
        cbar = ax.figure.colorbar(im, ax=ax)
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')

        fig.tight_layout()

        if save == True:
            fig.savefig(path)
        else:
            plt.show()
        
        plt.close()
        return
    
    def snapshotGIF(self, timesteps, period: int = 1, dim: int = 0, slice_index: int = 0, field: str = 'Ex'):
        matplotlib.use('agg')
        if self.grid == None:
            print("Visualizer not binded to grid.")
            return
        
        if timesteps > self.grid.maxTime:
            print("timesteps greater than allowed simulation time.")
            return 
        
        if dim < 0 or dim > 2:
            print("Invalid dimension, please choose integer 0 through 2.")
            return 
        
        if slice_index < 0 or slice_index >= self.grid.dims[dim] - 1:
            print("Invalid index for dimension, out of bounds.")
            return
        
        ID = uuid.uuid4()
        output_path = f'./grid-vis-output-im-{self.name}-{ID}/'

        if os.path.exists(output_path):
            os.system(f'rm -rf {output_path}')
        os.mkdir(output_path)
        self.savepaths.append(output_path)

        for i in range(timesteps):
            if i%period == 0:
                self.snapshot(dim=dim, slice_index=slice_index, field=field, path=output_path + f'{i//period}.png', save = True)
            self.grid.step()

        images = []
        for i in range(0, timesteps, period):
            # print(i)
            images.append(imageio.imread(output_path + f'{i//period}.png'))

        fullpath = output_path + '/result.gif'
        imageio.mimsave(fullpath, images)

        return fullpath



