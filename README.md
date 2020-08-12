# Scorpio_AMR_Loader
the external Scorpio frontend of YT

Usage
```python
import sys
sys.path.append("/Path/of/the/package")
from data_io import ScorpioLoader
variables = [1,1,1,1,1,1,1,1]  # variables in the simulation ['den','momx','momy','momz','bx','by','bz','ene']
dimensions = [100,100,150]     # initial resolution along x,y,z direction    
periodicity = [False]*3        # whether the boundary condition is periodic
flnm = './g0000.h5'            # filename of the simulation
loader = ScorpioLoader(flnm,variables, dimensions,periodicity)
ds = loader.load()             # this will give you a YT dataset, for more tutorial, check https://yt-project.org/
```
