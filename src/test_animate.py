import numpy as np
import matplotlib.pyplot as plt
import iris
import cartopy.crs as ccrs
from matplotlib import animation
import sys

# animation function
def animate(i):
    z = cubes[i].data
    print(i)
    cont = plt.contourf(lons, lats, z, levels=z_levels, cmap='RdBu', extend='both')
    plt.title(cubes[i].coord('time'))
    return cont



file= '/project/MJO_GCSS/Eqwaves_monitoring/processed_wave_data/z_wave_Kelvin_20191105_12Z_000.nc'
cubes = iris.load_cube(file)[:,0]
ntime = len(cubes.coord('time').points)
print(cubes)

ntime, nlon, nlat = cubes.shape
lons = cubes.coord('longitude').points
lats = cubes.coord('latitude').points
lons,lats = np.meshgrid(lons,lats)
print(lons.shape, cubes.data.shape)

proj = ccrs.PlateCarree(central_longitude=-180.0)

fig = plt.figure(figsize=(12, 5))
plt.subplot(211, projection=proj)
z_levels = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
plt.gca().coastlines()
#ax = plt.axes(xlim=(0, 360), ylim=(-24, 24))

anim = animation.FuncAnimation(fig, animate,
                           frames=ntime, interval=50)


plt.show()
#cube_iter = cubes.slices(('longitude', 'latitude'))
#ani = animate(cube_iter, qplt.contourf(cube_iter, vmin=-20), interval=50)
#plt.show()
