import os
import numpy as np
import iris
import iris.analysis.calculus as calculus
from mogreps import data_paths as mogreps_data_paths
def derivative(cube, axisname):
    R = 6378388.  # Radius of the earth
    deg2rad = 0.0174533  # pi/180.
    dcube = cube.copy()

    coord_names = np.array([c.var_name for c in cube.coords()])
    # print(coord_names)

    if axisname == 'latitude':
        lats = cube.coord('latitude').points
        axis_index = np.where(coord_names == 'latitude')[0][0]
        dlat = np.diff(lats) * deg2rad  # convert to radians
        dy = R * np.sin(dlat)  # constant at this latitude
        dcube = calculus.differentiate(cube, 'latitude')
        dcube /= iris.util.broadcast_to_shape(dy, dcube.shape, (axis_index,))

    if axisname == 'longitude':
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        axis_index = np.where(coord_names == 'latitude')[0][0]
        dlon = (lons[1] - lons[0]) * deg2rad  # convert to radians
        dx = np.array([R * np.cos(deg2rad * lat) * dlon for lat in lats])
        dcube = calculus.differentiate(cube, 'longitude')
        dcube /= iris.util.broadcast_to_shape(dx, dcube.shape, (axis_index,))

    return dcube


def compute_vort_div(forecast_date, mem=0):
    fc_dt = forecast_date
    print(fc_dt)
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    # hour is forecast hour
    digit3_mem = str('%03d' % mem)

    wave_names = np.array(['Kelvin', 'WMRG', 'R1', 'R2'])

    for wname in wave_names:
        print(wname)
        u_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
        'u_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                        str_month, str_day, str_hr,
                                        digit3_mem))

        v_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
        'v_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                        str_month, str_day, str_hr,
                                        digit3_mem))


        if os.path.exists(u_wave_data_file):
            print(u_wave_data_file)
            u_wave = iris.load_cube(u_wave_data_file)
        if os.path.exists(v_wave_data_file):
            v_wave = iris.load_cube(v_wave_data_file)

        # compute vorticity and divergence

        # Divergence
        div_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                        'div_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                        str_month, str_day, str_hr,
                                                                        digit3_mem))

        div = derivative(u_wave, 'longitude').regrid(u_wave, iris.analysis.Linear())
        div += derivative(v_wave, 'latitude').regrid(u_wave, iris.analysis.Linear())
        #print(div)
        # Write the file out
        iris.save(div, div_wave_data_file)
        print('Written %s' %div_wave_data_file)

        # Vorticity
        vort_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                          'vort_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                            str_month, str_day, str_hr,
                                                                            digit3_mem))

        vort = derivative(v_wave, 'longitude').regrid(u_wave, iris.analysis.Linear())
        vort -= derivative(u_wave, 'latitude').regrid(u_wave, iris.analysis.Linear())

        # Write the file out
        iris.save(vort, vort_wave_data_file)
        print('Written %s' % vort_wave_data_file)
