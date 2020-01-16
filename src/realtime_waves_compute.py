import os
import iris
import datetime
import numpy as np
from mogreps import data_paths as mogreps_data_paths


def read_data(filenames, latmax=None):
    ### read in u,v,z data from file could be replaced by other routine
    if latmax == None:
        latmax = 24

    u_file, v_file, z_file = filenames

    # Reading the data using Iris
    u = iris.load_cube(u_file, 'x_wind').intersection(latitude=(-latmax, latmax))
    v = iris.load_cube(v_file, 'y_wind').intersection(latitude=(-latmax, latmax))
    z = iris.load_cube(z_file, 'geopotential_height').intersection(latitude=(-latmax, latmax))

    # Changing the coordinate names to time and pressure
    #u.coord('t').long_name = 'time'
    #u.coord('Pressure').long_name = 'pressure'

    # copy the coords over to v and z
    #v.coords = z.coords = u.coords

    return u, v, z


def uz_to_qr(u, z, g_on_c):
    # transform u,z to q, r using q=z*(g/c) + u; r=z*(g/c) - u
    q = z * g_on_c + u
    r = z * g_on_c - u
    return q, r


def filt_project(qf, rf, vf, lats, y0, waves, pmin, pmax, kmin, kmax, c_on_g):
    '''
    Author: Steve J. Woolnough (November 2019)
    :param qf:
    :param rf:
    :param vf:
    :param lats:
    :param y0:
    :param waves:
    :param pmin:
    :param pmax:
    :param kmin:
    :param kmax:
    :param c_on_g:
    :return:
    '''

    # find size of arrays
    nf = qf.shape[0]
    nz = qf.shape[1]
    nlats = lats.size
    nk = qf.shape[3]

    # Find frequencies and wavenumbers corresponding to pmin,pmax and kmin,kmax in coeff matrices
    f = np.fft.fftfreq(nf, 0.25)
    fmin = np.where((f >= 1. / pmax))[0][0]
    fmax = (np.where((f > 1. / pmin))[0][0]) - 1
    f1p = fmin
    f2p = fmax + 1
    f1n = nf - fmax
    f2n = nf - fmin + 1
    k1p = kmin
    k2p = kmax + 1
    k1n = nk - kmax
    k2n = nk - kmin + 1

    ### note that need to adjust for the fact that array referencing doesn't include last point!!!!!!

    # Define the parobolic cylinder functions
    spi2 = np.sqrt(2 * np.pi)
    dsq = np.array([spi2, spi2, 2 * spi2, 6 * spi2])  # Normalization for the 1st 4 paroblic CF
    d = np.zeros([dsq.size, nlats])
    y = lats[:] / y0
    ysq = y ** 2
    d[0, :] = np.exp(-ysq / 4.0)
    d[1, :] = y * d[0, :]
    d[2, :] = (ysq - 1.0) * d[0, :]
    d[3, :] = y * (ysq - 3.0) * d[0, :]
    dlat = np.abs(lats[0] - lats[1])

    qf_Kel = np.zeros([nf, nz, nk], dtype=complex)
    qf_mode = np.zeros([dsq.size, nf, nz, nk], dtype=complex)
    rf_mode = np.zeros([dsq.size, nf, nz, nk], dtype=complex)
    vf_mode = np.zeros([dsq.size, nf, nz, nk], dtype=complex)

    # reorder the spectral coefficents to make the latitudes the last dimension
    qf = np.transpose(qf, [0, 1, 3, 2])
    rf = np.transpose(rf, [0, 1, 3, 2])
    vf = np.transpose(vf, [0, 1, 3, 2])

    for m in np.arange(dsq.size):
        if m == 0:  # For eastward moving Kelvin Waves
            qf_Kel[f1n:f2n, :, k1p:k2p] = np.sum(qf[f1n:f2n, :, k1p:k2p, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
            qf_Kel[f1p:f2p, :, k1n:k2n] = np.sum(qf[f1p:f2p, :, k1n:k2n, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
        # For westward moving waves
        qf_mode[m, f1n:f2n, :, k1n:k2n] = np.sum(qf[f1n:f2n, :, k1n:k2n, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
        qf_mode[m, f1p:f2p, :, k1p:k2p] = np.sum(qf[f1p:f2p, :, k1p:k2p, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
        rf_mode[m, f1n:f2n, :, k1n:k2n] = np.sum(rf[f1n:f2n, :, k1n:k2n, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
        rf_mode[m, f1p:f2p, :, k1p:k2p] = np.sum(rf[f1p:f2p, :, k1p:k2p, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
        vf_mode[m, f1n:f2n, :, k1n:k2n] = np.sum(vf[f1n:f2n, :, k1n:k2n, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)
        vf_mode[m, f1p:f2p, :, k1p:k2p] = np.sum(vf[f1p:f2p, :, k1p:k2p, :] \
                                                 * np.squeeze(d[m, :]) * dlat, axis=-1) / (dsq[m] * y0)

    uf_wave = np.zeros([waves.size, nf, nz, nlats, nk], dtype=complex)
    zf_wave = np.zeros([waves.size, nf, nz, nlats, nk], dtype=complex)
    vf_wave = np.zeros([waves.size, nf, nz, nlats, nk], dtype=complex)

    for w in np.arange(waves.size):
        if waves[w] == 'Kelvin':
            for j in np.arange(nlats):
                uf_wave[w, :, :, j, :] = \
                    0.5 * qf_Kel[:, :, :] * d[0, j]
                zf_wave[w, :, :, j, :] = \
                    0.5 * qf_Kel[:, :, :] * d[0, j] * c_on_g

        if waves[w] == 'WMRG':
            for j in np.arange(nlats):
                uf_wave[w, :, :, j, :] = \
                    0.5 * qf_mode[1, :, :, :] * d[1, j]
                zf_wave[w, :, :, j, :] = \
                    0.5 * qf_mode[1, :, :, :] * d[1, j] * c_on_g
                vf_wave[w, :, :, j, :] = \
                    vf_mode[0, :, :, :] * d[0, j]

        if waves[w] == 'R1':
            for j in np.arange(nlats):
                uf_wave[w, :, :, j, :] = \
                    0.5 * (qf_mode[2, :, :, :] * d[2, j] - rf_mode[0, :, :, :] * d[0, j])
                zf_wave[w, :, :, j, :] = \
                    0.5 * (qf_mode[2, :, :, :] * d[2, j] + rf_mode[0, :, :, :] * d[0, j]) * c_on_g
                vf_wave[w, :, :, j, :] = \
                    vf_mode[1, :, :, :] * d[1, j]

        if waves[w] == 'R2':
            for j in np.arange(nlats):
                uf_wave[w, :, :, j, :] = \
                    0.5 * (qf_mode[3, :, :, :] * d[3, j] - rf_mode[1, :, :, :] * d[1, j])
                zf_wave[w, :, :, j, :] = \
                    0.5 * (qf_mode[3, :, :, :] * d[3, j] + rf_mode[1, :, :, :] * d[1, j]) * c_on_g
                vf_wave[w, :, :, j, :] = \
                    vf_mode[2, :, :, :] * d[2, j]

    return uf_wave, zf_wave, vf_wave



def makes_5d_cube(data, wave_names, time_coord, pressure_coord,
                  lat_coord, lon_coord):
    # ===========================================================================
    # # Make a 3D cube of Latitude, wavenumber & frequency dimensions
    # ===========================================================================
    var_cube = iris.cube.Cube(data)
    var_cube.rename('wave_anomalies')
    wave_coord = iris.coords.DimCoord(range(len(wave_names)), long_name='wave_name')
    wave_coord.guess_bounds()
    # ['Kelvin', 'WMRG', 'R1', 'R2']
    var_cube.attributes = {'Kelvin': 0, 'WMRG': 1, 'R1': 2, 'R2': 3}
    wave_coord.attributes = {'Kelvin': 0, 'WMRG': 1, 'R1': 2, 'R2': 3}
    var_cube.add_dim_coord(wave_coord, 0)
    var_cube.add_dim_coord(time_coord, 1)
    var_cube.add_dim_coord(pressure_coord, 2)
    var_cube.add_dim_coord(lat_coord, 3)
    var_cube.add_dim_coord(lon_coord, 4)
    return var_cube


def write_data(wave_cube, forecast_date, var_name = 'u_wave', mem=0):
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    # hour is forecast hour
    digit3_mem = str('%03d' % mem)
    digit2_mem = str('%02d' % mem)

    # Realistically you will probably only want to write out say (T-4:T+7) so that
    # you can plot an animation of the last few days and the forecast
    # total of 45 time points
    write_out_times = 45

    wave_names = wave_cube.coord('wave_name').attributes
    # Writing files out for each wave separately
    for wname, index in wave_names.items():
        print(wname, index)
        var_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              '%s_%s_%s%s%s_%sZ_%s.nc' % (var_name, wname, str_year,
                                                          str_month, str_day, str_hr,
                                                          digit3_mem))
        iris.save(wave_cube[index, -write_out_times:], var_file_out,
                  netcdf_format="NETCDF3_CLASSIC")
        print(var_file_out)

def compute_wave_data(forecast_date,
                      ntimes_total=360,
                      ntimes_analysis=332,
                      ntimes_forecast=28,
                      mem=0):
    # hour is forecast hour
    digit3_mem = str('%03d' % mem)
    digit2_mem = str('%02d' % mem)

    # ---------------------------------------------------------
    # start of main routine

    # define some physical parameters
    # ----------------------------------
    g = 9.8
    beta = 2.3e-11
    radea = 6.371e6
    spd = 86400.
    ww = 2 * np.pi / spd

    # Define some parameters spefic to the methodology
    latmax = 24.    # +/- latitude range over which to process data.
    kmin = 2        # minimum zonal wavenumber
    kmax = 40       # maximum zonal wavenumber
    pmin = 2.0      # minimum period (in days)
    pmax = 30.0     # maximum period (in days)
    y0 = 6.0        # meridional trapping scale (degrees)
    wave_names = np.array(['Kelvin', 'WMRG', 'R1', 'R2'])  # List of wave types to output

    y0real = 2 * np.pi * radea * y0 / 360.0  # convert trapping scale to metres
    ce = 2 * y0real ** 2 * beta
    g_on_c = g / ce
    c_on_g = ce / g
    # ----------------------------------

    # read in 90 days of u,v,z data at 6 hourly time resolution and equatorial grid point
    # Methodology test and applied with 1 degree resolution data, but not dependent on it
    # For Met Office operational forecasts the data would be 83 days of analysis and 7 days of forecast

    # read data
    # Now read the forecasts
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    u_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'uwnd_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))
    v_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'vwnd_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))
    z_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'zhgt_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))

    u, v, z = read_data([u_file, v_file, z_file], latmax=24)
    lons = u.coord('longitude')
    lats = u.coord('latitude')
    press = u.coord('pressure')
    time_coord = u.coord('time')
    times = [time_coord.units.num2date(tp).strftime("%Y%m%d%h") for tp in time_coord.points]

    # convert u,z to q,r
    q, r = uz_to_qr(u.data, z.data, g_on_c)

    # Fourier transform in time and longitude
    qf = np.fft.fft2(q.data, axes=(0, -1))
    rf = np.fft.fft2(r.data, axes=(0, -1))
    vf = np.fft.fft2(v.data, axes=(0, -1))

    # Project onto individual wave modes
    uf_wave, zf_wave, vf_wave = filt_project(qf, rf, vf, lats.points, y0, wave_names, pmin, pmax, kmin, kmax, c_on_g)

    # Inverse Fourier transform in time and longitude
    u_wave = np.real(np.fft.ifft2(uf_wave, axes=(1, -1)))
    z_wave = np.real(np.fft.ifft2(zf_wave, axes=(1, -1)))
    v_wave = np.real(np.fft.ifft2(vf_wave, axes=(1, -1)))

    print(u_wave[:, 0, 0, 0, 0])
    print(z_wave[:, 0, 0, 0, 0])
    print(v_wave[:, 0, 0, 0, 0])
    print(u_wave.shape)

    # Make iris cubes before writing the data
    u_wave_cube = makes_5d_cube(u_wave, wave_names, time_coord, press, lats, lons)
    v_wave_cube = makes_5d_cube(v_wave, wave_names, time_coord, press, lats, lons)
    z_wave_cube = makes_5d_cube(z_wave, wave_names, time_coord, press, lats, lons)

    print(u_wave_cube)

    write_data(u_wave_cube, forecast_date, var_name='u_wave')
    write_data(v_wave_cube, forecast_date, var_name='v_wave')
    write_data(z_wave_cube, forecast_date, var_name='z_wave')


    # End of main routine
    ###-----------------------------------------------------------------

if __name__ == '__main__':
    yesterday = datetime.datetime(2019, 10, 29, 0, 0)
    compute_wave_data(yesterday)