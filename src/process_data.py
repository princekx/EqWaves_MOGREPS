import os
import sys
import iris
import numpy as np
import datetime
import iris.experimental.equalise_cubes as iex
from analysis import data_paths as analysis_data_paths
from analysis import retrieve_analysis_data as analysis_retrieval
from mogreps import data_paths as mogreps_data_paths
from src import reg_lat_long_grid

def cube_merger(cubelist, times_list, time_unit):
    # make attributes equal
    iex.equalise_attributes(cubelist)

    for n, uc in enumerate(cubelist):
        if uc.coord('time'):
            uc.remove_coord('time')
        if uc.coord('forecast_reference_time'):
            uc.remove_coord('forecast_reference_time')
        if uc.coord('forecast_period'):
            uc.remove_coord('forecast_period')
        # Add new time coordinate
        uc.add_aux_coord(iris.coords.AuxCoord(times_list[n], long_name='time', units=time_unit))
    return cubelist.merge()

def combine_data_efficient_method(analysis_dates, forecast_date,
                 ntimes_total=360,
                 ntimes_analysis=332,
                 ntimes_forecast=28,
                 mem=0, latmax=24):

    # Create a 1 x 1 degree grid for regridding
    ref_grid = reg_lat_long_grid.create_cube()
    ref_grid = ref_grid.intersection(latitude=(-latmax, latmax))
    print(ref_grid)

    # To set up a time axis
    # arranged such that negative values are for analyses and positive for forecasts
    # with the time unit as
    time_unit = 'hours since %s' %forecast_date
    times_list = np.arange((-ntimes_analysis + 1) * 6, (ntimes_forecast + 1) * 6, 6)

    # date is forecast date
    # hour is forecast hour
    digit3_mem = str('%03d' % mem)
    digit2_mem = str('%02d' % mem)

    fc_times = range(6, 174, 6)  # 6 hourly data

    # Read the analysis data first
    print('Reading analysis data...')
    ana_data_files = []
    for an_dt in analysis_dates:
        str_year, str_month, str_day, str_hr = str(an_dt.year), str('%02d' % an_dt.month), \
                                               str('%02d' % an_dt.day), str('%02d' % an_dt.hour)
        analysis_data_dir = os.path.join(analysis_data_paths.dirs('analysis_data_dir'),
                                         str_year, str_month, str_day)
        print(an_dt.strftime("%Y%m%d %H"))
        # analysis data file
        ana_data_files.append(os.path.join(analysis_data_dir, 'qg%sT000.pp' % (str_hr)))

    ucubes_an = iris.load_cube(ana_data_files, 'x_wind').regrid(ref_grid, iris.analysis.Linear())
    vcubes_an = iris.load_cube(ana_data_files, 'y_wind').regrid(ref_grid, iris.analysis.Linear())
    zcubes_an = iris.load_cube(ana_data_files, 'geopotential_height').regrid(ref_grid, iris.analysis.Linear())

    # Now read the forecasts
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    print('Reading forecast data...')
    fc_data_files = []
    for fcx in fc_times:
        # print('Member is %s' % mem)
        mogreps_data_dir = os.path.join(mogreps_data_paths.dirs('mog_forecast_data_dir'),
                                        str_year, str_month, str_day, digit3_mem)
        fc_data_files.append(os.path.join(mogreps_data_dir, 'qg%sT%s.pp' % (str_hr, str('%03d' % fcx))))

        print(fcx.strftime("%Y%m%d %H"))

    ucubes_fc = iris.load_cube(fc_data_files, 'x_wind').regrid(ref_grid, iris.analysis.Linear())
    vcubes_fc = iris.load_cube(fc_data_files, 'y_wind').regrid(ref_grid, iris.analysis.Linear())
    zcubes_fc = iris.load_cube(fc_data_files, 'geopotential_height').regrid(ref_grid, iris.analysis.Linear())

    ucubes = iris.cube.CubeList([ucubes_an, ucubes_fc])
    vcubes = iris.cube.CubeList([vcubes_an, vcubes_fc])
    zcubes = iris.cube.CubeList([zcubes_an, zcubes_fc])

    # make attributes equal
    iex.equalise_attributes(ucubes)
    iex.equalise_attributes(vcubes)
    iex.equalise_attributes(zcubes)

    # concatenate the data
    ucubes = ucubes.concatenate()
    vcubes = vcubes.concatenate()
    zcubes = zcubes.concatenate()

    '''
    ucubes = cube_merger(ucubes, times_list, time_unit)
    vcubes = cube_merger(vcubes, times_list, time_unit)
    zcubes = cube_merger(zcubes, times_list, time_unit)

    
    if os.path.exists(ana_data_file):
            print(ana_data_file)
            # u wind
            dum = iris.load_cube(ana_data_file, 'x_wind').intersection(latitude=(-latmax, latmax))
            ucubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # v wind
            dum = iris.load_cube(ana_data_file, 'y_wind').intersection(latitude=(-latmax, latmax))
            vcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # geopotential height (m)
            dum = iris.load_cube(ana_data_file, 'geopotential_height').intersection(latitude=(-latmax, latmax))
            zcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))
    
    
    # Now read the forecasts
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    for fcx in fc_times:
        # print('Member is %s' % mem)
        mogreps_data_dir = os.path.join(mogreps_data_paths.dirs('mog_forecast_data_dir'),
                                        str_year, str_month, str_day, digit3_mem)
        fc_data_file = os.path.join(mogreps_data_dir, 'qg%sT%s.pp' % (str_hr, str('%03d' % fcx)))
        if os.path.exists(fc_data_file):
            print(fc_data_file)
            # u wind
            dum = iris.load_cube(fc_data_file, 'x_wind').intersection(latitude=(-latmax, latmax))
            ucubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # v wind
            dum = iris.load_cube(fc_data_file, 'y_wind').intersection(latitude=(-latmax, latmax))
            vcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # geopotential height (m)
            dum = iris.load_cube(fc_data_file, 'geopotential_height').intersection(latitude=(-latmax, latmax))
            zcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))
    
    # Merge the cubes
    try:
        ucubes = cube_merger(ucubes, times_list, time_unit)
    except:
        print('Failed merge ucubes.')
    try:
        vcubes = cube_merger(vcubes, times_list, time_unit)
    except:
        print('Failed merge vcubes.')
    try:
        zcubes = cube_merger(zcubes, times_list, time_unit)
    except:
        print('Failed merge zcubes.')
    '''
    u_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'uwnd_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                       str_day, str_hr, digit3_mem))
    v_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'vwnd_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))
    z_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'zhgt_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))

    iris.save(ucubes, u_file_out, netcdf_format="NETCDF3_CLASSIC")
    iris.save(vcubes, v_file_out, netcdf_format="NETCDF3_CLASSIC")
    iris.save(zcubes, z_file_out, netcdf_format="NETCDF3_CLASSIC")
    print('%s Written.' % z_file_out)


def combine_data_standard_method(analysis_dates, forecast_date,
                                  ntimes_total=360,
                                  ntimes_analysis=332,
                                  ntimes_forecast=28,
                                  mem=0, latmax=24):
    # Create a 1 x 1 degree grid for regridding
    ref_grid = reg_lat_long_grid.create_cube()
    ref_grid = ref_grid.intersection(latitude=(-latmax, latmax))
    print(ref_grid)

    # to set up a time axis
    # arranged such that negative values are for analyses and positive for forecasts
    # with the time unit as
    time_unit = 'hours since %s' % forecast_date
    times_list = np.arange((-ntimes_analysis + 1) * 6, (ntimes_forecast + 1) * 6, 6)

    # date is forecast date
    # hour is forecast hour
    digit3_mem = str('%03d' % mem)
    digit2_mem = str('%02d' % mem)

    fc_times = range(6, 174, 6)  # 6 hourly data

    # Read the analysis data first
    print('Reading analysis data...')
    ucubes = []
    vcubes = []
    zcubes = []
    for an_dt in analysis_dates:
        str_year, str_month, str_day, str_hr = str(an_dt.year), str('%02d' % an_dt.month), \
                                               str('%02d' % an_dt.day), str('%02d' % an_dt.hour)
        analysis_data_dir = os.path.join(analysis_data_paths.dirs('analysis_data_dir'),
                                         str_year, str_month, str_day)
        print(an_dt.strftime("%Y%m%d %H"))
        # analysis data file
        ana_data_file = os.path.join(analysis_data_dir, 'qg%sT000.pp' % (str_hr))

        if os.path.exists(ana_data_file):

            #print(ana_data_file)
            # u wind
            dum = iris.load_cube(ana_data_file, 'x_wind').intersection(latitude=(-latmax, latmax))
            ucubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # v wind
            dum = iris.load_cube(ana_data_file, 'y_wind').intersection(latitude=(-latmax, latmax))
            vcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # geopotential height (m)
            dum = iris.load_cube(ana_data_file, 'geopotential_height').intersection(latitude=(-latmax, latmax))
            zcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))


    # Now read the forecasts
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    print('Reading forecast data...')
    # Now read the forecasts
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    for fcx in fc_times:
        # print('Member is %s' % mem)
        mogreps_data_dir = os.path.join(mogreps_data_paths.dirs('mog_forecast_data_dir'),
                                        str_year, str_month, str_day, digit3_mem)
        fc_data_file = os.path.join(mogreps_data_dir, 'qg%sT%s.pp' % (str_hr, str('%03d' % fcx)))
        if os.path.exists(fc_data_file):
            print(fc_data_file)
            # u wind
            dum = iris.load_cube(fc_data_file, 'x_wind').intersection(latitude=(-latmax, latmax))
            ucubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # v wind
            dum = iris.load_cube(fc_data_file, 'y_wind').intersection(latitude=(-latmax, latmax))
            vcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

            # geopotential height (m)
            dum = iris.load_cube(fc_data_file, 'geopotential_height').intersection(latitude=(-latmax, latmax))
            zcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

    ucubes = cube_merger(iris.cube.CubeList(ucubes), times_list, time_unit)
    vcubes = cube_merger(iris.cube.CubeList(vcubes), times_list, time_unit)
    zcubes = cube_merger(iris.cube.CubeList(zcubes), times_list, time_unit)


    u_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'uwnd_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))
    v_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'vwnd_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))
    z_file_out = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                              'zhgt_%s%s%s_%sZ_%s.nc' % (str_year, str_month,
                                                         str_day, str_hr, digit3_mem))

    iris.save(ucubes, u_file_out, netcdf_format="NETCDF3_CLASSIC")
    iris.save(vcubes, v_file_out, netcdf_format="NETCDF3_CLASSIC")
    iris.save(zcubes, z_file_out, netcdf_format="NETCDF3_CLASSIC")
    print('%s Written.' % z_file_out)

def generate_bias_correction_data(forecast_date,
                                  ndays=30,
                                  ntimes_total=360,
                                  ntimes_analysis=332,
                                  ntimes_forecast=28,
                                  mem=0, latmax=24):
    # Create a 1 x 1 degree grid for regridding
    ref_grid = reg_lat_long_grid.create_cube()
    ref_grid = ref_grid.intersection(latitude=(-latmax, latmax))
    print(ref_grid)
    # Generate a list of data for
    ucubes = []
    vcubes = []
    zcubes = []

    for nd in range(ndays):
        # going back 38 days to allow 30 days + 7 days of forecasts + 1 day of slack
        #analysis_dates = [forecast_date
        #                  - datetime.timedelta(days=(ndays+ntimes_forecast/4+1) - nd)
        #                  + datetime.timedelta(hours=i * 6) for i in range(ntimes_total)]
        analysis_minus38days_date = forecast_date - \
                                    datetime.timedelta(days=(ndays + ntimes_forecast / 4 + 1) - nd)

        print(forecast_date, analysis_minus38days_date)

        # check and retrieve if all data are available
        # 1. Retrieve analysis
        # for each of the 30 previous dates, go back 83 days and retrieve data
        analysis_dates = sorted([analysis_minus38days_date - datetime.timedelta(hours=i * 6)
                                   for i in range(ntimes_total)])

        # Data retrieval for any missing dates
        for an_dt in analysis_dates:
            print(an_dt)
            analysis_retrieval.retrieve_analysis_data(an_dt)

        # Read the analysis data first
        data_files_to_read = []
        for an_dt in analysis_dates:
            str_year, str_month, str_day, str_hr = str(an_dt.year), str('%02d' % an_dt.month), \
                                                   str('%02d' % an_dt.day), str('%02d' % an_dt.hour)
            analysis_data_dir = os.path.join(analysis_data_paths.dirs('analysis_data_dir'),
                                             str_year, str_month, str_day)
            # analysis data file
            ana_data_file = os.path.join(analysis_data_dir, 'qg%sT000.pp' % (str_hr))
            data_files_to_read.append(ana_data_file)


        '''
            if os.path.exists(ana_data_file):
                print(ana_data_file)
                # u wind
                dum = iris.load_cube(ana_data_file, 'x_wind').intersection(latitude=(-latmax, latmax))
                ucubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

                # v wind
                dum = iris.load_cube(ana_data_file, 'y_wind').intersection(latitude=(-latmax, latmax))
                vcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))

                # geopotential height (m)
                dum = iris.load_cube(ana_data_file, 'geopotential_height').intersection(latitude=(-latmax, latmax))
                zcubes.append(dum.regrid(ref_grid, iris.analysis.Linear()))
        '''
        print(ucubes)
        #sys.exit()