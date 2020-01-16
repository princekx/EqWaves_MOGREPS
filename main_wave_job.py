#!/opt/scitools/environments/default/2019_02_27/bin/python
import sys, os
import glob
sys.path.append('/home/h03/hadpx/MJO/Monitoring/eqwaves')
import datetime
from analysis import retrieve_analysis_data as analysis_retrieval
from mogreps import retrieve_mogreps_data as mogreps_retrieval
from mogreps import data_paths as mogreps_data_paths
from src import process_data
from src import realtime_waves_compute as wave
from src import compute_div_vort as div_vort
from src import plot_waves
from src import generate_html

def check_if_done(fc_dt):
    '''
    Checks if the calculation for the date has already been done.

    :param fc_dt: forecast date as datetime object
    :return: logical
    '''
    str_year, str_month, str_day, str_hr = str(fc_dt.year), \
                                           str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), \
                                           str('%02d' % fc_dt.hour)
    wave_data_files = glob.glob(os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                             '*_wave_*_%s%s%s_%sZ_*.nc' % (str_year, str_month,
                                                                           str_day, str_hr)))
    if len(wave_data_files) < 20:
        return False
    else:
        return True


def do_main(forecast_dt, ntimes_total=360,
            ntimes_analysis=332,
            ntimes_forecast=28):
    '''
    Total combined data length = 90 days (360 points i.e, 83 x 4 analyses + 7 x 4 forecasts)
    that is 332 points of analyses and 28 forecast points @ 6hrly output
    start date minus 332 hours ( = 83 days x 4 samples a day)

    :param forecast_dt: forecast date as datetime object
    :param ntimes_total: total number of time points in the data
                        ntimes_total = ntimes_analysis + ntimes_forecast
    :param ntimes_analysis: total number of analysis data attached
    :param ntimes_forecast: number of forecast points
    :return:
    '''

    # Check if the computation is already done
    if not check_if_done(forecast_dt):

        print('Computations are not/partially done for %s. Redoing...' % forecast_dt)
        # 1. Retrieve analysis
        analysis_dts = sorted([forecast_dt - datetime.timedelta(hours=i * 6)
                               for i in range(ntimes_analysis)])
        for fdt in analysis_dts:
            analysis_retrieval.retrieve_analysis_data(fdt)

        # 2. Retrieve MOGREPS
        # Only forecast_dt is needed here as it will by default download
        # all 7 days of forecast data
        mogreps_retrieval.retrieve_mogreps_forecast_data(forecast_dt)

        # 3. Combine data
        process_data.combine_data_standard_method(analysis_dts, forecast_dt)
        # TODO Below is an efficient method, but need more work to get it working
        # process_data.combine_data_efficient_method(analysis_dts, forecast_dt)

        # 4. Compute waves
        wave.compute_wave_data(forecast_dt)

        # 5. TODO: Post process
        # process_data.generate_bias_correction_data(forecast_dt, ntimes_total=360)

        # 6. Compute vorticity and divergence
        div_vort.compute_vort_div(forecast_dt)
    else:
        print('Computations are done for %s. Just skipping to plots...' % forecast_dt)

    # 6. Plot wave data
    # This generates static or animated plots.
    plot_waves.plot_waves_maps(forecast_dt, static_plots=False, animations=True)

    # 7. Generate an HTML file
    generate_html.write_html_page(forecast_dt)

if __name__ == '__main__':
    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    hours = [0, 6, 12, 18]
    for hour in hours:
        #yesterday = datetime.datetime(2019, 12, 18, hour, 0)
        yesterday = datetime.datetime(yesterday.year,
                                      yesterday.month,
                                      yesterday.day, hour, 0)
        do_main(yesterday)

    ''' 
    # Do the analysis
    d1 = datetime.date(2018, 12, 22)  # start date
    d2 = datetime.date(2019, 01, 06)  # end date
    delta = d2 - d1         # timedelta

    for i in range(delta.days + 1):
        date = d1 + datetime.timedelta(days=i)
        do_main(date)
	'''
