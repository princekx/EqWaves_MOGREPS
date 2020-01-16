def dirs(key):
    self = {}
    self['queryfiles'] = '/net/home/h03/hadpx/MJO/Monitoring/eqwaves/mogreps/queryfiles'
    self['mog_moose_dir'] = 'moose:/opfc/atm/mogreps-g/prods/'
    self['mog_forecast_data_dir'] = '/scratch/hadpx/Eqwaves_monitoring/raw_data/mogreps'
    self['mog_forecast_out_dir'] = '/scratch/hadpx/Eqwaves_monitoring/processed_wave_data/mogreps'
    #self['mog_forecast_ascii_dir'] = '/project/MJO_GCSS/Eqwaves_monitoring/processed_data/mogreps/mjo_data_ascii'
    self['analy_out_dir'] = '/scratch/hadpx/Eqwaves_monitoring/processed_data/analysis'
    self['analy_in_dir'] = '/scratch/hadpx/Eqwaves_monitoring/processed_data/analysis'
    self['plot_dir'] = '/scratch/hadpx/Eqwaves_monitoring/processed_wave_data/mogreps/plots/'
    return self[key]
