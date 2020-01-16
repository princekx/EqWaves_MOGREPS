def dirs(key):
    self = {}
    # analysis query files : Give full paths
    self['queryfiles'] = '/home/h03/hadpx/MJO/Monitoring/eqwaves/analysis/queryfiles'
    #self['mog_forecast_data_dir'] = '/project/local_global/fra7/MJO/forecast/mjo_data_mogg/'
    self['analysis_data_dir'] = '/project/MJO_GCSS/Eqwaves_monitoring/raw_data/analysis/'
    self['analysis_data_out_dir'] = '/project/MJO_GCSS/Eqwaves_monitoring/processed_wave_data/analysis/'
    self['analysis_moose_dir'] = 'moose:/opfc/atm/global/prods/'
    return self[key]