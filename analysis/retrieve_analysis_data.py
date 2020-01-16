import os
import sys
import datetime
from analysis import data_paths


def retrieve_analysis_data(forecast_dt):
    str_year, str_month, str_day, str_hr = str(forecast_dt.year), str('%02d' % forecast_dt.month), \
                                           str('%02d' % forecast_dt.day), str('%02d' % forecast_dt.hour)
                                                                            
    date_label = '%s%s%s_%s' % (str_year, str_month, str_day, str_hr)
    print(('Doing date : %s' % date_label))
    query_files_dir = data_paths.dirs('queryfiles')

    #fc_times = np.arange(0, 168, 12)
    fcx = 0 # analysis data

    moosedir = os.path.join(data_paths.dirs('analysis_moose_dir'), '%s.pp' % (str_year))

    remote_data_dir = os.path.join(data_paths.dirs('analysis_data_dir'), str_year, str_month, str_day)
    if not os.path.exists(remote_data_dir):
        os.makedirs(remote_data_dir)
    # for fc in fc_times:
    #     fcx = fc.copy()
    #     if fc == 0:
    #         fcx = 3
    #
    fct = str('%03d' % fcx)
    print(('Forecast is %s ' % fct))

    # Filename on moose:
    # .calc.pp extention not needed for files after PS41 on 2018/09/25
    if datetime.datetime(int(str_year), int(str_month), int(str_day)) < datetime.datetime(2018, 9, 25):
        filemoose = 'prods_op_gl-mn_%s%s%s_%s_%s.calc.pp' % (str_year, str_month, str_day, str_hr, fct)
    else:
        filemoose = 'prods_op_gl-mn_%s%s%s_%s_%s.pp' % (str_year, str_month, str_day, str_hr, fct)

    # Data to be stored locally as:
    outfile = 'qg%sT%s.pp' %(str_hr, str('%03d' % fcx))

    file_moose = os.path.join(moosedir, filemoose)
    #print(file_moose)

    try:
        #os.system('moo ls %s' %file_moose)
        query_file = os.path.join(query_files_dir, 'all_vars')
        local_query_file1 = os.path.join(query_files_dir, 'localquery')
        #print(local_query_file1)
        # Replace the fctime and filemoose in query file
        replacements = {'fctime':fct, 'filemoose':filemoose}

        with open(query_file) as query_infile, open(local_query_file1, 'w') as query_outfile:
            for line in query_infile:
                for src, target in replacements.items():
                    line = line.replace(src, target)
                query_outfile.write(line)

        print(local_query_file1)

        # do the retrieval
        # To check whether file is present and is not empty,
        outfile_path = os.path.join(remote_data_dir, outfile)
        if os.path.exists(outfile_path) and os.path.getsize(outfile_path) > 0:
            print(('%s exist. Skipping retrieval.' % outfile_path))
        else:
            print(outfile)
            command = '/opt/moose-client-wrapper/bin/moo select %s %s %s' % (local_query_file1, moosedir, outfile_path)
            print(command)
            # Execute command
            os.system(command)
    except:
        print(('%s not returned. Check file on moose' % file_moose))
        sys.exit()

             
                
            
if __name__ == '__main__':
    today = datetime.date.today()
    #yesterday = today - datetime.timedelta(days=1)
    #yesterday = today

    # Do the analysis
    d1 = datetime.date(2019, 1, 1)  # start date
    d2 = datetime.date(2019, 1, 1)  # end date
    delta = d2 - d1  # timedelta

    for i in range(delta.days + 1):
        yesterday = d1 + datetime.timedelta(days=i)

        retrieve_analysis_data(yesterday)

