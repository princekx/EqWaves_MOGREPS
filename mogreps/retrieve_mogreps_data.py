import os
import sys
import datetime
import numpy as np
from mogreps import data_paths


def retrieve_mogreps_forecast_data(forecast_dt, mem=0):
    str_year, str_month, str_day, str_hr = str(forecast_dt.year), str('%02d' % forecast_dt.month), \
                                           str('%02d' % forecast_dt.day), str('%02d' % forecast_dt.hour)

    date_label = '%s%s%s_%s' % (str_year, str_month, str_day, str_hr)
    print('Doing date : %s' % date_label)
    query_files_dir = data_paths.dirs('queryfiles')

    print(query_files_dir)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!!! BUG ALERT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Bug in MOGREPS data archival. It only archives 12 members. 
    # Correct this once the full 36 members are available on moose.
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # This is the correct path for prod run
    moosedir = data_paths.dirs('mog_moose_dir') + '%s%s.pp' % (str_year, str_month)

    # Use this for temporary use for PS39 bug
    # moosedir = data_paths.dirs('mog_moose_dir')# + '%s%s.pp' % (str_year, str_month)

    fc_times = np.arange(6, 174, 6)  # 6 hourly data

    # Ensemble member

    digit3_mem = str('%03d' % mem)
    digit2_mem = str('%02d' % mem)

    print('Member is %s' % mem)
    remote_data_dir = os.path.join(data_paths.dirs('mog_forecast_data_dir'), str_year, str_month, str_day, digit3_mem)

    if not os.path.exists(remote_data_dir):
        print('Making dir: %s' % remote_data_dir)
        os.makedirs(remote_data_dir)
        # os.system('ls %s' %remote_data_dir)

    # os.system('moo ls %s' %file_moose)
    query_file = os.path.join(query_files_dir, 'mogg_query')
    local_query_file1 = os.path.join(query_files_dir, 'local_query1')

    for fc in fc_times:
        print(fc)
        fcx = fc.copy()
        # if fc == 0:
        #    fcx = 3

        fct = str('%03d' % fcx)
        fc_3d = str('%03d' % fc)
        print('Forecast is %s ' % fct)
        filemoose = 'prods_op_mogreps-g_%s%s%s_%s_%s_%s.pp' % (str_year, str_month, str_day, str_hr, digit2_mem, fct)
        # outfile = 'englaa_pd_H%s_T%s.pp' % (hr_2c, fc_3d)
        outfile = 'qg%sT%s.pp' % (str_hr, str('%03d' % fcx))

        file_moose = os.path.join(moosedir, filemoose)
        print(file_moose)

        try:
            # Replace the fctime and filemoose in query file
            replacements = {'fctime': fct, 'filemoose': filemoose}

            with open(query_file) as query_infile, open(local_query_file1, 'w') as query_outfile:
                for line in query_infile:
                    for src, target in replacements.items():
                        line = line.replace(src, target)
                    query_outfile.write(line)

            # do the retrieval
            remote_data_dir_outfile = '%s/%s' % (remote_data_dir, outfile)

            if not os.path.exists(remote_data_dir_outfile):
                command = '/opt/moose-client-wrapper/bin/moo select %s %s %s' % (local_query_file1, \
                                                                                 moosedir, \
                                                                                 remote_data_dir_outfile)
                print(command)
                os.system(command)
            else:
                print('%s found. Skipping retrieval...' % remote_data_dir_outfile)
        except:
            print('%s not returned. Check file on moose' % file_moose)
            sys.exit()
    else:
        print('%s has files. Skipping retrieval...' % remote_data_dir)


if __name__ == '__main__':
    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    # yesterday = datetime.date(2016, 4, 24)

    # Do the analysis yesterday
    retrieve_mogreps_forecast_data(yesterday)
