import os
import sys
import glob
import datetime 
from analysis import data_paths
from analysis import retrieve_analysis_data as retrieve

def retrieve_analysis_83days(yesterday):
    # to retrieve data for a range of selected dates:
    d1 = yesterday - datetime.timedelta(days=83)  # start date
    d2 = yesterday  # end date
    #d2 = today  # end date
    print(d1, d2)
    delta = d2 - d1         # timedelta
    
    anal_dir = data_paths.dirs('analysis_data_out_dir')
    for i in range(delta.days + 1):
        date = d1 + datetime.timedelta(days=i)
        anal_work_dir = os.path.join(anal_dir,
                                     str(date.year),
                                     str('%02d' % date.month),
                                     str('%02d' % date.day))
        print(anal_work_dir)
        print('Reading data from %s' % anal_work_dir)
        files = glob.glob('%s/*.nc' % anal_work_dir)

        if not files:
            print('Files missing. Will retrieve...')
            retrieve.retrieve_analysis_data(date)


if __name__ == '__main__':
    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    retrieve_analysis_83days(yesterday)