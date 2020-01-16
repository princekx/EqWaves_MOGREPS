#!/opt/scitools/environments/default/2019_02_27/bin/python
import sys, os
import glob
sys.path.append('/home/h03/hadpx/MJO/Monitoring/eqwaves')
import datetime
from mogreps import data_paths as mogreps_data_paths

def write_html_page(forecast_date, mem=0):
    fc_date_label = forecast_date.strftime("%Y%m%d%H")
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    # hour is forecast hour
    digit3_mem = str('%03d' % mem)
    wave_names = ['Kelvin', 'WMRG', 'R1', 'R2']
    levels = ['850', '200']
    fig_dir = os.path.join(mogreps_data_paths.dirs('plot_dir'))

    wave_gifs_names = {}
    for wave_name in wave_names:
        for level in levels:
            wave_gifs_names['%s_%s' %(wave_name, level)] = os.path.join('./',
                                                                        'anim_wave_%s_%shPa_%s%s%s_%sZ_%s.gif'
                                                                        % (wave_name, int(level),
                                                                           str_year, str_month, str_day,
                                                                           str_hr, digit3_mem))

    print(wave_gifs_names)


    html = """
    <html>
    <head>
    </head>
    <body lang="en-US" dir="ltr">
    <h1>Equatorial waves monitoring</h1>
    <p>Monitoring equatorial waves from MOGREPS forecasts. Only the
    deterministic forecasts are used. The full methodology for the
    computation is given in:<br/>
    Ferrett, S., S.Woolnough, K. Hodges, C. Holloway, J. Methven and G-Y Yang, 2019: Linking Extreme
    Precipitation in Southeast Asia to Equatorial Waves , QJRMS.</p>
    
    <p>
    Geopotential height anomalies are plotted as shaded contours (-10 to 10 with steps of 2 metres). 
    Line contours are for wind convergence for Kelvin and WMRG waves and Vorticity for R1 and R2 waves. 
    Contours are at 2.5 and 5 * 1e6 s-1.
    </p>
    
    <p>Latest available forecast : <font size="+2"> %s </font> </p> \n""" % forecast_date

    html += """<p> <font size="+2" > <a href="./"> All other animations </a>  </font> </p> """

    table = """<table>\n
    <tr>
    <td> Wave </td>
    <td> 850 hPa </td>
    <td> 200 hPa </td>"""

    for wave_name in wave_names:
        table += "<tr>\n"
        table += "<td> %s </td>\n" % wave_name
        wave_gif_name = wave_gifs_names['%s_850' %(wave_name)]
        table += "<td> <a href=\""+wave_gif_name+"\""+">\n"
        table += """<img border="0" src=\""""+wave_gif_name+"\""+""" width="800" >\n"""
        table += "</a></td>\n"

        wave_gif_name = wave_gifs_names['%s_200' % (wave_name)]
        table += "<td> <a href=\"" + wave_gif_name + "\"" + ">\n"
        table += """<img border="0" src=\"""" + wave_gif_name + "\"" + """ width="800" >\n"""
        table += "</a></td>\n"
        table += "</tr>\n"
    table += "</table>"""

    html += table
    html += """
    </body>
    </html>
    """
    print(html)
    file_html = open(os.path.join(fig_dir, 'wave_static_page.html'), 'w')
    file_html.write(html)
    file_html.close()
    print('%s is written. ' %file_html.name)


if __name__ == '__main__':
    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    hours = [0]  # [0, 6, 12, 18]
    for hour in hours:
        yesterday = datetime.datetime(2020, 1, 14, hour, 0)
        #yesterday = datetime.datetime(yesterday.year,
        #                              yesterday.month,
        #                              yesterday.day, hour, 0)
        write_html_page(yesterday)
