import numpy as np
import iris
import glob
import os
import sys
from skimage import measure
import itertools
from bokeh.layouts import row, widgetbox, column
from bokeh.models import ColumnDataSource, HoverTool, Select, \
    LinearColorMapper, ColorBar, CheckboxGroup, Button, GeoJSONDataSource,\
    BasicTicker, PrintfTickFormatter
from bokeh.plotting import figure, show, output_file, curdoc
from bokeh.palettes import Spectral11, RdBu

def get_contour_data(cube, levels=None, cmap=Spectral11, line_width=2):
    lons = cube.coord('longitude').points
    lats = cube.coord('latitude').points

    xs = []
    ys = []
    line_colors = []
    line_dashes = []
    line_widths = []

    for level, color in zip(levels, itertools.cycle(cmap)):
        #print(level, color)
        contours = measure.find_contours(cube.data, level)
        for contour in contours:
            xs.append(contour[:, 1] + min(lons))
            ys.append(contour[:, 0] + min(lats))
            line_widths.append(3 + 3 * abs(level) / abs(cube.data.max()))
            if level > 0:
                line_colors.append('crimson')
                line_dashes.append('dashed')
            else:
                line_colors.append('darkgreen')
                line_dashes.append('solid')
                # p.line(x, y, color='grey', line_width=2)
    contour_data = {'contour_xs': xs, 'contour_ys': ys,
                    'contour_line_colors': line_colors,
                    'contour_line_widths': line_widths}
    return contour_data

def read_file_data_times():
    #print(date_select.value)
    z_cube_file = os.path.join(wave_data_dir, 'z_wave_%s_%s.nc' %(wave_select.value, date_select.value))
    if os.path.exists(z_cube_file):
        print(z_cube_file)
        cube = iris.load_cube(z_cube_file)
        time_coord = cube.coord('time')
        cube_times = [time_coord.units.num2date(tp) for tp in time_coord.points]
        fc_times = [t.strftime("%Y%m%d_%HZ") for t in cube_times]
        ddate_split = '_'.join(date_select.value.split('_')[:2])
        index_t = fc_times.index(ddate_split)
        ntimes = len(cube_times)
        indices = [str(i-index_t) for i in range(ntimes)]
        indices_hours = [str(6*(i - index_t)) for i in range(ntimes)]

        # get level info
        pressures = [str(p) for p in cube.coord('pressure').points]
    else:
        indices = ['999']
        indices_hours = ['999']
        fc_times = ['999']
        pressures = ['999']
    return indices, indices_hours, fc_times, pressures

def get_dates_times_mems(wave_data_dir):
    dum_files = glob.glob(os.path.join(wave_data_dir, '*z_wave*Kelvin*.nc'))
    dates_times = np.unique([dum_file.split('Kelvin_')[1].split('.nc')[0] for dum_file in dum_files])
    members = np.unique([dum_file.split('Kelvin_')[1].split('_')[2].split('.')[0] for dum_file in dum_files])
    dates_times = [str(date) for date in dates_times]
    return dates_times, members


def previous_date():
    index = all_dates_times.index(date_select.value)
    date_select.value = all_dates_times[index + 1]

def next_date():
    index = all_dates_times.index(date_select.value)
    date_select.value = all_dates_times[index - 1]

def previous_lead():
    index_l = indices_hours.index(lead_time_select.value)
    lead_time_select.value = indices_hours[index_l - 1]

def next_lead():
    index_l = indices_hours.index(lead_time_select.value)
    lead_time_select.value = indices_hours[index_l + 1]

def open_files(cube_date, cube_wave_name, cube_pressure, cube_lead_time):
    # Z cube
    shade_cube_file = os.path.join(wave_data_dir, 'z_wave_%s_%s.nc' % (cube_wave_name, cube_date))
    shade_cube = iris.load_cube(shade_cube_file).intersection(longitude=(-180, 180))
    #print(shade_cube.coord('time'))
    constraints = iris.Constraint(pressure=float(cube_pressure),
                                  time=lambda cell: cell.point
                                                    == shade_cube.coord('time').units.num2date(
                                      int(cube_lead_time)))
    if cube_wave_name in ['Kelvin', 'WMRG']:
        contour_cube_file = os.path.join(wave_data_dir,
                                      'div_wave_%s_%s.nc' % (cube_wave_name,
                                                             cube_date))
        if os.path.exists(contour_cube_file):
            contour_cube = iris.load_cube(contour_cube_file).intersection(longitude=(-180, 180))
            contour_cube *= -10e6 # to make it convergence

    else:
        contour_cube_file = os.path.join(wave_data_dir,
                                      'vort_wave_%s_%s.nc' % (cube_wave_name,
                                                             cube_date))
        if os.path.exists(contour_cube_file):
            contour_cube = iris.load_cube(contour_cube_file).intersection(longitude=(-180, 180))
            contour_cube *= 10e6

    shade_cube = shade_cube.extract(constraints)
    contour_cube = contour_cube.extract(constraints)

    return shade_cube, contour_cube

def make_plot(data_source, label, shade_levels=[-10, -8, -4, -2, -1, 1, 2, 4, 8, 10],
                       contour_levels=[5]):
    # Make plot
    with open("countries.geo.json", "r") as f:
        countries = GeoJSONDataSource(geojson=f.read())

    x_range = (0, 180)  # could be anything - e.g.(0,1)
    y_range = (-24, 24)
    plot = figure(plot_height=400, plot_width=1100, x_range=x_range, y_range=y_range,
                  tools=["pan, reset, save, wheel_zoom, hover"], title='World Countries',
                  x_axis_label='Longitude', y_axis_label='Latitude', aspect_scale=3)

    color_mapper_z = LinearColorMapper(palette="Spectral11", low=-10, high=10)
    color_bar = ColorBar(color_mapper=color_mapper_z, major_label_text_font_size="15pt",
                         ticker=BasicTicker(desired_num_ticks=len(shade_levels)),
                         label_standoff=6, border_line_color=None, orientation="horizontal",
                         location=(0, 0))

    plot.image('image', x=-180, y=-30,
               dw=360, dh=60,
               source=data_source,
               color_mapper=color_mapper_z)
    plot.patches("xs", "ys", color=None, line_color="grey", source=countries, alpha=0.5)
    #
    plot.multi_line('contour_xs', 'contour_ys', line_color='contour_line_colors',
                    line_width='contour_line_widths', source=data_source, alpha=0.7)
    plot.title.text = label
    plot.add_layout(color_bar, 'below')
    return plot

def update_data(attr, old, new):
    print(attr, old, new)
    print('Update data files')
    print(date_select.value,
          wave_select.value,
          pressure_select.value,
          lead_time_select.value)
    #shade_cube, contour_cube = open_files()
    shade_cube, contour_cube = open_files(date_select.value, wave_select.value,
                                          pressure_select.value, lead_time_select.value)


    src = {'image': [shade_cube.data]}
    label = '%s %s hPa Init Time: %s, Lead:%sH\n' % (wave_select.value,
                                                     pressure_select.value,
                                                     date_select.value,
                                                     lead_time_select.value)
    if selected_wave in ['Kelvin', 'WMRG']:
        label += 'Z (shades, m [-10:2:10]), Convergence (contours, x 10e6 s-1)'
    if selected_wave in ['R1', 'R2']:
        label += 'Z (shades, m [-10:2:10]), Vorticity (contours, x 10e6 s-1) '

    # update title
    plot.title.text = label
    data_source.data.update(src)

    # Update contours
    contour_src = get_contour_data(contour_cube, contour_levels)
    data_source.data.update(contour_src)



wave_data_dir = '/project/MJO_GCSS/Eqwaves_monitoring/processed_wave_data/'
wave_names = ['Kelvin', 'WMRG', 'R1', 'R2']
contour_levels = [2.5, 5, 10]

selected_wave = wave_names[0]
#wave_data_dir = mogreps_data_paths.dirs('mog_forecast_out_dir')
all_dates_times, members = get_dates_times_mems(wave_data_dir)
all_dates_times = all_dates_times[::-1]
selected_date = all_dates_times[0]
print(selected_date, selected_wave)

lead_times = [str(t) for t in range(-96, 174, 6)]
pressures = ['850', '200']

selected_lead_time = '0'
selected_pressure = pressures[0]
shade_cube, contour_cube = open_files(selected_date, selected_wave,
                                      selected_pressure, selected_lead_time)


# set up a drop down menu of available dates
date_select = Select(value=selected_date, title='Forecast:', options=all_dates_times)
wave_select = Select(value=selected_wave, title='Wave name:', options=wave_names)

# read the times in the file
indices, indices_hours, fc_times, pressures = read_file_data_times()

pressure_select = Select(value=pressures[0], title='Pressure level:', options=pressures)
lead_time_select = Select(value='0', title='Lead time:', options=indices_hours)

# Main data source frame
data_source = ColumnDataSource(data={'image': [shade_cube.data]})

# injecting contour information
contour_src = get_contour_data(contour_cube, contour_levels)
data_source.data.update(contour_src)

label = '%s %s hPa Init Time: %s, Lead:%sH\n' % (wave_select.value,
                                                 pressure_select.value,
                                                 date_select.value,
                                                 lead_time_select.value)
if selected_wave in ['Kelvin', 'WMRG']:
    label += 'Z (shades, m [-10:2:10]), Convergence (contours, x 10e6 s-1)'
if selected_wave in ['R1', 'R2']:
    label += 'Z (shades, m [-10:2:10]), Vorticity (contours, x 10e6 s-1) '
print(label)

print(data_source)
# Initial plot
plot = make_plot(data_source, label, contour_levels=contour_levels)

# Buttons for date select
previous_date_button = Button(label="Previous", button_type="success")
next_date_button = Button(label="Next", button_type="success")
previous_date_button.on_click(previous_date)
next_date_button.on_click(next_date)

# buttons for lead time
previous_lead_button = Button(label="Previous", button_type="success")
next_lead_button = Button(label="Next", button_type="success")

previous_lead_button.on_click(previous_lead)
next_lead_button.on_click(next_lead)

# update plot for change of any of the widgets
for w in [date_select, wave_select, pressure_select, lead_time_select]:
    w.on_change('value', update_data)




controls = column(date_select, previous_date_button, next_date_button,
                  wave_select, pressure_select,
                  lead_time_select, previous_lead_button, next_lead_button)


# Set up layouts and add to document
#controls = widgetbox(date_select)

#show(row(controls, plot))
#show(row(plot))

curdoc().add_root(row(controls, plot))
curdoc().title = "Equatorial Waves"
