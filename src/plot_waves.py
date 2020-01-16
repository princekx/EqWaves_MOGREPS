import os
import matplotlib
# use the Agg environment to generate an image rather than outputting to screen
matplotlib.use('Agg')

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from mogreps import data_paths as mogreps_data_paths
from src import compute_div_vort as div_vort


def make_animated_gif(shade_cubes, contour_cubes,
                      shade_levels=[-10, -8, -6, -4, -2, 2, 4, 6, 8, 10],
                      contour_levels=[3],
                      label=None,
                      fig_name=None):
    t = shade_cubes.coord('time').points
    ntime = len(shade_cubes.coord('time').points)

    fig = plt.figure(figsize=(12, 2.5), dpi=100)

    ims = []

    for i in range(ntime):
        im = iplt.contourf(shade_cubes[i], levels=shade_levels, cmap='RdBu', extend='both')
        imc = iplt.contour(contour_cubes[i], levels=contour_levels, cmap='RdBu_r')
        plt.gca().coastlines()
        #plt.colorbar(im, orientation='horizontal')
        plt.tight_layout()

        time_label = label + 'Valid:'+str(shade_cubes[i].coord('time')).split('[')[1].split(']')[0]
        time_label += ' [T=%sh]' % t[i]
        #print(time_label)
        text = plt.gca().text(0, 30, time_label, ha='center')
        add_arts = im.collections + imc.collections
        ims.append(add_arts + [text])

    # generate the animation object
    ani = animation.ArtistAnimation(fig, ims, interval=500, blit=False, repeat_delay=1000)

    # Save the animation to a gif file
    ani.save(fig_name, dpi=100, writer='imagemagick')
    print('Generated the animation %s' % fig_name)


def plot_anim_figures(forecast_date, digit3_mem, div_wave, vort_wave, z_wave, wname):
    fc_date_label = forecast_date.strftime("%Y%m%d%H")
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    # div_wave = div_wave.intersection(longitude=(50, 180))
    # vort_wave = vort_wave.intersection(longitude=(50, 180))
    # z_wave = z_wave.intersection(longitude=(50, 180))

    time_points = z_wave.coord('time').points
    time_coord = div_wave.coord('time')
    cube_dates = [time_coord.units.num2date(tp) for tp in time_coord.points]
    # print(time_coord, forecast_date)
    index = np.where(cube_dates == np.datetime64(forecast_date))[0][0]

    pressure_levels = div_wave.coord('pressure').points
    print(len(cube_dates))

    fig_dir = os.path.join(mogreps_data_paths.dirs('plot_dir'))
    if not os.path.exists(fig_dir):
        print('Making dir: %s' % fig_dir)
        os.makedirs(fig_dir)

    for p, plevel in enumerate(pressure_levels):
        k = 1
        fig_name = os.path.join(fig_dir, 'anim_wave_%s_%shPa_%s%s%s_%sZ_%s.gif' % (wname,
                                                                                   int(plevel),
                                                                                   str_year,
                                                                                   str_month,
                                                                                   str_day,
                                                                                   str_hr,
                                                                                   digit3_mem))
        if wname in ['Kelvin', 'WMRG']:
            # plot convergence
            z_levels = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
            conv_levels = [2.5, 5]

            label = '%s %s hPa Init Time: %s \n' % (wname, plevel, forecast_date)
            label += 'Z (shades, m [-10:2:10]), ' \
                     'Convergence (contours, [2.5, 5] x 10e6 s-1) \n'

            shade_cubes = z_wave.extract(iris.Constraint(pressure=plevel))
            contour_cubes = -1e6 * div_wave.extract(iris.Constraint(pressure=plevel))
            make_animated_gif(shade_cubes, contour_cubes, shade_levels=z_levels,
                              contour_levels=conv_levels, label=label,
                              fig_name=fig_name)
        if wname in ['R1', 'R2']:
            # plot convergence
            z_levels = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
            vort_levels = [2.5, 5]

            label = '%s %s hPa Init Time: %s \n' % (wname, plevel, forecast_date)
            label += 'Z (shades, m [-10:2:10]), ' \
                     'Vorticity (contours, [2.5, 5] x 10e6 s-1) \n'
            shade_cubes = z_wave.extract(iris.Constraint(pressure=plevel))
            contour_cubes = 1e6 * vort_wave.extract(iris.Constraint(pressure=plevel))
            make_animated_gif(shade_cubes, contour_cubes, shade_levels=z_levels,
                              contour_levels=vort_levels, label=label,
                              fig_name=fig_name)


def plot_static_dailymean_figures(forecast_date, div_wave, vort_wave, z_wave, wname):
    fc_date_label = forecast_date.strftime("%Y%m%d%H")

    div_wave = div_wave.intersection(longitude=(50, 180))
    vort_wave = vort_wave.intersection(longitude=(50, 180))
    z_wave = z_wave.intersection(longitude=(50, 180))

    print(div_wave.coord('time').points)
    time_coord = div_wave.coord('time')
    cube_dates = [time_coord.units.num2date(tp) for tp in time_coord.points]
    # print(time_coord, forecast_date)
    index = np.where(cube_dates == np.datetime64(forecast_date))[0][0]
    averaging_dt = 4  # for daily means

    pressure_levels = div_wave.coord('pressure').points
    print(len(cube_dates))

    fig_dir = os.path.join(mogreps_data_paths.dirs('plot_dir'))
    if not os.path.exists(fig_dir):
        print('Making dir: %s' % fig_dir)
        os.makedirs(fig_dir)

    for plevel in pressure_levels:
        k = 1
        for i in range(index + 1, len(cube_dates) - averaging_dt, averaging_dt):
            print(i - index, i - index + averaging_dt)
            # print(cube_dates[i], cube_dates[i+averaging_dt])
            time_label = '%s %s hPa Init Time: %s \n Valid Time:%s to %s \n' % (wname, plevel, forecast_date,
                                                                                cube_dates[i],
                                                                                cube_dates[i + averaging_dt])
            fig_label = '%s_%s_%s_D%s_24hr_mean.png' % (wname, int(plevel), fc_date_label, k)
            print(fig_label)
            div_mean = div_wave[i:i + averaging_dt].collapsed('time', iris.analysis.MEAN)
            div_to_plot = div_mean.extract(iris.Constraint(pressure=plevel))

            vort_mean = vort_wave[i:i + averaging_dt].collapsed('time', iris.analysis.MEAN)
            vort_to_plot = vort_mean.extract(iris.Constraint(pressure=plevel))

            z_mean = z_wave[i:i + averaging_dt].collapsed('time', iris.analysis.MEAN)
            z_to_plot = z_mean.extract(iris.Constraint(pressure=plevel))

            k += 1
            #
            if wname in ['Kelvin', 'WMRG']:
                # plot convergence
                z_levels = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
                conv_levels = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
                time_label += 'Z (shades, m), Convergence (contours, x 10e6 s-1)'
                cb = iplt.contourf(z_to_plot, levels=z_levels, cmap='RdBu', extend='both')

                cb2 = iplt.contour(-1. * div_to_plot * 1e6, levels=conv_levels, cmap='RdBu_r')
                # qplt.contourf(vort_to_plot*1e6)#, levels=np.linspace(-10, 10, 11), cmap='RdBu', extend='both')
                print(time_label)
                plt.gca().coastlines()
                plt.text(111.5, 25, time_label, ha='center')
                plt.colorbar(cb, orientation='horizontal')
                # plt.colorbar(cb2, orientation='horizontal')

            if wname in ['R1', 'R2']:
                # plot convergence
                z_levels = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
                vort_levels = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
                time_label += 'Z (shades, m), Vorticity (contours, x 10e6 s-1)'
                cb = iplt.contourf(z_to_plot, levels=z_levels, cmap='RdBu', extend='both')

                cb2 = iplt.contour(vort_to_plot * 1e6, levels=vort_levels, cmap='RdBu_r')
                # qplt.contourf(vort_to_plot*1e6)#, levels=np.linspace(-10, 10, 11), cmap='RdBu', extend='both')
                print(time_label)
                plt.gca().coastlines()
                plt.text(111.5, 25, time_label, ha='center')
                plt.colorbar(cb, orientation='horizontal')
                # plt.colorbar(cb2, orientation='horizontal')

            plt.savefig(os.path.join(fig_dir, fig_label))
            plt.close()


def plot_waves_maps(forecast_date, mem=0, static_plots=True, animations=True):
    fc_dt = forecast_date
    str_year, str_month, str_day, str_hr = str(fc_dt.year), str('%02d' % fc_dt.month), \
                                           str('%02d' % fc_dt.day), str('%02d' % fc_dt.hour)

    # hour is forecast hour
    digit3_mem = str('%03d' % mem)

    wave_names = np.array(['Kelvin', 'WMRG', 'R1', 'R2'])

    wave_figures_dict = {}

    for wname in wave_names:
        print(wname)
        u_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                        'u_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                        str_month, str_day, str_hr,
                                                                        digit3_mem))

        v_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                        'v_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                        str_month, str_day, str_hr,
                                                                        digit3_mem))

        z_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                        'z_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                        str_month, str_day, str_hr,
                                                                        digit3_mem))

        if os.path.exists(u_wave_data_file):
            print(u_wave_data_file)
            u_wave = iris.load_cube(u_wave_data_file)
        if os.path.exists(v_wave_data_file):
            v_wave = iris.load_cube(v_wave_data_file)
        if os.path.exists(z_wave_data_file):
            z_wave = iris.load_cube(z_wave_data_file)

        # compute vorticity and divergence
        # Divergence
        div_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                          'div_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                            str_month, str_day, str_hr,
                                                                            digit3_mem))

        # Vorticity
        vort_wave_data_file = os.path.join(mogreps_data_paths.dirs('mog_forecast_out_dir'),
                                           'vort_wave_%s_%s%s%s_%sZ_%s.nc' % (wname, str_year,
                                                                              str_month, str_day, str_hr,
                                                                              digit3_mem))
        if not os.path.exists(div_wave_data_file):
            div, vort = div_vort.compute_vort_div(fc_dt)

        if os.path.exists(div_wave_data_file):
            div_wave = iris.load_cube(div_wave_data_file)
        if os.path.exists(vort_wave_data_file):
            vort_wave = iris.load_cube(vort_wave_data_file)

        if static_plots:
            # Plot static figures
            plot_static_dailymean_figures(fc_dt, div_wave, vort_wave, z_wave, wname)

        if animations:
            # plot animations
            plot_anim_figures(fc_dt, digit3_mem, div_wave, vort_wave, z_wave, wname)
