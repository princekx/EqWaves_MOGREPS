import iris

def main():
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Analysis data for 2019/09/19/18Z is retrieved but is corrupt. So we are
    # replacing that with average of 2019/09/19/12Z and 2019/09/20/00Z
    data1 = iris.load('/project/MJO_GCSS/Eqwaves_monitoring/raw_data/analysis/2019/09/19/qg12T000.pp')
    data2 = iris.load('/project/MJO_GCSS/Eqwaves_monitoring/raw_data/analysis/2019/09/20/qg00T000.pp')

    new_data =[]
    for i in range(len(data1)):
        data_mean = data1[i].copy()
        data_mean.data = 0.5 * (data1[i].data + data2[i].data)
        data_mean.long_name = data1[i].name()
        data_mean.attributes.name = data1[i].name()
        data_mean.remove_coord('time')

        tc = data1[i].coord('time')
        time_point = tc.points[0]
        new_tc = iris.coords.DimCoord([time_point + 6], long_name='time', units=tc.units)
        data_mean.add_aux_coord(new_tc)
        print(data_mean.coord('time'))
        new_data.append(data_mean)

    data_out_file = '/project/MJO_GCSS/Eqwaves_monitoring/raw_data/analysis/2019/09/19/qg18T000.pp'
    iris.save(new_data, data_out_file)
    print(data_out_file)

if __name__ == '__main__':
    main()