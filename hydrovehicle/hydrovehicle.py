#from __future__ import print_function
from __future__ import division
import argparse
import utils
import pickle
import json
import numpy as np
import hydroengine as hyd
import rasterio as rio
import matplotlib.pyplot as plt



def main(argc):

    init_date = argc.init_date

    with rio.open(argc.precip) as pp:
        pp_affine = pp.affine
        pp_nodata = int(pp.nodata)
        pp_data = (pp.read()).clip(min=0)

    with rio.open(argc.tmin) as tmin:
        tmin_affine = tmin.affine
        tmin_nodata = int(tmin.nodata)
        tmin_data = tmin.read()
        tmin_data[tmin_data != tmin_nodata] -= 273.15

    with rio.open(argc.tmax) as tmax:
        tmax_affine = tmax.affine
        tmax_nodata = int(tmax.nodata)
        tmax_data = tmax.read()
        tmax_data[tmax_data != tmax_nodata] -= 273.15

    hbv_pars = {}
    with open(argc.params) as json_file:
        pars = json.load(json_file)
        for key, value in pars.items():
            with rio.open(value, 'r') as src:
                hbv_pars[key] = src.read(1)

    # retrieve latitude of cells
    lon, lat = tmax_affine * np.indices(tmin_data[0, :, :].shape)

    # retrieve adjacency matrix
    graph = utils.ParseNetwork(argc.network_file)
    adj_net = graph.conn_matrix
    num_links = len(adj_net.index)

    # initiate hydrologic engine using pickled states from previous run
    if argc.restart:
        swe = pickle.load(open('swe.pickled', 'rb'))
        # swe = np.zeros_like(pp_data[0, :, :])
        pond = pickle.load(open('pond.pickled', 'rb'))
        sm = pickle.load(open('sm.pickled', 'rb'))
        soils = pickle.load(open('soils.pickled', 'rb'))
        Q = pickle.load(open('streamflows.pickled', 'rb'))
    else:  # initiate hydrologic engine with empty storages
        swe = np.zeros_like(pp_data[0, :, :])
        pond = np.zeros_like(pp_data[0, :, :])
        sm = np.zeros_like(pp_data[0, :, :])
        soils = []
        Q = np.zeros(num_links)


    #adj_net = np.array([[0, 0, 1], [1, 0, 0], [0, 0, 0]])
    rr = hyd.HBV(86400, swe, pond, sm, soils, **hbv_pars)
    mc = hyd.Routing(adj_net, 86400)

    ro_ts = []
    Q_ts = []
    qold = np.zeros(num_links)
    e = np.array(graph.get_parameter('e'))
    ks = np.array(graph.get_parameter('ks'))

    for i in np.arange(pp_data.shape[0]):
        print "Calcuating time step ", i + 1
        # Calculate potential evapotranspiration
        pet = hyd.hamon_pe((tmin_data[i, :, :] + tmax_data[i, :, :]) * 0.5, lat, i)
        runoff = rr.run_time_step(pp_data[i, :, :], tmax_data[i, :, :], tmin_data[i, :, :], pet, argc.basin_shp,
                                  affine=tmax_affine, nodata=pp_nodata)
        # runoff[1] = runoff[-1] = 0
        print "runoff ", runoff, np.sum(runoff)
        Q = mc.muskingum_routing(Q, ks, e, np.array(runoff), qold)
        qold = np.array(runoff)  # np.insert(runoff, 3, 0)
        #print "Q", Q, np.sum(Q)
        ro_ts.append(runoff)
        Q_ts.append(Q)

    rr.pickle_current_states()
    # pickle current streamflows
    pickle.dump(Q, open("streamflows.pickled", "wb"))
    print Q_ts
    utils.WriteOutputTimeSeries(adj_net, init_date).write_json(Q_ts)

    # # plt.imshow(np.clip(et, a_min=0, a_max=10000))
    # # plt.show()
    # plt.plot(Q_ts)
    # # plt.plot(ro_ts)
    # plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hydrologic engine')
    parser.add_argument('init_date', help='Simulation start date (mm/dd/yy)')
    parser.add_argument('precip', help='NetCDF file with daily precipitation (mm/day)')
    parser.add_argument('tmin', help='NetCDF file with minimum daily temperature (K)')
    parser.add_argument('tmax', help='file with maximum daily temperature(K)')
    parser.add_argument('params', help='json dictionary with names of parameter files (see documentation)')

    parser.add_argument('network_file', help='stream network, shapefile or geojson format')
    parser.add_argument('basin_shp', help='shapefile with subcatchments for each node')

    parser.add_argument('--restart', dest='restart', action='store_true')

    args = parser.parse_args()
    main(args)
