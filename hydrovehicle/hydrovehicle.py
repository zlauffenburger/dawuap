from __future__ import division
import argparse
import pickle
import json
import networkx as nx
import numpy as np
import hydroengine as hyd
import rasterio as rio
import matplotlib.pyplot as plt


# def geojson_to_graphx(in_geojson):
#     reaches = []
#     with open(in_geojson) as src:
#         feats = json.load(src)
#
#     for edges in feats['features']:
#         reaches.append([tuple(edges['geometry']['coordinates'][0]), tuple(edges['geometry']['coordinates'][-1])])
#
#     return reaches, feats
#
# reaches, feats = geojson_to_graphx('/Users/marcomaneta/Documents/DataSandBox/MontanaHydrology/Network.geojson')

def main(argc):

    Test_hbv = {
        'pp_temp_thres': 0,
        'p_base': 5,
        'ddf': 20,
        'soil_max_wat': 500.0,
        'soil_beta': 3,
        'aet_lp_param': 0.5,
    }

    with rio.open(argc.precip) as pp:
        pp_affine = pp.affine
        pp_data = (pp.read()).clip(min=0)

    with rio.open(argc.tmin) as tmin:
        tmin_affine = tmin.affine
        tmin_data = tmin.read() - 273.15

    with rio.open(argc.tmax) as tmax:
        tmax_affine = tmax.affine
        tmax_data = tmax.read() - 273.15

    # retrieve latitude of cells
    lon, lat = tmax_affine * np.indices(tmin_data[0, :, :].shape)

    # initiate hydrologic engine using pickled states from previous run
    if argc.restart:
        swe = pickle.load(open('swe.pickled', 'rb'))
        swe = np.zeros_like(pp_data[0, :, :])
        pond = pickle.load(open('pond.pickled', 'rb'))
        sm = pickle.load(open('sm.pickled', 'rb'))
        soils = pickle.load(open('soils.pickled', 'rb'))
    else:  # initiate hydrologic engine with empty storages
        swe = np.zeros_like(pp_data[0, :, :])
        pond = np.zeros_like(pp_data[0, :, :])
        sm = np.zeros_like(pp_data[0, :, :])
        soils = []


    # retrieve adjacency matrix
    #graph = geojson_to_graphx(argc.network_geojson)




    adj_net = np.array([[0,0,1],[1,0,0],[0,0,0]])
    rr = hyd.HBV(86400, swe, pond, sm, soils, **Test_hbv)
    mc = hyd.routing(adj_net, 86400)




    ro_ts = []
    Q_ts = []
    Q = np.zeros((3))
    qold = np.zeros((3))
    e = np.zeros((3)) + 0.4
    ks = np.zeros((3)) + 864000

    for i in np.arange(pp_data.shape[0]):

        print "Calcuating time step ", i+1
        # Calculate potential evapotranspiration
        pet = hyd.hamon_pe((tmin_data[i, :, :]+tmax_data[i, :, :])*0.5, lat, i)
        runoff = rr.run_time_step(pp_data[i, :, :], tmax_data[i, :, :], tmin_data[i, :, :], pet, argc.basin_shp, affine=tmax_affine)
        runoff[1] = runoff[-1] = 0
        print "runoff", runoff, np.sum(runoff)
        Q = mc.muskingum_routing(Q, ks, e, np.array(runoff), qold)
        qold = np.array(runoff) # np.insert(runoff, 3, 0)
        print "Q", Q, np.sum(Q)
        ro_ts.append(runoff)
        Q_ts.append(Q)

    #plt.imshow(pet)
    #plt.colorbar()
    #plt.show()



    rr.pickle_current_states()


    #plt.imshow(np.clip(et, a_min=0, a_max=10000))
    #plt.show()
    plt.plot(Q_ts)
    #plt.plot(ro_ts)
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hydrologic engine')
    parser.add_argument('precip',  help='file with daily precipitation')
    parser.add_argument('tmin', help='file with minimum daily temperature')
    parser.add_argument('tmax', help='file with maximum daily temperature')

    parser.add_argument('network_geojson', help='file with network topology in geojson format')
    parser.add_argument('basin_shp', help='shapefile with subcatchments for each node')

    parser.add_argument('--restart', dest='restart', action='store_true')

    args = parser.parse_args()
    main(args)
