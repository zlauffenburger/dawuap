#from __future__ import print_function
from __future__ import division
import argparse
import datetime
from dateutil.parser import parse
import utils
import pickle
import json
import numpy as np
import hydroengine as hyd
import econengine as econ
from .coupling import HydroEconCoupling
import rasterio as rio
import tqdm


def main(argc):

    init_date = argc.init_date

    pp = utils.RasterParameterIO(argc.precip)
    pp_affine = pp.affine
    pp_nodata = pp.nodata
    pp_data = pp.array.clip(min=0)

    tmin = utils.RasterParameterIO(argc.tmin)
    tmin_affine = tmin.affine
    tmin_nodata = tmin.nodata
    tmin_data = tmin.array
    tmin_data[tmin_data != tmin_nodata] -= 273.15

    tmax = utils.RasterParameterIO(argc.tmax)
    tmax_affine = tmax.affine
    tmax_nodata = tmax.nodata
    tmax_data = tmax.array
    tmax_data[tmax_data != tmax_nodata] -= 273.15

    hbv_pars = {}
    with open(argc.params) as json_file:
        pars = json.load(json_file)
        for key, value in pars.items():
            hbv_pars[key] = utils.RasterParameterIO(value, 1).array

    # Create base raster as template to write tiff outputs
    # base_map = utils.RasterParameterIO()

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

    # Open water user object
    with open('Farms.json') as json_farms:
        farms = json.load(json_farms)

    # Open scenarios object
    with open('Scenarios.json') as json_scenarios:
        scenarios = json.load(json_scenarios)

    # retrieve the list of farms in the json input
    lst_farms = farms['farms']

    """
    Creates the rainfall-runoff model object, the routing object and the model coupling objects 
    """
    rr = hyd.HBV(86400, swe, pond, sm, soils, **hbv_pars)
    mc = hyd.Routing(adj_net, 86400)
    hydecon = HydroEconCoupling(mc, lst_farms)

    # simulates all users with loaded scenarios
    hydecon.simulate_all_users(scenarios)


    ro_ts = []
    Q_ts = []
    qold = np.zeros(num_links)
    e = np.array(graph.get_parameter('e'))
    ks = np.array(graph.get_parameter('ks'))

    with tqdm.tqdm(total=pp_data.shape[0], unit='days') as pbar:

        for i in np.arange(pp_data.shape[0]):
            #print "Calcuating time step ", i + 1

            cur_date = (parse(init_date) + i * datetime.timedelta(seconds=rr.dt)).strftime("%Y%m%d")

            # Calculate potential evapotranspiration
            pet = hyd.hamon_pe((tmin_data[i, :, :] + tmax_data[i, :, :]) * 0.5, lat, i)
            runoff = rr.run_time_step(pp_data[i, :, :], tmax_data[i, :, :], tmin_data[i, :, :], pet, argc.basin_shp,
                                      affine=tmax_affine, nodata=pp_nodata)
            # runoff[1] = runoff[-1] = 0

            Q = mc.muskingum_routing(Q, ks, e, np.array(runoff), qold)
            qold = np.array(runoff)  # np.insert(runoff, 3, 0)

            ro_ts.append(runoff)
            Q_ts.append(Q)

            # write to drive the states for the current time step

            latlon = np.array((lat.ravel(), lon.ravel()))
            rr.write_current_states(cur_date, ".tif", pp.write_array_to_geotiff)
            pbar.update()

        rr.pickle_current_states()
        # pickle current streamflows
        pickle.dump(Q, open("streamflows.pickled", "wb"))

        utils.WriteOutputTimeSeries(adj_net, init_date).write_json(Q_ts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hydrologic engine')
    parser.add_argument('init_date', help='Simulation start date (mm/dd/yy)')
    parser.add_argument('precip', help='NetCDF file with daily precipitation (mm/day)')
    parser.add_argument('tmin', help='NetCDF file with minimum daily temperature (K)')
    parser.add_argument('tmax', help='file with maximum daily temperature(K)')
    parser.add_argument('params', help='json dictionary with names of parameter files (see documentation)')

    parser.add_argument('network_file', help='stream network, shapefile or geojson format')
    parser.add_argument('basin_shp', help='shapefile with subcatchments for each node')

    #parser.add_argument('--farms', type=str,  dest='restart', action='store_true')

    parser.add_argument('--restart', dest='restart', action='store_true')



    args = parser.parse_args()
    main(args)
