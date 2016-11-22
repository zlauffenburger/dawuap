from __future__ import division
import argparse
import numpy as np
import hydroengine as hyd
import rasterio as rio
import matplotlib.pyplot as plt

def main(argc):

    Test_hbv = {
        'pp_temp_thres': 2,
        'p_base': 10,
        'ddf': 10,
        'soil_max_wat': 50.0,
        'soil_beta': 0.5,
        'aet_lp_param': 0.5,
    }

    with rio.open(argc.precip) as pp:
        pp_affine = pp.affine
        pp_data = pp.read()
        pp.transform

    with rio.open(argc.tmin) as tmin:
        tmin_affine = tmin.affine
        tmin_data = tmin.read() - 273.15

    with rio.open(argc.tmax) as tmax:
        tmax_affine = tmax.affine
        tmax_data = tmax.read() - 273.15

    # retrieve latitude of cells
    lon, lat = tmax_affine * np.indices(tmin_data[0, :, :].shape)

    # initiate hydrologic engine with empty storages
    swe = np.zeros_like(pp_data[0, :, :])
    pond = np.zeros_like(pp_data[0, :, :])
    sm = np.zeros_like(pp_data[0, :, :])
    rr = hyd.HBV(1, swe, pond, sm, **Test_hbv)





    for i in np.arange(1):#(pp_data.shape[0]):

        # Calculate potential evapotranspiration
        pet = hyd.hamon_pe((tmin_data[i, :, :]+tmax_data[i, :, :])*0.5, lat, i)
        runoff = rr.run_time_step(pp_data[i, :, :], tmax_data[i, :, :], tmin_data[i, :, :], pet, argc.basin_shp, affine=tmax_affine)
        print runoff

    #plt.imshow(np.clip(et, a_min=0, a_max=10000))
    #plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hydrologic engine')
    parser.add_argument('precip',  help='file with daily precipitation')
    parser.add_argument('tmin', help='file with minimum daily temperature')
    parser.add_argument('tmax', help='file with maximum daily temperature')

    parser.add_argument('basin_shp', help='shapefile with subcatchments for each node')

    args = parser.parse_args()
    main(args)
