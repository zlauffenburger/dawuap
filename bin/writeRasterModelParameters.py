import argparse
import numpy as np
import utils
import os, sys


def main(argc):
    """This script writes default model raster parameters to a collection of geotiff files.
    """

    fn_json_param = 'param_files_test.json'
    fn_temp_thres = 'pp_temp_thres.tif'
    fn_ddf = 'ddf.tif'
    fn_soil_max_wat = 'soil_max_wat.tif'
    fn_soil_beta = 'soil_beta.tif'
    fn_aet_lp_param = 'aet_lp_param.tif'

    rasterFile = utils.ModelRasterDatasetHBV(argc.basemap,
                                             fn_temp_thres,
                                             fn_ddf,
                                             fn_soil_max_wat,
                                             fn_soil_beta,
                                             fn_aet_lp_param)

    rasterFile.base_map.write_array_to_geotiff(fn_temp_thres,  np.ones(rasterFile.base_map.shape[:2]) * 2.0)
    rasterFile.base_map.write_array_to_geotiff(fn_ddf, np.ones(rasterFile.base_map.shape[:2]) * 0.02)
    rasterFile.base_map.write_array_to_geotiff(fn_soil_max_wat, np.ones(rasterFile.base_map.shape[:2]) * 50.0 )
    rasterFile.base_map.write_array_to_geotiff(fn_soil_beta, np.ones(rasterFile.base_map.shape[:2]) * 0.5)
    rasterFile.base_map.write_array_to_geotiff(fn_aet_lp_param, np.ones(rasterFile.base_map.shape[:2]) * 0.5)

    rasterFile.write_parameter_input_file(fn_json_param)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Writes default raster model parameters to geotiff files and produces json'
                    ' dictionary file with raster parameters ' 
                    ' ',
        epilog="""Default operation is to read a basemap raster file, typically a climate input file, and writes a 
                                    suite of default raster parameter files and the associated json dictionary needed by
                                    hydrovehicle.  

                                     Example: {} ./data/precip.nc""".format(
            os.path.basename(sys.argv[0])))
    parser.add_argument('basemap', type=str, help="Path to a base map with desired projection and raster dimensions. "
                                                  "Any format supported by rasterio.")

    args = parser.parse_args()
    main(args)
