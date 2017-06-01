from __future__ import division
import rasterio as rio
import numpy as np
import json


class ReadRaster(object):
    """
    Base class to read raster datasets, retrieve and keep metadata, and write geotiffs
    """

    def __init__(self, fn_base_raster):
        self.fn_base_raster = fn_base_raster
        # profile and shape are updated with call to _read_base_raster()
        self.profile = None
        self.shape = None
        self.affine = None
        self.array = None
        self._read_base_raster()

    def _read_base_raster(self):
        with rio.open(self.fn_base_raster, 'r') as src:
            self.profile = src.profile
            self.shape = src.shape
            self.affine = src.affine
            self.array = src.read()

    def _write_array_to_geotiff(self, fn_out, np_array):
        """
        Function to write numpy arrays as tiff file using metadata from a template GeoTiff map.

        :param fn: filename of output tiff file, string
        :param array: numpy array to be converted to tiff
        :return: None
        """
        if not isinstance(np_array, np.ndarray) or np_array.shape != self.shape:
            raise ValueError("Shape mismatch!. Template shape is %s. Arrays shape is %s" % (src.shape, array.shape))
        self.profile.update(driver='GTiff', count=1, dtype='float64')
        with rio.open(fn_out, 'w', **self.profile) as dst:
            dst.write(np_array.astype('float64'), 1)


class RasterParameterIO(ReadRaster):
    """
    Class to manipulate raster model parameters
    """
    def __init__(self, fn_base_raster):
        super(RasterParameterIO, self).__init__(fn_base_raster)

    def write_array_to_geotiff(self, fn_out, np_array):
        """
            Wrapper of function defined in parent class to write numpy arrays
             as tiff file using metadata from a template GeoTiff map.

            :param fn_out: filename of output tiff file, string
            :param np_array: numpy array to be converted to tiff
            :return: None
            """
        self._write_array_to_geotiff(fn_out, np_array)


class ModelRaterDatasetHBV(object):

    def __init__(self,
                 fn_base_map,
                 fn_temp_thres=None,
                 fn_ddf=None,
                 fn_soil_max_wat=None,
                 fn_soil_beta=None,
                 fn_aet_lp_param=None):
        """
        Initializes the class holding the HBV's raster parameters. A base map, typically
        a climate input map, is required. The rest of the object is initialized
        with the parameter maps. If the maps do not exist it initializes the maps
        with a default value

        :param fn_temp_thres:
        :param fn_ddf:
        :param fn_soil_max_wat:
        :param fn_soil_beta:
        :param fn_aet_lp_param:
        :return:
        """

        self.fn_base_map = fn_base_map
        self.fn_temp_thres = fn_temp_thres
        self.fn_ddf = fn_ddf
        self.fn_soil_max_wat = fn_soil_max_wat
        self.fn_soil_beta = fn_soil_beta
        self.fn_aet_lp_param = fn_aet_lp_param

        self.base_map = RasterParameterIO(self.fn_base_map)

    def write_parameter_input_file(self, fn_out):
        """
        Writes a json dictionary with the parameter name keys and the path to files
        holding the parameter map. This file is needed to run the hydrovehicle
        :param fn_out: filename to write the json dictionary
        :return: None
        """
        self.filenamedic = {
            'pp_temp_thres': self.fn_temp_thres,
            'ddf': self.fn_ddf,
            'soil_max_wat': self.fn_soil_max_wat,
            'soil_beta': self.fn_soil_beta,
            'aet_lp_param': self.fn_aet_lp_param
        }
        with open(fn_out, 'w') as src:
            json.dump(self.filenamedic, src)





def write_structured_parameter_array(filenames, shape):
    """
    Utility to write a json dictionary with HBV model raster parameters.
    The GeoTiff files must exist

    :param filenames:
    :param shape:
    :return:
    """
    tp = zip(filenames.keys(), ['float32' for i in range(len(filenames))], [shape for i in range(len(filenames))])
    pars = np.empty(len(filenames), dtype=tp)

    for name, items in filenames.items():
        with rio.open(items, 'r') as src:
            pars[name] = src.read()

    #np.zeros(filenames.items(), dtype=[('', '', )])