from __future__ import print_function
import numpy as np
from context import utils
import nose
import filecmp


class TestRasterIO(object):
    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)

    @classmethod
    def teardown_class(cls):
        print("tearing down class " + cls.__name__)

    def setup(self):
        self.fn_netcdf = '../tests/test_data/precip.nc'


    def teardown(self):
        pass

    def test_init_netcdf(self):
        a = utils.RasterParameterIO(self.fn_netcdf)
        nose.tools.assert_equal(a.fn_base_raster, self.fn_netcdf)

    def test_write_array_as_geotiff(self):
        a = utils.RasterParameterIO(self.fn_netcdf)
        a.write_array_to_geotiff('../tests/test_data/aet_lp_param.tif', np.ones((131, 301)) * 0.5)

class TestRasterDatasetHBV(object):
    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)

    @classmethod
    def teardown_class(cls):
        print("tearing down class " + cls.__name__)

    def setup(self):
        self.fn_netcdf = '../tests/test_data/precip.nc'
        self.fn_json_param = '../tests/test_data/param_files_test.json'
        self.fn_temp_thres = '../tests/test_data/pp_temp_thres.tif'
        self.fn_ddf = '../tests/test_data/ddf.tif'
        self.fn_soil_max_wat = '../tests/test_data/soil_max_wat.tif'
        self.fn_soil_beta = '../tests/test_data/soil_beta.tif'
        self.fn_aet_lp_param = '../tests/test_data/aet_lp_param.tif'

    def teardown(self):
        pass

    def test_init(self):
        a = utils.ModelRasterDatasetHBV(self.fn_netcdf,
                                        self.fn_temp_thres,
                                        self.fn_ddf,
                                        self.fn_soil_max_wat,
                                        self.fn_soil_beta,
                                        self.fn_aet_lp_param)
        nose.tools.assert_equal(a.fn_base_map, self.fn_netcdf)
        nose.tools.assert_equal(a.fn_temp_thres, self.fn_temp_thres)
        nose.tools.assert_equal(a.fn_ddf, self.fn_ddf)
        nose.tools.assert_equal(a.fn_soil_max_wat, self.fn_soil_max_wat)
        nose.tools.assert_equal(a.fn_soil_beta, self.fn_soil_beta)
        nose.tools.assert_equal(a.fn_aet_lp_param, self.fn_aet_lp_param)
        nose.tools.assert_greater(a.base_map.shape[0], 0)
        nose.tools.assert_greater(a.base_map.shape[1], 0)

    def test_write_parameter_input_file(self):
        a = utils.ModelRasterDatasetHBV(self.fn_netcdf,
                                        self.fn_temp_thres,
                                        self.fn_ddf,
                                        self.fn_soil_max_wat,
                                        self.fn_soil_beta,
                                        self.fn_aet_lp_param)
        a.write_parameter_input_file('../tests/test_data/param_files.json')
        filecmp.cmp('../tests/test_data/param_files.json', self.fn_json_param)

    # def test_write_structured_parameter_array(self):
    #     import json
    #     filenamedic = {
    #         'pp_temp_thres': '../tests/test_data/pp_temp_thres.tif',
    #         #'p_base': '../tests/test_data/p_base.tif',
    #         'ddf': '../tests/test_data/ddf.tif',
    #         'soil_max_wat': '../tests/test_data/soil_max_wat.tif',
    #         'soil_beta': '../tests/test_data/soil_beta.tif',
    #         'aet_lp_param': '../tests/test_data/aet_lp_param.tif'
    #     }
    #     with open('test_data/param_files.json', 'w') as src:
    #         json.dump(filenamedic, src)
    #
    #     #utils.write_structured_parameter_array(filenamedic, (131, 301))
    #


   # def test_add_rr_model_parameters_to_shapefile(self):
   #     utils.add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLiteLatLon.shp')

    # def test_add_rr_model_parameters_to_shapefile2(self):
    #     utils.add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLiteLatLon.shp', 'test_data/test')

