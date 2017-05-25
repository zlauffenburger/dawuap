from __future__ import print_function
from unittest import TestCase
import numpy as np
from hydroutils import utils



class TestParseNetwork(object):

    @classmethod
    def setup_class(klass):
        print "SETUP!"

    @classmethod
    def teardown_class(klass):
        print "TEAR DOWN!"

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_init(self):
        a = utils.parse_network('test_data/mt_network.geojson')
        nose.tools.assert_equal('test_data/mt_network.geojson', a.fn_vector)

    def test_parse_network(self):
        pass


# Test_hbv = np.zeros(6, dtype='3int8, float32, (2,3)float64')
# Test_hbv = {
#     'pp_temp_thres': 0,
#     'p_base': 5,
#     'ddf': 20,
#     'soil_max_wat': 500.0,
#     'soil_beta': 3,
#     'aet_lp_param': 0.5,
# }


class TestUtils(TestCase):
    def test_write_structured_parameter_array(self):
        import json
        filenamedic = {
            'pp_temp_thres': '../tests/test_data/pp_temp_thres.tif',
            #'p_base': '../tests/test_data/p_base.tif',
            'ddf': '../tests/test_data/ddf.tif',
            'soil_max_wat': '../tests/test_data/soil_max_wat.tif',
            'soil_beta': '../tests/test_data/soil_beta.tif',
            'aet_lp_param': '../tests/test_data/aet_lp_param.tif'
        }
        with open('test_data/param_files.json', 'w') as src:
            json.dump(filenamedic, src)

        #utils.write_structured_parameter_array(filenamedic, (131, 301))

    def test_write_array_as_tiff(self):
        utils.write_array_as_tiff('test_data/aet_lp_param.tif', 'test_data/precip.nc', np.ones((131, 301))*0.5)

   # def test_add_rr_model_parameters_to_shapefile(self):
   #     utils.add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLiteLatLon.shp')

    def test_add_rr_model_parameters_to_shapefile2(self):
        utils.add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLiteLatLon.shp', 'test_data/test')

