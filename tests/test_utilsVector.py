from __future__ import print_function
import numpy as np
import pandas as pd
from context import utils
import pickle
import nose
import filecmp
import os
import fiona

class TestParseNetwork(object):

    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)


    @classmethod
    def teardown_class(cls):
        print ("tearing down class " + cls.__name__)

    def setup(self):
        self.network_geojson = '../tests/test_data/mt_network.geojson'
        self.network_shp = '../tests/test_data/mt_network.shp'
        self.watsheds_shp = '../tests/test_data/mt_subsheds.shp'
        self.watsheds_geojson = '../tests/test_data/mt_subsheds.geojson'

    def teardown(self):
        pass

    def test_init_geojson(self):
        a = utils.ParseNetwork(self.network_geojson)
        nose.tools.assert_equal(a.fn_vector, self.network_geojson)

    def test_init_shp(self):
        a = utils.ParseNetwork(self.network_shp)
        nose.tools.assert_equal(a.fn_vector, self.network_shp)

    def test_calc_connectivity_matrix(self):
        # Load pre-saved connectivity matrix
        with open('../tests/test_data/conn_matrix.pickle', 'r') as f:
            conmat = pickle.load(f)

        a = utils.ParseNetwork(self.network_shp)
        pd.DataFrame.equals(conmat, a.conn_matrix)

    def test_get_parameter(self):
        a = utils.ParseNetwork(self.network_shp)
        area = a.get_parameter('AreaC')
        print(str(area))
        isinstance(area, list)


class TestParameterIO(object):
    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)

    @classmethod
    def teardown_class(cls):
        print("tearing down class " + cls.__name__)

    def setup(self):
        self.network_geojson = '../tests/test_data/mt_network.geojson'
        self.network_shp = '../tests/test_data/mt_network.shp'
        self.watsheds_shp = '../tests/test_data/mt_subsheds.shp'
        self.watsheds_geojson = '../tests/test_data/mt_subsheds.geojson'

    def teardown(self):
        pass

    def test_init_geojson(self):
        a = utils.VectorParameterIO(self.network_geojson)
        nose.tools.assert_equal(a.fn_vector, self.network_geojson)

    def test_init_shp(self):
        a = utils.VectorParameterIO(self.network_shp)
        nose.tools.assert_equal(a.fn_vector, self.network_shp)

    def test_write_dataset_shp(self):
        outfile = r'../tests/test_data/output.shp'
        a = utils.VectorParameterIO(self.network_shp)
        a.write_dataset(outfile)
        assert filecmp.cmp(outfile, a.fn_vector, shallow=False)

    def test_write_dataset_geojson(self):
        outfile = r'../tests/test_data/output.geojson'
        a = utils.VectorParameterIO(self.network_geojson)
        a._write_fiona_object(outfile)
        assert filecmp.cmp(outfile, a.fn_vector, shallow=False)


class TestModelVectorDatasets(object):
    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)

    @classmethod
    def teardown_class(cls):
        print("tearing down class " + cls.__name__)

    def setup(self):
        self.network_geojson = '../tests/test_data/mt_network.geojson'
        self.network_shp = '../tests/test_data/mt_network.shp'
        self.watsheds_shp = '../tests/test_data/mt_subsheds.shp'
        self.watsheds_geojson = '../tests/test_data/mt_subsheds.geojson'

    def teardown(self):
        pass

    def test_write_muskingum_parameters_shp(self):
        outfile = r'../tests/test_data/mt_network_par.shp'
        params = []
        #  create parameter file for testing
        with fiona.open(self.network_shp, 'r') as src:
            for item in src:
                prop = item['properties']
                d = {'ARCID': prop['ARCID'],
                     'e': 0.4,
                     'ks': 86400}
                params.append(d)

        a = utils.ModelVectorDatasets(fn_network=self.network_shp)
        a.write_muskingum_parameters(outfile, params)
        # check if the parameters are in the output files
        with fiona.open(outfile, 'r') as src:
            for item in src:
                prop = item['properties']
                nose.tools.assert_almost_equals(prop['e'], 0.4)

    def test_write_muskingum_parameters_geojson(self):
        outfile = r'../tests/test_data/mt_network_par.geojson'
        params = []
        #  create parameter file for testing
        with fiona.open(self.network_geojson, 'r') as src:
            for item in src:
                prop = item['properties']
                d = {'ARCID': prop['ARCID'],
                     'e': 0.4,
                     'ks': 86400}
                params.append(d)

        a = utils.ModelVectorDatasets(fn_network=self.network_geojson)
        a.write_muskingum_parameters(outfile, params)
        # check if the parameters are in the output files
        with fiona.open(outfile, 'r') as src:
            for item in src:
                prop = item['properties']
                nose.tools.assert_almost_equals(prop['e'], 0.4)

    def test_write_hbv_parameters_shp_default(self):
        outfile = r'../tests/test_data/output2.shp'
        params = []
        a = utils.ModelVectorDatasets(fn_subsheds=self.watsheds_shp)
        a.write_hvb_parameters(outfile, params)
        #check if the parameters are in the output files
        with fiona.open(outfile, 'r') as src:
            for item in src:
                prop = item['properties']
                nose.tools.assert_almost_equals(prop['hbv_ck0'], 10)

    def test_write_hbv_parameters_shp(self):
        outfile = r'../tests/test_data/output2.shp'
        params = []
        #  create parameter file for testing
        with fiona.open(self.watsheds_shp, 'r') as src:
            for item in src:
                prop = item['properties']
                d = {'GRIDCODE': prop['GRIDCODE'],
                     'hbv_ck0': 12,
                     'hbv_ck1': 50,
                     'hbv_ck2': 10000,
                     'hbv_hl1': 50,
                     'hbv_perc': 50,
                     'hbv_pbase': 5}
                params.append(d)

        a = utils.ModelVectorDatasets(fn_subsheds=self.watsheds_shp)
        a.write_hvb_parameters(outfile, params)
        # check if the parameters are in the output files
        with fiona.open(outfile, 'r') as src:
            for item in src:
                prop = item['properties']
                nose.tools.assert_almost_equals(prop['hbv_ck0'], 12)


    # @nose.tools.raises(Exception)
    # def test_write_hbv_parameters_exception(self):
    #     outfile = r'../tests/test_data/output2.shp'
    #     params = []
    #     a = utils.ModelVectorDatasets(fn_network=self.watsheds_shp)
    #     a.write_hvb_parameters(outfile, params)

# Test_hbv = np.zeros(6, dtype='3int8, float32, (2,3)float64')
# Test_hbv = {
#     'pp_temp_thres': 0,
#     'p_base': 5,
#     'ddf': 20,
#     'soil_max_wat': 500.0,
#     'soil_beta': 3,
#     'aet_lp_param': 0.5,
# }





        #utils.write_structured_parameter_array(filenamedic, (131, 301))

    # def test_write_array_as_tiff(self):
    #     utils.write_array_as_tiff('test_data/aet_lp_param.tif', 'test_data/precip.nc', np.ones((131, 301))*0.5)

   # def test_add_rr_model_parameters_to_shapefile(self):
   #     utils.add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLiteLatLon.shp')

    # def test_add_rr_model_parameters_to_shapefile2(self):
    #     utils.add_rr_model_parameters_to_shapefile('test_data/HUC8_NetworkLiteLatLon.shp', 'test_data/test')

