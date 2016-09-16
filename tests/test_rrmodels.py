from unittest import TestCase

from .context import hydroengine
import numpy as np
import nose
import rasterstats as rst
import rasterio as rio
import matplotlib.pyplot as plt
import json
from descartes import PolygonPatch


class Test_hbv(TestCase):

    @classmethod
    def setup_class(klass):
        print "setting up class " + klass.__name__
        Test_hbv.args = {
            'image_res': 30,
            'catch_area': 50000,
            'pp_temp_thres': 2,
            'p_base': 10,
            'ddf': 0.02,
            'soil_max_wat': 50.0,
            'soil_beta': 0.5,
            'aet_lp_param': 0.5,
            'storage_parameter_1': 50,
            'storage_parameter_1': 50,
            'surface_conductance': 50,
            'mid_layer_conductance': 50,
            'lower_layer_conductance': 50,
        }
        Test_hbv.swe_o = np.zeros((20, 30))
        Test_hbv.pond_o = np.zeros_like(Test_hbv.swe_o)
        Test_hbv.sm_o = np.zeros_like(Test_hbv.swe_o)
        Test_hbv.stw1_o = np.zeros_like(Test_hbv.swe_o)
        Test_hbv.stw2_o = np.zeros_like(Test_hbv.swe_o)


    @classmethod
    def teardown_class(klass):
        print "TEAR DOWN class " + klass.__name__

    def setUp(self):
        print "setting up " + self.__class__.__name__
        self.ohbv = hydroengine.hbv(Test_hbv.swe_o, Test_hbv.pond_o, Test_hbv.sm_o, Test_hbv.stw1_o, Test_hbv.stw2_o, **Test_hbv.args)

    def tearDown(self):
        Test_hbv.swe_o = np.zeros((20, 30))
        Test_hbv.pond_o = np.zeros_like(Test_hbv.swe_o)

    #@with_setup(setup, teardown)
    def TestInitialization(self):
        np.testing.assert_array_equal(Test_hbv.swe_o, self.ohbv.swe)
        nose.tools.assert_equals(self.ohbv.p_base, 10)
        nose.tools.assert_equals(self.ohbv.t_thres, 2)
        nose.tools.assert_equals(self.ohbv.ddf, 0.02)

    #@with_setup(setup, teardown)
    def test_snowpack_allswe(self):
        # All rain becomes swe
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        t_max = np.ones_like(Test_hbv.swe_o) * -3
        t_min = np.ones_like(Test_hbv.swe_o) * -4
        self.ohbv.snowpack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(Test_hbv.swe_o)*0.03, 'swe')
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(Test_hbv.swe_o) * 0.0, 'pond')

    #@with_setup(setup, teardown)
    def test_snowpack_allrain(self):
        # All rain becomes ponded water
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        t_max = np.ones_like(Test_hbv.swe_o) * 5
        t_min = np.ones_like(Test_hbv.swe_o) * 3
        self.ohbv.snowpack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(Test_hbv.swe_o) * 0.03)

    #@with_setup(setup, teardown)
    def test_snowpack_mixed(self):
        # 50/50 rain and snow
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        t_max = np.ones_like(Test_hbv.swe_o) * 3
        t_min = np.ones_like(Test_hbv.swe_o) * 1
        self.ohbv.snowpack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(Test_hbv.swe_o)*0.03/2, 'equal swe', verbose=True)
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(Test_hbv.swe_o) * 0.03/2, 'equal rain', verbose=True)

    def test_snowpack_allmelts(self):
        # 50/50 rain and snow but ridiculous melt coefficient to melt all snowpack
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        t_max = np.ones_like(Test_hbv.swe_o) * 4
        t_min = np.ones_like(Test_hbv.swe_o) * 1
        self.ohbv.ddf= 100
        self.ohbv.snowpack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(Test_hbv.swe_o) * 0.0, 'zero swe', verbose=True)
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(Test_hbv.swe_o) * 0.03, 'all rain', verbose=True)

    # Tests for soil_processes module
    def test_soil_processes_allrunsoff(self):
        self.ohbv.sm = np.zeros_like(Test_hbv.swe_o)+50.0
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        pet = np.ones_like(Test_hbv.swe_o) * 0.0
        self.ohbv.pond = precip
        self.ohbv.soil_processes(pet)
        np.testing.assert_array_equal(self.ohbv.ovlnd_flow, precip, 'all rain', verbose=True)

    def test_soil_processes_allinfiltrates(self):
        self.ohbv.sm = np.zeros_like(Test_hbv.swe_o) + 0.0
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        pet = np.ones_like(Test_hbv.swe_o) * 0.0
        self.ohbv.pond = precip
        self.ohbv.soil_processes(pet)
        np.testing.assert_array_equal(self.ohbv.ovlnd_flow, np.zeros_like(Test_hbv.swe_o), 'all rain runs off', verbose=True)
        np.testing.assert_array_equal(self.ohbv.sm, precip, 'all rain infiltrates', verbose=True)

    """
    def test_precipitation_excess(self):
        #precip = rio.open('./tests/test_data/PRCP201301_thematic.tif')
        precip = rio.open('test_data/DEM_64m_1992.tif')
        affine = precip.affine
        runoff = precip.read()
        self.ohbv.runoff = runoff[0, :, :]
        assert(isinstance(self.ohbv.runoff, np.ndarray))

        self.ohbv.precipitation_excess('test_data/WBDHU8_MT.shp', affine)

        #ro_map = json.dumps(self.ohbv.stw1)
        geom = self.ohbv.stw1[0]['geometry']


        BLUE = '#6699cc'
        fig = plt.figure()
        ax = fig.gca()
        ax.add_patch(PolygonPatch(geom, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2))
        ax.axis('scaled')
        plt.show()
     """
    def test_runoff(self):
        self.ohbv.runoff()





