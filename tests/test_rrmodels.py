from __future__ import print_function
from context import hydroengine
import numpy as np
import cPickle as pickle
import nose
import rasterio as rio
#import matplotlib.pyplot as plt
import json


class TestHBV(object):

    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)
        cls.args = {
            'image_res': 30,
            'catch_area': 50000,
            'pp_temp_thres': 2,
            'p_base': 10,
            'ddf': 0.02,
            'soil_max_wat': np.ones((20, 30))*50.0,
            'soil_beta': 0.5,
            'aet_lp_param': 0.5,
            'storage_parameter_1': 50,
            'storage_parameter_1': 50,
            'surface_conductance': 50,
            'mid_layer_conductance': 50,
            'lower_layer_conductance': 50,
        }
        cls.swe_o = np.zeros((20, 30))
        cls.pond_o = np.zeros_like(TestHBV.swe_o)
        cls.sm_o = np.zeros_like(TestHBV.swe_o)
        cls.aet_o = np.zeros_like(TestHBV.swe_o)


    @classmethod
    def teardown_class(klass):
        print("TEAR DOWN class " + klass.__name__)

    def setUp(self):
        print("setting up " + self.__class__.__name__)
        TestHBV.ohbv = hydroengine.HBV(86400, TestHBV.swe_o, TestHBV.pond_o, TestHBV.sm_o, **TestHBV.args)
        TestHBV.ohbv.aet = self.aet_o

    def tearDown(self):
        TestHBV.swe_o = np.zeros((20, 30))
        TestHBV.pond_o = np.zeros_like(TestHBV.swe_o)

    #@with_setup(setup, teardown)
    def TestInitialization(self):
        np.testing.assert_array_equal(TestHBV.swe_o, self.ohbv.swe)
        nose.tools.assert_equals(self.ohbv.t_thres, 2)
        nose.tools.assert_equals(self.ohbv.ddf, 0.02)

    #@with_setup(setup, teardown)
    def test_snowpack_allswe(self):
        # All rain becomes swe
        precip = np.ones_like(TestHBV.swe_o) * 0.03
        t_max = np.ones_like(TestHBV.swe_o) * -3
        t_min = np.ones_like(TestHBV.swe_o) * -4
        self.ohbv.snow_pack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(TestHBV.swe_o)*0.03, 'swe')
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(TestHBV.swe_o) * 0.0, 'pond')

    #@with_setup(setup, teardown)
    def test_snowpack_allrain(self):
        # All rain becomes ponded water
        precip = np.ones_like(TestHBV.swe_o) * 0.03
        t_max = np.ones_like(TestHBV.swe_o) * 5
        t_min = np.ones_like(TestHBV.swe_o) * 3
        self.ohbv.snow_pack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(TestHBV.swe_o) * 0.03)

    #@with_setup(setup, teardown)
    def test_snowpack_mixed(self):
        # 50/50 rain and snow
        precip = np.ones_like(TestHBV.swe_o) * 0.03
        t_max = np.ones_like(TestHBV.swe_o) * 3
        t_min = np.ones_like(TestHBV.swe_o) * 1
        self.ohbv.snow_pack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(TestHBV.swe_o)*0.03/2, 'equal swe', verbose=True)
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(TestHBV.swe_o) * 0.03/2, 'equal rain', verbose=True)

    def test_snowpack_allmelts(self):
        # 50/50 rain and snow but ridiculous melt coefficient to melt all snowpack
        precip = np.ones_like(TestHBV.swe_o) * 0.03
        t_max = np.ones_like(TestHBV.swe_o) * 4
        t_min = np.ones_like(TestHBV.swe_o) * 1
        self.ohbv.ddf= 100
        self.ohbv.snow_pack(precip, t_max, t_min)
        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(TestHBV.swe_o) * 0.0, 'zero swe', verbose=True)
        np.testing.assert_array_equal(self.ohbv.pond, np.ones_like(TestHBV.swe_o) * 0.03, 'all rain', verbose=True)

    # Tests for soil_processes module
    def test_soil_processes_allrunsoff(self):
        self.ohbv.sm = np.zeros_like(TestHBV.swe_o)+50.0
        precip = np.ones_like(TestHBV.swe_o) * 0.03
        pet = np.ones_like(TestHBV.swe_o) * 0.0
        self.ohbv.pond = precip
        self.ohbv.soil_processes(pet)
        np.testing.assert_array_equal(self.ohbv.ovlnd_flow, precip, 'all rain', verbose=True)

    def test_soil_processes_allinfiltrates(self):
        self.ohbv.sm = np.zeros_like(TestHBV.swe_o) + 0.0
        precip = np.ones_like(TestHBV.swe_o) * 0.03
        pet = np.ones_like(TestHBV.swe_o) * 0.0
        self.ohbv.pond = precip
        self.ohbv.soil_processes(pet)
        np.testing.assert_array_equal(self.ohbv.ovlnd_flow, np.zeros_like(TestHBV.swe_o), 'all rain runs off', verbose=True)
        np.testing.assert_array_equal(self.ohbv.sm, precip, 'all rain infiltrates', verbose=True)


    def test_precipitation_excess_nosoilfile(self):

        precip = rio.open('test_data/PRCP201301_thematic.tif')
        #precip = rio.open('test_data/DEM_64m_1992.tif')
        affine = precip.transform
        runoff = precip.read()
        self.ohbv.ovlnd_flow = runoff[0, :, :]
        assert(isinstance(self.ohbv.ovlnd_flow, np.ndarray))

        self.ohbv.precipitation_excess('test_data/HUC8_NetworkLite.shp', affine)
        print(self.ohbv.soils[0][1].Qall)
        print(self.ohbv.runoff)

    def test_precipitation_excess_soilfile(self):
        precip = rio.open('test_data/PRCP201301_thematic.tif')
        # precip = rio.open('test_data/DEM_64m_1992.tif')
        affine = precip.transform
        runoff = precip.read()
        self.ohbv.ovlnd_flow = runoff[0, :, :]
        assert (isinstance(self.ohbv.ovlnd_flow, np.ndarray))

        self.ohbv.soils = pickle.load(open('soils.pickled', "rb"))
        self.ohbv.precipitation_excess('test_data/HUC8_NetworkLite.shp', affine)
        print(self.ohbv.soils[0][1].Qall)
        print(self.ohbv.runoff)
        # #ro_map = json.dumps(self.ohbv.stw1)
        # geom = self.ohbv.stw1[0]['geometry']
        # geom2 = self.ohbv.stw1[1]['geometry']
        # geom3 = self.ohbv.stw1[2]['geometry']
        #
        #
        #
        # BLUE = '#6699cc'
        # fig = plt.figure()
        # ax = fig.gca()
        # ax.add_patch(PolygonPatch(geom, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2))
        # ax.add_patch(PolygonPatch(geom2, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2))
        # ax.add_patch(PolygonPatch(geom3, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2))
        # ax.axis('scaled')
        # plt.show()

    def test_pickle_current_states(self):
        self.ohbv.pickle_current_states()
        sm = pickle.load(open("sm.pickled","rb"))
        swe = pickle.load(open("swe.pickled", "rb"))
        pond = pickle.load(open("pond.pickled", "rb"))
        nose.tools.assert_equals(self.ohbv.sm.all(), sm.all())
        nose.tools.assert_equals(self.ohbv.swe.all(), swe.all())
        nose.tools.assert_equals(self.ohbv.pond.all(), pond.all())

    def test_write_current_states(self):

        self.ohbv.write_current_states("08122006", "txt", lambda s, data: np.savetxt(s, data))
        swe = np.loadtxt("swe_08122006.txt")
        pond = np.loadtxt("pond_08122006.txt")
        melt = np.loadtxt("melt_08122006.txt")
        aet = np.loadtxt("aet_08122006.txt")

        nose.tools.assert_equals(self.ohbv.swe.all(), swe.all())
        nose.tools.assert_equals(self.ohbv.pond.all(), pond.all())
        nose.tools.assert_equals(self.ohbv.melt.all(), melt.all())
        nose.tools.assert_equals(self.ohbv.aet.all(), aet.all())










