from unittest import TestCase

from .context import hydroengine
import numpy as np
import nose
import matplotlib.pyplot as plt


class Test_hbv(TestCase):

    @classmethod
    def setup_class(klass):
        print "setting up class " + klass.__name__
        Test_hbv.args = {
            'pp_temp_thres': 2,
            'p_base': 10,
            'ddf': 0.02
        }
        Test_hbv.swe_o = np.zeros((20, 30))
        Test_hbv.pond_o = np.zeros_like(Test_hbv.swe_o)


    @classmethod
    def teardown_class(klass):
        print "TEAR DOWN class " + klass.__name__

    def setUp(self):
        print "setting up " + self.__class__.__name__
        self.ohbv = hydroengine.hbv(Test_hbv.swe_o, Test_hbv.pond_o, **Test_hbv.args)

    def tearDown(self):
        Test_hbv.swe_o = np.zeros((20, 30))
        Test_hbv.pond_o = np.zeros_like(Test_hbv.swe_o)
        pass

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

