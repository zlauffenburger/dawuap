from unittest import TestCase

from .context import hydroengine
import numpy as np
import nose
import matplotlib.pyplot as plt


class Test_hbv(object):
    @classmethod
    def setup_class(klass):
        print "SETUP!"
        args = {
            'pp_temp_thres': 2,
            'p_base': 10,
            'ddf': 0.02
        }
        Test_hbv.swe_o = np.zeros((20, 30))
        Test_hbv.pond_o = np.zeros_like(Test_hbv.swe_o)

        Test_hbv.ohbv = hydroengine.hbv(Test_hbv.swe_o, Test_hbv.pond_o, **args)

    @classmethod
    def teardown_class(klass):
        print "TEAR DOWN!"

    def TestInitialization(self):
        np.testing.assert_array_equal(Test_hbv.swe_o, self.ohbv.swe)
        nose.tools.assert_equals(self.ohbv.p_base, 10)
        nose.tools.assert_equals(self.ohbv.t_thres, 2)
        nose.tools.assert_equals(self.ohbv.ddf, 0.02)

    def test_snowpack(self):
        precip = np.ones_like(Test_hbv.swe_o) * 0.03
        t_max = np.ones_like(Test_hbv.swe_o) * -2
        t_min = np.ones_like(Test_hbv.swe_o) * -3
        Test_hbv.ohbv.snowpack(precip, t_max, t_min)

        np.testing.assert_array_equal(self.ohbv.swe, np.ones_like(Test_hbv.swe_o)*0.03)

