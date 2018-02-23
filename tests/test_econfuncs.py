from econengine import econfuncs
import numpy as np
import nose


class TestEconfuncs(object):

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_rho(self):
        nose.tools.assert_almost_equal(econfuncs.rho(0.5), (0.5-1)/0.5)

class TestFarm(object):

    @classmethod
    def setup_class(cls):
        print "SETUP!"
        cls.id = 1
        cls.name = "Test Farm"
        cls.eta = [.35, 0.29, 0.29, 1.33, 0.38, 0.38, 0.35, 0.35]
        cls.sigma = 0.5

        cls.obs_land = [0.1220, 0.0250, 0.0078, 0.0328, 0.1636, 0.0051, 0.0189, 0.6247]
        cls.et0 = 5
        cls.obs_water = np.array([0, 0, 27.25, 0, 0, 20., 0, 0]) + cls.et0
        cls.xbar = np.array([cls.obs_land, cls.obs_water]).T

        cls.ybar = [35, 2.2, 5.4, 1.7, 30, 110, 36, 36]
        cls.p = [5.12, 134, 134, 10.5, 4.76, 4.76, 5.93, 3.88]

        cls.qbar = np.array(cls.obs_land) * np.array(cls.ybar)
        cls.ybar_w = [0.06, 0.21, 0.21, 0.26, 0.06, 0.06, 0.06, 0.06]
        cls.a = econfuncs.Farm(TestFarm.id,
                           TestFarm.name,
                           TestFarm.eta,
                           TestFarm.sigma, 6, 5, 6, 7)





    @classmethod
    def teardown_class(cls):
        print "TEAR DOWN!"

    def setup(self):
        pass

    def teardown(self):
        pass

    @nose.tools.raises(ValueError)
    def test_check_calibration_criteria1_raiseValueError(self):

        nose.tools.assert_raises(ValueError,
                                 TestFarm.a._check_calibration_criteria(
                                     np.array(TestFarm.obs_land),
                                     np.array(TestFarm.ybar_w)*10,
                                 np.array(TestFarm.qbar),
                                 np.array(TestFarm.p)))

    @nose.tools.raises(ValueError)
    def test_check_calibration_criteria2_raiseValueError(self):

        nose.tools.assert_raises(ValueError,
                                 TestFarm.a._check_calibration_criteria(
                                     np.array(TestFarm.obs_land),
                                     np.array(TestFarm.ybar_w)*10,
                                     np.array(TestFarm.qbar),
                                     np.array(TestFarm.p)))

    def test_eta_sim(self):
        # tests function eta for calibrated delta values such that simulated eta with given
        # delta produces observed eta
        delta = [0.276, 0.227, 0.225, 0.620,0.304, 0.276, 0.261, 0.387]
        np.testing.assert_allclose(TestFarm.a._eta_sim(
            np.array(delta),
            np.array(TestFarm.obs_land),
            np.array(TestFarm.ybar_w),
            np.array(TestFarm.qbar),
            np.array(TestFarm.p)),
            TestFarm.eta)

    def test_y_bar_w_sim(self):
        # tests function eta for calibrated delta values such that simulated eta with given
        # delta produces observed eta
        beta = [[0.418, 0.582],
                [0.016, 0.984],
                [0.002, 0.998],
                [0.217, 0.783],
                [0.448, 0.552],
                [0.126, 0.874],
                [0.402,0.598],
                [0.522, 0.478]]
        delta = [0.276, 0.227, 0.225, 0.620,0.304, 0.276, 0.261, 0.387]
        np.testing.assert_allclose(TestFarm.a._y_bar_w_sim(
            np.array(beta),
            np.array(delta),
            np.array(TestFarm.xbar)),
            TestFarm.ybar_w)





