from econengine import econfuncs
import numpy as np
import json
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
    def setup_class(self):
        print "SETUP!"
        self.id = 1
        self.name = "Test Farm"

        # reference values to adimensionalize inputs
        self.refet = 25 # Reference et for transpiration
        self.refprices = np.array([5.33, 112, 112, 121.85, 6.86, 6.86, 8.02, 6.20])
        self.refyields = np.array([52.56, 2.11, 2.11, 1.57, 29.01, 63.63, 29, 70.37])
        self.acrev = self.refprices * self.refyields # reference per-acre revenue
        #

        self.eta = np.array([.35, 0.29, 0.29, 1.33, 0.38, 0.38, 0.35, 0.35])
        self.sigma = 0.5

        self.deltas = np.array([0.276, 0.227, 0.225, 0.574, 0.299, 0.276, 0.261, 0.507])
        self.betas = np.array([[0.947, 0.053],
                               [0.288, 0.712],
                               [0.053, 0.947],
                               [0.858, 0.142],
                               [0.952, 0.048],
                               [0.782, 0.218],
                               [0.944, 0.056],
                               [0.974, 0.026]])

        self.mus = np.array([0.153, 0.082, 0.056, 0.327, 0.306, 0.038, 0.070, 0.427])
        self.first_stage_lambda = np.array([-0.165])
        self.lambdas_land = np.array([[-0.076, -0.636, -1.443, -0.503, -0.270, -0.411, -0.202, 0.024],
                                 [0, 0, -6.482, 0, 0, -2.688, 0, 0]]).T

        self.obs_land = np.array([0.1220, 0.0250, 0.0078, 0.0328, 0.1636, 0.0051, 0.0189, 0.6247])

        self.et0 = 5
        self.obs_water = (np.array([0, 0, 27.25, 0, 0, 20., 0, 0]) + self.et0) * self.obs_land
        self.obs_water /= self.refet

        self.xbar = np.array([self.obs_land, self.obs_water]).T

        self.prices = np.array([5.82, 125, 125, 111.72, 4.80, 4.80, 6.88, 4.59])
        self.prices /=self.refprices # normalize prices

        self.costs = np.array([[111.56, 193.95, 390.02, 187.38, 120.80, 365.33, 135.13, 135.13],
                               [0, 0, 65.67, 0, 0, 48.20, 0, 0]]).T

        self.costs[:, 0] /= self.acrev
        self.costs[:, 1] *= self.refet/self.acrev

        self.ybar = np.array([35, 2.2, 5.4, 1.7, 30, 110, 36, 36])
        self.ybar /= self.refyields

        self.qbar = self.obs_land * np.array(self.ybar)
        self.ybar_w = np.array([0.06, 0.21, 0.21, 0.26, 0.06, 0.06, 0.06, 0.06])

        params = {
            'eta': self.eta,
            'sigma': self.sigma,
            'deltas': self.deltas,
            'betas': self.betas,
            'mus': self.mus,
            'first_stage_lambda': self.first_stage_lambda,
            'lambdas': self.lambdas_land,
            'costs': self.costs,
            'obs_allocation': self.xbar
        }

        with open('test_data/Farms.json') as json_farms:
            farms = json.load(json_farms)

        farm1 = farms['farms'][0]

        self.a = econfuncs.Farm(**farm1)

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
                                     np.array(TestFarm.sigma),
                                     np.array(TestFarm.eta),
                                     np.array(TestFarm.obs_land),
                                     np.array(TestFarm.ybar_w)*10,
                                 np.array(TestFarm.qbar),
                                 np.array(TestFarm.prices)))

    @nose.tools.raises(ValueError)
    def test_check_calibration_criteria2_raiseValueError(self):

        nose.tools.assert_raises(ValueError,
                                 TestFarm.a._check_calibration_criteria(
                                     np.array(TestFarm.sigma),
                                     np.array(TestFarm.eta),
                                     np.array(TestFarm.obs_land),
                                     np.array(TestFarm.ybar_w)*10,
                                     np.array(TestFarm.qbar),
                                     np.array(TestFarm.prices)))

    def test_eta_sim(self):
        # tests function eta for calibrated delta values such that simulated eta with given
        # delta produces observed eta
        delta = TestFarm.deltas

        np.testing.assert_allclose(TestFarm.a._eta_sim(
            np.array(TestFarm.sigma),
            np.array(delta),
            np.array(TestFarm.xbar),
            np.array(TestFarm.ybar_w),
            np.array(TestFarm.qbar),
            np.array(TestFarm.prices)),
            TestFarm.eta, rtol=1e-2)

    def test_y_bar_w_sim(self):
        # tests function eta for calibrated delta values such that simulated eta with given
        # delta produces observed eta
        beta = TestFarm.betas
        delta = TestFarm.deltas
        np.testing.assert_allclose(TestFarm.a._y_bar_w_sim(
            np.array(TestFarm.sigma),
            np.array(beta),
            np.array(delta),
            np.array(TestFarm.xbar)),
            TestFarm.ybar_w, rtol=1e-2)

    def test_prod_func(self):
        # tests production function for calibrated mu values such that simulated production with given
        # mu produces observed agricultural production
        beta = self.betas
        delta = self.deltas
        mu = self.mus

        np.testing.assert_allclose(TestFarm.a.production_function(
            np.array(TestFarm.sigma),
            np.array(beta),
            np.array(delta),
            np.array(mu),
            np.array(TestFarm.xbar)),
            TestFarm.ybar * TestFarm.xbar[:, 0], rtol=1e-2)

    def test_first_stage_lambda_land_lhs(self):
        # test optimality condition for land shadow price (lambda).
        # The function should return zero for lambda value that minimizes SSE
        # of constrained net revenue function

        import scipy.optimize as opt

        p = self.prices
        qbar = self.qbar
        delta = self.deltas
        ybarw = self.ybar_w
        xbar = self.xbar
        c = self.costs

        fun = lambda lam: np.sum((p * qbar * (delta - ybarw) - (c[:, 0] + lam) * xbar[:, 0])**2 +
                                 (p * qbar * ybarw - c[:, 1] * xbar[:, 1])**2)

        def jac(lam):
            return np.sum(-2 * xbar[:, 0] * (p * qbar * (delta - ybarw) - (c[:, 0] + lam) * xbar[:, 0]))

        lambda_opt = opt.minimize(fun, np.array([0]), jac=jac, method='Newton-CG')

        print lambda_opt

        np.testing.assert_allclose(
            TestFarm.a._first_stage_lambda_land_lhs(lambda_opt['x'], p, c, delta, qbar, ybarw, xbar),
            np.sum(2 * xbar[:, 0] * p * qbar * ybarw), rtol=1e-2)

    def test_set_reference_observations(self):

        observs = {
            'eta': self.eta,
            'ybar': self.ybar,
            'xbar': self.xbar,
            'ybar_w': self.ybar_w,
            'prices': self.prices,
            'costs': self.costs
        }

        pmp = TestFarm.a._set_reference_observations(**observs)
        res = pmp()
        print res
        nose.tools.assert_equals(res.success, True)










