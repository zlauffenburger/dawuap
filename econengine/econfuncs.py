from water_user import WaterUser
import numpy as np
import json
import scipy.optimize as sci


def rho(sigma):
    """returns elasticity of substitution rho

    :param sigma: substitution elasticity by crop and technology
    :return: rho parameter
    """
    return (sigma - 1) / sigma


class Farm(WaterUser):
    """Class representing economic behavior of farms"""

    def __init__(self, fname=None, **kwargs):

        if isinstance(fname, str):
            with open(fname) as json_data:
                kwargs = json.load(json_data)

        self.crop_list = kwargs.get('crop_list')
        self.input_list = kwargs.get('input_list')

        self.sigmas = np.asarray(kwargs['parameters'].get('sigmas'))

        self.deltas = np.asarray(kwargs['parameters'].get('deltas'))
        self.betas = np.asarray(kwargs['parameters'].get('betas'))
        self.mus = np.asarray(kwargs['parameters'].get('mus'))
        self.first_stage_lambda = np.asarray(kwargs['parameters'].get('first_stage_lambda'))
        self.lambdas_land = np.asarray(kwargs['parameters'].get('lambdas_land'))

        #self.prices = None
        #self.costs = kwargs.get('costs')  # land and water costs. Array with one row per crop. First column land second column water

        self.landsim = kwargs['simulated_states'].get('used_land')
        self.watersim = kwargs['simulated_states'].get('used_water')
        self.etasim = kwargs['simulated_states'].get('supply_elasticity_eta')
        self.ysim =kwargs['simulated_states'].get('yields')
        self.ysim_w = kwargs['simulated_states'].get('yield_elasticity_water')

        super(Farm, self).__init__(kwargs.get("id"), kwargs.get("name"))

    def _check_calibration_criteria(self, sigmas, eta, xbar, ybar_w, qbar, p):

        # Check calibration criteria 1
        if ((eta - ybar_w/(1 - ybar_w)) < 0).any():
            raise ValueError('calibration criteria 1'
                             'for farm %i with name %s failed' % (self.id, self.name))

        # Check calibration criteria 2
        b = xbar**2/(p * qbar)
        psi = sigmas * ybar_w / (eta * (1 - ybar_w))
        ind = np.arange(len(b))
        cc2 = b * eta * (1 - psi) - [np.sum(b[ind != i] *
                                                 eta[ind != i] * (1 + (1 / eta[ind != i]))**2 *
                                                 (1 + psi[ind != i] - ybar_w[ind != i])) for i in ind]
        if (cc2 > 0).any():
            raise ValueError('calibration criteria 2'
                             'for farm %i with name %s failed' % (self.id, self.name))

    def _eta_sim(self, sigmas, delta, xbar, ybar_w, qbar, p):
        """
        Simulated exogenous supply elasticities (eta) as a function of observed and prescribed parameters

        :param delta: CES production function parameter, 1D array
        :param xbar: Observed resource allocations, 2D array
        :param ybar_w: observed
        :param qbar:
        :param p:
        :return:
        """

        b = xbar[:, 0]**2 / (p * qbar)
        num = b / (delta * (1 - delta))
        dem = np.sum(num + (sigmas * b * ybar_w / (delta * (delta - ybar_w))))
        return delta / (1 - delta) #* (1 - (num/dem))

    def _y_bar_w_sim(self, sigmas, beta, delta, xbar):
        """
        Simualted yield elasticity (ybar_w) with respect to water

        :return:
        """
        r = rho(sigmas)
        num = beta[:, -1] * xbar[:, -1]**r
        den = np.diag(np.dot(beta, xbar.T**r))
        return delta * num/den

    def production_function(self, sigmas, beta, delta, mu, xbar):
        """
        Constant elasticity of substitution production function

        :return:
        """
        r = rho(sigmas)
        beta = beta.clip(min=0, max=1)
        return mu * np.diag(np.dot(beta, xbar.T**r))**(delta/r)

    def _first_stage_lambda_land_lhs(self, lambda_land, prices, costs, delta, qbar, y_bar_w, xbar):
        """
        First order optimality condition for the calibration of the land shadow value.
        Shadow value is calibrated to observed states when the function returns 0

        :param beta:
        :param delta:
        :param mu:
        :param xbar:
        :return:
        """
        #qbar  = #self.production_function(beta, delta, mu, xbar)
        #yw = #self._y_bar_w_sim(beta, delta, xbar)
        lambda_land = np.asarray(lambda_land)
        prices = np.asarray(prices)
        delta = np.asarray(delta)
        qbar = np.asarray(qbar)
        y_bar_w = np.asarray(y_bar_w)
        xbar = np.asarray(xbar)

        condition = -2. * (costs[:, 0] + lambda_land) * xbar[:, 0]**2 + 2 * xbar[:, 0] * prices * qbar * delta

        return np.sum(condition)

    def _lambda_land_water_lhs(self, lambda_land, first_stage_lambda, costs, xbar):

        l = np.array([first_stage_lambda, 0])
        rhs_lambdas = (costs + lambda_land + l) * xbar

        return rhs_lambdas

    def _convex_sum_constraint(self, betas):

        return betas.sum(axis=1)

    def _observed_activity(self, prices, eta, ybar_w, ybar, xbar):

        qbar = ybar * xbar[:, 0]

        return np.hstack((eta, ybar_w, qbar, np.ones_like(prices), np.sum(2 * xbar[:, 0] * prices * qbar * ybar_w),
                               -prices*qbar*ybar_w, prices*qbar*ybar_w))

    def set_reference_observations(self, **kwargs):

        eta = kwargs['eta']
        ybar = kwargs['ybar']
        xbar = kwargs['xbar']
        ybar_w = kwargs['ybar_w']

        prices = kwargs['prices']
        costs = kwargs['costs']

        qbar = ybar * xbar[:, 0]

        try:
            self._check_calibration_criteria(self.sigmas,
                                             eta,
                                             xbar[:, 0],
                                             ybar_w,
                                             qbar,
                                             prices
                                             )
        except ValueError as e:
            print "Flag raised for inconsistent observations with message: ", e
            print "NEW OBSERVATIONS NOT INCORPORATED INTO FARM... "
            return None

        def func(pars):

            first_stage_lambda = pars[-1] # first stage lambda always the last parameter
            pars2 = pars[:-1].reshape(-1, prices.size).T
            deltas = pars2[:, 0]
            betas = pars2[:, 1:3]
            mus = pars2[:, 3]
            lambdas = pars2[:, 4:]
            rhs = self._observed_activity(prices, eta, ybar_w, ybar, xbar)

            lhs = np.hstack((
                self._eta_sim(deltas, xbar, ybar_w, qbar, prices),
                self._y_bar_w_sim(betas, deltas, xbar),
                self.production_function(betas, deltas, mus, xbar),
                self._convex_sum_constraint(betas),
                self._first_stage_lambda_land_lhs(first_stage_lambda, prices, costs, deltas, qbar, ybar_w, xbar),
                self._lambda_land_water_lhs(lambdas, first_stage_lambda, costs, xbar).T.flatten()))

            return lhs - rhs

        def calibrate():
            x = np.hstack((self.deltas, self.betas.T.flatten(), self.mus, self.lambdas_land.T.flatten(), self.first_stage_lambda))
            return sci.root(func, x, method='lm')

        return calibrate

    def simulate(self):
        pass

    def calibrate(self, xbar, ybar_w):
        self.xbar = np.array(xbar)
        self.ybar_w = np.array(ybar_w)

