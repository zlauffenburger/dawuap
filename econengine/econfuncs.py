from water_user import WaterUser
import numpy as np
import scipy.optimize as sci


def rho(sigma):
    """returns elasticity of substitution rho

    :param sigma: substitution elasticity by crop and technology
    :return: rho parameter
    """
    return (sigma - 1) / sigma


class Farm(WaterUser):
    """Class representing economic behavior of farms"""

    def __init__(self, identifier, name, **kwargs):

        self.sigma = kwargs.get('sigma')
        self.rho = rho(self.sigma)

        self.deltas = kwargs.get('deltas')
        self.betas = kwargs.get('betas')
        self.mus = kwargs.get('mus')
        self.first_stage_lambda = kwargs.get('first_stage_lambda')
        self.lambdas = kwargs.get('lambdas')

        self.prices = None
        self.costs = kwargs.get('costs')  # land and water costs. Array with one row per crop. First column land second column water

        self.xbar = kwargs.get('obs_allocation')
        self.eta = kwargs.get('eta')
        self.ybar = None
        self.ybar_w = None

        self.qbar = None

        super(Farm, self).__init__(identifier, name)

    def _check_calibration_criteria(self, xbar, ybar_w, qbar, p):

        # Check calibration criteria 1
        if ((self.eta - ybar_w/(1 - ybar_w)) < 0).any():
            raise ValueError('calibration criteria 1'
                             'for farm %i with name %s failed' % (self.id, self.name))

        # Check calibration criteria 2
        b = xbar**2/(p * qbar)
        psi = self.sigma*ybar_w / (self.eta * (1 - ybar_w))
        ind = np.arange(len(b))
        cc2 = b * self.eta * (1 - psi) - [np.sum(b[ind != i] *
                                                 self.eta[ind != i] * (1 + (1 / self.eta[ind != i]))**2 *
                                                 (1 + psi[ind != i] - ybar_w[ind != i])) for i in ind]
        if (cc2 > 0).any():
            raise ValueError('calibration criteria 2'
                             'for farm %i with name %s failed' % (self.id, self.name))

    def _eta_sim(self, delta, xbar, ybar_w, qbar, p):
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
        dem = np.sum(num + (self.sigma*b*ybar_w/(delta * (delta - ybar_w))))
        return delta / (1 - delta) #* (1 - (num/dem))

    def _y_bar_w_sim(self, beta, delta, xbar):
        """
        Simualted yield elasticity (ybar_w) with respect to water

        :return:
        """
        r = rho(self.sigma)
        num = beta[:, -1] * xbar[:, -1]**r
        den = np.diag(np.dot(beta, xbar.T**r))
        return delta * num/den

    def production_function(self, beta, delta, mu, xbar):
        """
        Constant elasticity of substitution production function

        :return:
        """
        r = rho(self.sigma)
        beta = beta.clip(min=0, max=1)
        return mu * np.diag(np.dot(beta, xbar.T**r))**(delta/r)

    def _first_stage_lambda_land_lhs(self, lambda_land, prices, delta, qbar, y_bar_w, xbar):
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

        condition = -2. * (self.costs[:, 0] + lambda_land) * xbar[:, 0]**2 + 2 * xbar[:, 0] * prices * qbar * delta

        return np.sum(condition)

    def _lambda_land_water_lhs(self, lambda_land, first_stage_lambda, xbar):

        l = np.array([first_stage_lambda, 0])
        rhs_lambdas = (self.costs + lambda_land + l) * xbar

        return rhs_lambdas

    def _convex_sum_constraint(self, betas):

        return betas.sum(axis=1)

    def _observed_activity(self, prices, ybar_w, ybar, xbar):

        qbar = ybar * xbar[:, 0]

        return np.hstack((self.eta, ybar_w, qbar, np.ones_like(prices), np.sum(2 * xbar[:, 0] * prices * qbar * ybar_w),
                               -prices*qbar*ybar_w, prices*qbar*ybar_w))

    def set_reference_observations(self, **kwargs):

        try:
            self._check_calibration_criteria(kwargs['xbar'][:,0],
                                             kwargs['ybar_w'],
                                             kwargs['xbar'][:,0] * kwargs['ybar'],
                                             kwargs['prices']
                                             )
        except ValueError as e:
            print "Flag raised for inconsistent observations with message: ", e
            print "NEW OBSERVATIONS NOT INCORPORATED INTO FARM... "
            return None

        self.eta = kwargs['eta']
        self.ybar = kwargs['ybar']
        self.xbar = kwargs['xbar']
        self.ybar_w = kwargs['ybar_w']

        self.prices = kwargs['prices']
        self.costs = kwargs['costs']

        self.qbar = self.ybar * self.xbar[:, 0]

        def func(pars):

            first_stage_lambda = pars[-1] # first stage lambda always the last parameter
            pars2 = pars[:-1].reshape(-1, self.prices.size).T
            deltas = pars2[:, 0]
            betas = pars2[:, 1:3]
            mus = pars2[:, 3]
            lambdas = pars2[:, 4:]
            rhs = self._observed_activity(self.prices, self.ybar_w, self.ybar, self.xbar)

            lhs = np.hstack((
                self._eta_sim(deltas, self.xbar, self.ybar_w, self.qbar, self.prices),
                self._y_bar_w_sim(betas, deltas, self.xbar),
                self.production_function(betas, deltas, mus, self.xbar),
                self._convex_sum_constraint(betas),
                self._first_stage_lambda_land_lhs(first_stage_lambda, self.prices, deltas, self.qbar, self.ybar_w, self.xbar),
                self._lambda_land_water_lhs(lambdas, first_stage_lambda, self.xbar).T.flatten()))

            return lhs - rhs

        def calibrate():
            x = np.hstack((self.deltas, self.betas.T.flatten(), self.mus, self.lambdas.T.flatten(), self.first_stage_lambda))
            return sci.root(func, x, method='lm')

        return calibrate

    def simulate(self):
        pass

    def calibrate(self, xbar, ybar_w):
        self.xbar = np.array(xbar)
        self.ybar_w = np.array(ybar_w)

