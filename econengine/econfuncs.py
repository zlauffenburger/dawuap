from __future__ import division
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
        self.crop_id = np.asanyarray(kwargs.get('crop_id'), dtype=np.int)
        self.irr_eff = np.asanyarray(kwargs.get('irrigation_eff'), dtype=np.float)
        self.irr = np.asanyarray(kwargs.get('irrigation_mask'), dtype=bool)

        self.ref_et = np.asarray(kwargs['normalization_refs'].get('reference_et'))
        self.ref_prices = np.asarray(kwargs['normalization_refs'].get('reference_prices'))
        self.ref_yields = np.asarray(kwargs['normalization_refs'].get('reference_yields'))
        self.ref_land = np.asarray(kwargs['simulated_states'].get('used_land')).sum()

        self.sigmas = np.asarray(kwargs['parameters'].get('sigmas'))
        if len(self.sigmas) == 1:
            self.sigmas = np.repeat(self.sigmas, len(self.crop_list))

        self.deltas = np.asarray(kwargs['parameters'].get('deltas'))
        self.betas = np.asarray(kwargs['parameters'].get('betas'))
        self.mus = np.asarray(kwargs['parameters'].get('mus'))
        self.first_stage_lambda = np.asarray(kwargs['parameters'].get('first_stage_lambda'))
        self.lambdas_land = np.asarray(kwargs['parameters'].get('lambdas_land'))

        #self.prices = None
        #self.costs = kwargs.get('costs')  # land and water costs. Array with one row per crop. First column land second column water

        self._landsim = np.asarray(kwargs['simulated_states'].get('used_land')) / self.ref_land
        self._watersim = np.asarray(kwargs['simulated_states'].get('used_water')) / self.ref_et
        self.etasim = np.asarray(kwargs['simulated_states'].get('supply_elasticity_eta'))
        self._ysim = np.asarray(kwargs['simulated_states'].get('yields')) / self.ref_yields
        self.ysim_w = np.asarray(kwargs['simulated_states'].get('yield_elasticity_water'))

        # This to be filled with information provided during simulations
        self.crop_start_date = None
        self.crop_cover_date = None
        self.crop_end_date = None


        super(Farm, self).__init__(kwargs.get("id"), kwargs.get("source_id"), kwargs.get("name"))

    @property
    def landsim(self):
        return self._landsim * self.ref_land

    @landsim.setter
    def landsim(self, value):
        self._landsim = value / self.ref_land

    @property
    def watersim(self):
        """Returns un-normalized total simulated water use"""
        return self._watersim * self.ref_et

    @watersim.setter
    def watersim(self, value):
        self._watersim = value / self.ref_et

    @property
    def ysim(self):
        """Returns un-normalized simulated crop yields"""
        return self._ysim * self.ref_yields

    @ysim.setter
    def ysim(self, value):
        """Sets simulated yields"""
        self._ysim = value / self.ref_yields

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

    @staticmethod
    def _eta_sim(sigmas, delta, xbar, ybar_w, qbar, p):
        """
        Simulates exogenous supply elasticities (eta) as a function
         of observed land and water allocations and parameters

        :param delta: CES production function returns-to-scale parameter, 1D array
        :param xbar: Observed resource allocations, 2D array (ncrops x nresources)
        :param ybar_w: Yield elasticity with respect to water use, 1D array
        :param qbar: Observed total production for each crop
        :param p: Observed crop prices received by farmer for each crop, 1D array
        :return: vector of simulated supply elasticities of the same shape as delta
        """

        b = xbar[:, 0]**2 / (p * qbar)
        num = b / (delta * (1 - delta))
        dem = np.sum(num + (sigmas * b * ybar_w / (delta * (delta - ybar_w))))
        return delta / (1 - delta) * (1 - (num/dem))

    @staticmethod
    def _y_bar_w_sim(sigmas, beta, delta, xbar):
        """
        Simulates yield elasticity with respect to water (ybar_w) as a function of observed
        land and water allocation and parameters.

        :param sigmas: Elasticity of substitution parameter.
        :param beta: CES shares parameter, 2D array (ncrops x nresources).
        :param delta: CES production function returns-to-scale parameter, 1D array.
        :return: Vector of simulated yield elasticities with respect to water
         of the same shape as delta
        """
        r = rho(sigmas)
        num = beta[:, -1] * xbar[:, -1]**r
        den = np.diag(np.dot(beta, xbar.T**r))
        return delta * num/den

    @staticmethod
    def production_function(sigmas, beta, delta, mu, xbar, et0=[0]):
        """
        Constant elasticity of substitution production function

        :param sigmas: Elasticity of substitution parameter.
        :param beta: CES shares parameter, 2D array (ncrops x nresources).
        :param delta: CES production function returns-to-scale parameter, 1D array.
        :param mu: CES productivity parameter, 1D array.
        :param xbar: Resource allocation, 2D array (ncrops, nresources)
        :param et0: evapotranspiration, defaults to zero
        :return: vector of crop production with same shape as delta
        """
        r = rho(sigmas)
        beta = beta.clip(min=0, max=1)
        x = xbar.clip(min=0.0001)
        # adds evaporatranspiration to crop water if irrigation and et are dissagregated
        x[:, -1] = x[:, -1] + np.asarray(et0)
        return mu * np.diag(np.dot(beta, x.T**r))**(delta/r)

    @staticmethod
    def _first_stage_lambda_land_lhs(lambda_land, prices, costs, delta, qbar, y_bar_w, xbar):
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

    @staticmethod
    def _lambda_land_water_lhs(lambda_land, first_stage_lambda, deltas, prices, costs, qbar, xbar):

        fstStgLbda = np.array([first_stage_lambda, 0])
        # this following term only applies to land, hence 0 on the water column
        p_qbar_delta = np.asarray(prices * qbar * deltas)[:, np.newaxis]
        p_qbar_delta = np.append(p_qbar_delta, np.zeros_like(p_qbar_delta), axis=1)

        lhs_lambdas = (costs + lambda_land + fstStgLbda) * xbar - p_qbar_delta

        return lhs_lambdas

    @staticmethod
    def _convex_sum_constraint(betas):
        """sum the columns of CES production function share parameters"""

        return betas.sum(axis=1)

    @staticmethod
    def _observed_activity(prices, eta, ybar_w, ybar, xbar):
        """ produce the rhs of optimality equation by stacking
         the vectors of observations in a 1D vector"""

        qbar = ybar * xbar[:, 0]

        return np.hstack((eta, ybar_w, qbar, np.ones_like(prices), np.sum(2 * xbar[:, 0] * prices * qbar * ybar_w),
                               -prices*qbar*ybar_w, prices*qbar*ybar_w))

    def _set_reference_observations(self, **kwargs):

        eta = np.array(kwargs['eta'])
        ybar = np.array(kwargs['ybar']) / self.ref_yields
        landbar = np.array(kwargs['obs_land'])
        waterbar = np.array(kwargs['obs_water']) / self.ref_et
        xbar = np.array([landbar, waterbar]).T
        ybar_w = np.array(kwargs['ybar_w'])

        prices = np.array(kwargs['prices']) / self.ref_prices
        costs = np.array(kwargs['costs'])

        costs[:, 0] /= (self.ref_prices * self.ref_yields)
        costs[:, 1] *= self.ref_et/(self.ref_prices * self.ref_yields)

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

            sigmas = self.sigmas

            first_stage_lambda = pars[-1] # first stage lambda always the last parameter
            pars2 = pars[:-1].reshape(-1, prices.size).T
            deltas = pars2[:, 0]
            betas = pars2[:, 1:3]
            mus = pars2[:, 3]
            lambdas = pars2[:, 4:]
            rhs = self._observed_activity(prices, eta, ybar_w, ybar, xbar)

            lhs = np.hstack((
                self._eta_sim(sigmas, deltas, xbar, ybar_w, qbar, prices),
                self._y_bar_w_sim(sigmas, betas, deltas, xbar),
                self.production_function(sigmas, betas, deltas, mus, xbar),
                self._convex_sum_constraint(betas),
                self._first_stage_lambda_land_lhs(first_stage_lambda, prices, costs, deltas, qbar, ybar_w, xbar),
                self._lambda_land_water_lhs(lambdas, first_stage_lambda, deltas, prices, costs, qbar, xbar).T.flatten()))

            return lhs - rhs

        def calibrate(solve_pmp_program=True):
            x = np.hstack((self.deltas, self.betas.T.flatten(),
                           self.mus, self.lambdas_land.T.flatten(), self.first_stage_lambda))
            if solve_pmp_program:
                return sci.root(func, x, method='lm')
            else:
                return func(x)

        return calibrate

    def write_farm_dict(self, fname):
        """Dumps farm information to a dictionary and returns it and writes it to fname"""

        farm_dic = {
            "id": str(self.id),
            "source_id": str(self.source_id),
            "name": str(self.name),
            "crop_list": self.crop_list,
            "input_list": self.input_list,
            "irrigation_mask": self.irr.tolist(),
            "crop_id": self.crop_id.tolist(),
            "irrigation_eff": self.crop_id.tolist(),
            "parameters": {
                "sigmas": self.sigmas.tolist(),
                "deltas": self.deltas.tolist(),
                "betas": self.betas.tolist(),
                "mus": self.mus.tolist(),
                "first_stage_lambda": self.first_stage_lambda.tolist(),
                "lambdas_land": self.lambdas_land.tolist()
            },
            "constraints": {
                "land:": [-1],
                "water": [-1]
            },
            "simulated_states": {
                "used_land": self.landsim.tolist(),
                "used_water": self.watersim.tolist(),
                "supply_elasticity_eta": self.etasim.tolist(),
                "yields": self.ysim.tolist(),
                "yield_elasticity_water": self.ysim_w.tolist()
            },
            "normalization_refs": {
                "reference_et": self.ref_et.tolist(),
                "reference_prices": self.ref_prices.tolist(),
                "reference_yields": self.ref_yields.tolist(),
                "reference_land": self.ref_land.tolist()
            }
        }

        with open(fname, 'w') as file:
            file.write(json.dumps(farm_dic))

        return farm_dic

    def simulate(self, **kwargs):
        """Simulates resource allocation given given the current set of function parameters in the class.

        Parameters
        ==========

        :param kwargs:
            Dictionary with lists or arrays of prices, costs and constraints to production.

        :Example:

        ::

            observs = {
            'farm_id': 107,
            'evapotranspiration': [5., 5.]
            'prices': [5.82, 125],
            'costs': [111.56, 193.95],
            'land_constraint': 100,
            'water_constraint': 100,
            'crop_start_date': ["5/15/2014", "5/15/2014", "5/15/2014", "5/15/2014", "5/15/2014",
                                "5/15/2014", "5/15/2014", "5/15/2014"],
            'crop_cover_date': ["7/02/2014", "7/02/2014", "7/02/2014", "7/02/2014", "7/02/2014",
                                "7/02/2014", "7/02/2014", "7/02/2014"],
            'crop_end_date': ["8/25/2014", "8/25/2014", "8/25/2014", "8/25/2014", "8/25/2014",
                              "8/25/2014","8/25/2014", "8/25/2014"],
            }

            Farm_obj.simulate(**observs)

        """

        et0 = np.array(kwargs['evapotranspiration'])/self.ref_et
        prices = np.array(kwargs['prices'])/self.ref_prices
        if prices.ndim < 2:
            prices = prices[:, np.newaxis]

        costs = np.array(kwargs['costs'])

        costs[:, 0] /= (self.ref_prices * self.ref_yields)
        costs[:, 1] *= self.ref_et/(self.ref_prices * self.ref_yields)
        L = np.array(kwargs['land_constraint'])
        W = np.array(kwargs['water_constraint'])
        LW = np.hstack((L, W))

        self.crop_start_date = np.array(kwargs['crop_start_date'])
        self.crop_cover_date = np.array(kwargs['crop_cover_date'])
        self.crop_end_date = np.array(kwargs['crop_end_date'])

        def func(res):

            num_crops = len(self.crop_list)
            num_inputs = len(self.input_list)

            # parse inputs vector into convenient variables
            # first input resources for each crop (matrix x)
            # then simulated production (q)
            # finally lagrange multipliers
            x = res[:num_crops * num_inputs].reshape(num_inputs, num_crops).T
            xstar = x.copy()
            xstar[:, -1] += et0
            q = self.production_function(self.sigmas, self.betas, self.deltas, self.mus, x, et0)[:, np.newaxis]
            lbdas = res[(x.size ): (x.size ) + num_inputs]
            psi = res[((x.size ) + lbdas.size):][:, np.newaxis]

            # build the left hand side of system of equations
            r = rho(self.sigmas)[:, np.newaxis]
            num = prices * self.deltas[:, np.newaxis] * q * self.betas * (xstar ** r)
            xstar = xstar.clip(min=0.0001)
            den = np.diag(np.dot(self.betas, (xstar**r).T))
            drevdx = num / den[:, np.newaxis]

            pluset = np.zeros_like(x)
            pluset[:, 0] += et0

            xc = x[:, -1].copy()
            #xc[self.irr] = 0

            lhs = np.hstack((drevdx.T.flatten(),  (x + pluset).sum(axis=0), (psi*xc[:, np.newaxis]).flatten()))

            # build right hand side of system of equations

            dcdx = (costs + self.lambdas_land + lbdas + psi) * xstar
            #qbar = self.production_function(self.sigmas, self.betas, self.deltas, self.mus, x, et0)

            rhs = np.hstack((dcdx.T.flatten(), LW, np.zeros_like(xc)))

            #print np.sum(lhs - rhs)
            return lhs - rhs

        # prepare initial guesses

        q0 = self._ysim * self._landsim
        lam = np.zeros(len(self.input_list))
        lam_irr = np.zeros(len(self.crop_list))
        x0 = np.hstack((self._landsim, self._watersim, lam, lam_irr))
        output = sci.root(func, x0, method='lm')

        return output

    def calibrate(self, **kwargs):
        """Calibrates the economic model of agricultural production.

        Parameters
        ==========

        :param kwargs:
            Dictionary with lists or arrays of observed agricultural activity:

        :Example:

        ::

            observs = {
            'eta': [.35, 0.29],
            'ybar': [35, 2.2],
            'obs_land': [0.1220, 0.4.],
            'obs_water': [0.0250, 0.035],
            'ybar_w': [0.06, 0.21],
            'prices': [5.82, 125],
            'costs': [111.56, 193.95]}

            Farm_obj.calibrate(**observs)


        :return:
        """

        # set reference observations, solve nonlinear program and return results
        res = self._set_reference_observations(**kwargs)()

        if res.success:

            # retrieve optimized parameters, parse them and update member parameter variables
            pars = res['x']

            # Assign optimal values to Farm object
            self.first_stage_lambda = pars[-1]  # first stage lambda always the last parameter
            pars2 = pars[:-1].reshape(-1, len(self.crop_list)).T
            self.deltas = pars2[:, 0]
            self.betas = pars2[:, 1:3]
            self.mus = pars2[:, 3]
            self.lambdas_land = pars2[:, 4:]

        return res
