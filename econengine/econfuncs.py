from water_user import WaterUser
import numpy as np


def rho(sigma):
    """returns elasticity of substitution rho

    :param sigma: substitution elasticity by crop and technology
    :return: rho parameter
    """
    return (sigma - 1) / sigma


class Farm(WaterUser):
    """Class representing economic behavior of farms"""

    def __init__(self, identifier, name, eta, sigma, share, prices, costs,  L):
        self.eta = np.array(eta)
        self.rho = rho(sigma)
        self.sigma = sigma

        self.land_ref = np.array(L)
        super(Farm, self).__init__(identifier, name)

    def _check_calibration_criteria(self, xbar, ybar_w, qbar, p):

        # Check calibration criteria 1
        if ((self.eta - ybar_w/(1 - ybar_w)) < 0).any():
            raise ValueError('calibration criteria 1'
                             'for farm %i with name %s failed' % (self.id, self.name))

        # Check calibration criteria 2
        b = np.sqrt(xbar)/(p * qbar)
        psi = self.sigma*ybar_w / (self.eta * (1 - ybar_w))
        ind = np.arange(len(b))
        cc2 = b * self.eta * (1 - psi) - [np.sum(b[ind != i] *
                                                 self.eta[ind != i] * np.sqrt(1 + (1 / self.eta[ind != i])) *
                                                 (1 + psi[ind != i] - ybar_w[ind != i])) for i in ind]
        if (cc2 > 0).any():
            print cc2
            raise ValueError('calibration criteria 2'
                             'for farm %i with name %s failed' % (self.id, self.name))

    def simulate(self):
        pass

    def calibrate(self, xbar, ybar_w):
        self.xbar = np.array(xbar)
        self.ybar_w = np.array(ybar_w)

