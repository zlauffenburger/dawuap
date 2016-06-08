from __future__ import division
import numpy as np


class hbv(object):
    """Implementation of the HVB rainfall-runoff model

    """

    def __init__(self):
        self.p_base = 10

        pass


    def snowpack(self, incid_precip, temp, params):

        # if(tmax < t_thres_):
        #     self.snow = incid_precip
        # else if(tmax)
        #
        # if(temp > self.Ts):
        #
        #     self.snow = self.ddf * (temp - self.Ts)
        # else:
        #     self.melt = 0;
        pass

    def excess_precip_to_runoff(self):
        u = lambda i: np.array(- ((self.p_base - 2 * i + 2) * np.abs(self.p_base - 2 * i + 2) + (2 * i - self.p_base) * np.abs(
            2 * i - self.p_base) - 4 * self.p_base)/(2 * self.p_base ^ 2)).clip(min=0)

        v = np.array([0,0,0,0, 2,3,0,2.3,1.4,2,0,0,0,0,0,0.2,1.2,2.2,1.5,0,0,0,0])
        Q = np.convolve(v,u(np.arange(1,10)))
        print np.trapz(v)
        print np.trapz(Q)
        print np.trapz(u(np.arange(1,10)))


        return Q

