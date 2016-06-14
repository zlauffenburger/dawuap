from __future__ import division
import numpy as np


class hbv(object):
    """Implementation of the HVB rainfall-runoff model

    """

    def __init__(self, swe_o, pond_o, **params):

        self.t_thres = params['pp_temp_thres']
        self.p_base = params['p_base']
        self.ddf = params['ddf']
        self.swe = swe_o
        self.pond = pond_o

    def snowpack(self, incid_precip, t_max, t_min):

        rain = 0
        if (t_max <= self.t_thres):
            swe = incid_precip
            rain = 0
        elif (t_min > self.t_thres):
            rain = incid_precip
            swe = 0
        else:
            swe = incid_precip * np.array(self.t_thres - t_min / (t_max - t_min)).clip(min=0)

        self.swe += swe

        tp = 0.5 * (t_max + t_min)

        self.melt[np.greater(tp, self.t_thres)]

        if (tp > self.t_thres):
            if (self.swe > 0):
                self.melt = self.ddf * (tp - self.t_thres)
                self.swe -= self.melt
                if (self.melt > self.swe):
                    self.melt = self.swe
                    self.swe = 0

        else:
            self.melt = 0;

        self.pond += self.melt + rain

    def excess_precip_to_runoff(self):
        def u(i):
            return np.array(
                - ((self.p_base - 2 * i + 2) * np.abs(self.p_base - 2 * i + 2) + (2 * i - self.p_base) * np.abs(
                    2 * i - self.p_base) - 4 * self.p_base) / (2 * self.p_base ^ 2)).clip(min=0)

        v = np.array([0, 0, 0, 0, 2, 3, 0, 2.3, 1.4, 2, 0, 0, 0, 0, 0, 0.2, 1.2, 2.2, 1.5, 0, 0, 0, 0])
        Q = np.convolve(v, u(np.arange(1, 10)))
        print np.trapz(v)
        print np.trapz(Q)
        print np.trapz(u(np.arange(1, 10)))

        return Q
