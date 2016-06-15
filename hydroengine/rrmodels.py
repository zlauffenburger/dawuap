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
        self.melt = np.zeros_like(self.swe)


    def snowpack(self, incid_precip, t_max, t_min):

        swe = np.zeros_like(self.swe)
        rain = np.zeros_like(self.swe)

        ind_allswe = np.less_equal(t_max, self.t_thres)
        ind_allrain = np.greater(t_min, self.t_thres)
        ind_mixed = np.logical_not(np.logical_or(ind_allswe, ind_allrain))

        swe[ind_allswe] = incid_precip[ind_allswe]
        rain[ind_allrain] = incid_precip[ind_allrain]
        swe[ind_mixed] = (incid_precip * (np.array((self.t_thres - t_min) / (t_max - t_min)).clip(min=0)))[ind_mixed]
        rain[ind_mixed] = (incid_precip - swe)[ind_mixed]

        self.swe += swe

        # Begin snowmelt routine

        tp = 0.5 * (t_max + t_min)

        ind_melt = np.greater(tp, self.t_thres)
        ind_swe = np.greater(self.swe, 0)

        # if average temp is above threshold and there is swe on the ground, produce melt
        ind = np.logical_and(ind_melt, ind_swe)
        self.melt[ind] = (self.ddf * (tp - self.t_thres))[ind]

        # if the amount of potential melt is larger than available swe, cap melt to available swe
        self.melt[np.greater(self.melt, self.swe)] = self.swe[np.greater(self.melt, self.swe)]

        self.swe -= self.melt

        self.pond += self.melt + rain

    def soil_processes(self):
        pass

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
