from __future__ import division
import numpy as np
import rasterstats as rst
import rasterio as rio


class hbv(object):
    """Implementation of the classic HVB rainfall-runoff model

    """

    def __init__(self, swe_o, pond_o, sm_o, stw1_o, stw2_o, **params):

        # Geomtry information
        #self.pixel_area = params['image_res']*params['image_res']
        #self.catchment_area = params['catch_area']
        #self.t_step = params['time_step']

        # snow paramters
        self.t_thres = params['pp_temp_thres']
        self.p_base = params['p_base']
        self.ddf = params['ddf']

        # soil paramters
        self.fcap = params['soil_max_wat']
        self.beta = params['soil_beta']
        self.lp = params['aet_lp_param']


        # snow state and flux variables
        self.swe = swe_o
        self.pond = pond_o
        self.melt = np.zeros_like(self.swe)

        # soil state and flux variables
        self.sm = sm_o
        self.delta_sm = np.zeros_like(self.swe)
        self.runoff = np.zeros_like(self.swe)
        self.stw1 = stw1_o
        self.stw2 = stw2_o

        self.aet = np.zeros_like(self.swe);


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

    def soil_processes(self, pot_et):

        # calculate AET
        self.aet = (pot_et * self.sm/(self.fcap*self.lp)).clip(min=0.0, max=1.0)

        # For cell where soil moisture is already at capacity, runoff is all the ponded water plus excess water in soil
        ind_sm_geq_fcap = np.greater_equal(self.sm, self.fcap)
        self.runoff[ind_sm_geq_fcap] = self.pond[ind_sm_geq_fcap] + (self.sm - self.fcap)[ind_sm_geq_fcap]
        self.sm[ind_sm_geq_fcap] = self.fcap

        # in all other cells calculate the portion of ponded water that goes into the soil storage
        ind_sm_lt_fcap = np.logical_not(ind_sm_geq_fcap)
        self.delta_sm[ind_sm_lt_fcap] = (self.pond * (1 - np.power((self.sm/self.fcap), self.beta)))[ind_sm_lt_fcap]
        self.sm[ind_sm_lt_fcap] += self.delta_sm[ind_sm_lt_fcap]
        self.runoff[ind_sm_lt_fcap] = (self.pond - self.delta_sm)[[ind_sm_lt_fcap]]

        # Check if cell exceed storage capacity after adding delta_sm
        ind_sm_geq_fcap = np.greater_equal(self.sm, self.fcap)
        self.runoff[ind_sm_geq_fcap] += (self.sm - self.fcap)[ind_sm_geq_fcap]
        self.sm[ind_sm_geq_fcap] = self.fcap

        # if there is sufficient soil moisture to satisfy aet, reduce sm
        ind_sm_gt_aet = np.greater(self.sm, self.aet)
        self.sm[ind_sm_gt_aet] -= self.aet[ind_sm_gt_aet]
        # otherwise take all water storage and limit aet
        ind_sm_leq_aet = np.logical_not(ind_sm_gt_aet)
        self.aet[ind_sm_leq_aet] = self.sm[ind_sm_leq_aet]
        self.sm[ind_sm_leq_aet] = 0.0

    def discharge(self, shp_poly, affine=None):

        # this is for testing purposes
        affine = rio.Affine
        src = rio.open('./tests/test_data/DEM_64m_1992.tif')
        img = src.read()



        # Add runoff to intermediate tank
        self.stw1 = np.sum(self.runoff * self.pixel_area)/self.catchment_area






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
