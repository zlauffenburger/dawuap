from __future__ import division
from abc import ABCMeta, abstractmethod
import numpy as np
import rasterstats as rst
import rasterio as rio


class rrmodel(object):
    """
    This is a base class that defines the basic interface of rainfall runoff models
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def runoff(self):
        pass

class hbv(rrmodel):
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

        # self.hl1 = params['storage_parameter_1']
        # self.hl2 = params['storage_parameter_2']
        #
        # self.ck0 = params['surface_conductance']
        # self.ck1 = params['mid_layer_conductance']
        # self.ck2 = params['lower_layer_conductance']

        self.soils = {}

        # snow state and flux variables
        self.swe = swe_o
        self.pond = pond_o
        self.melt = np.zeros_like(self.swe)

        # soil state and flux variables
        self.sm = sm_o
        self.delta_sm = np.zeros_like(self.swe)
        self.ovlnd_flow = np.zeros_like(self.swe)
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
        self.ovlnd_flow[ind_sm_geq_fcap] = self.pond[ind_sm_geq_fcap] + (self.sm - self.fcap)[ind_sm_geq_fcap]
        self.sm[ind_sm_geq_fcap] = self.fcap

        # in all other cells calculate the portion of ponded water that goes into the soil storage
        ind_sm_lt_fcap = np.logical_not(ind_sm_geq_fcap)
        self.delta_sm[ind_sm_lt_fcap] = (self.pond * (1 - np.power((self.sm/self.fcap), self.beta)))[ind_sm_lt_fcap]
        self.sm[ind_sm_lt_fcap] += self.delta_sm[ind_sm_lt_fcap]
        self.ovlnd_flow[ind_sm_lt_fcap] = (self.pond - self.delta_sm)[[ind_sm_lt_fcap]]

        # Check if cell exceed storage capacity after adding delta_sm
        ind_sm_geq_fcap = np.greater_equal(self.sm, self.fcap)
        self.ovlnd_flow[ind_sm_geq_fcap] += (self.sm - self.fcap)[ind_sm_geq_fcap]
        self.sm[ind_sm_geq_fcap] = self.fcap

        # if there is sufficient soil moisture to satisfy aet, reduce sm
        ind_sm_gt_aet = np.greater(self.sm, self.aet)
        self.sm[ind_sm_gt_aet] -= self.aet[ind_sm_gt_aet]
        # otherwise take all water storage and limit aet
        ind_sm_leq_aet = np.logical_not(ind_sm_gt_aet)
        self.aet[ind_sm_leq_aet] = self.sm[ind_sm_leq_aet]
        self.sm[ind_sm_leq_aet] = 0.0


    def precipitation_excess(self, shp_wtshds, affine=None, stats=['mean']):

        # builds a geojson object with required statistics for each catchment
        self.stw1 = rst.zonal_stats(shp_wtshds, self.ovlnd_flow, nodata=-32768, affine=affine, geojson_out=True,
                                    prefix='runoff_', stats=stats)

        soil_layers = {}

        for i in range(len(self.stw1)):
            props = self.stw1[i]['properties']
            if props['runoff_mean'] > props['hbv_hl1']:
                soil_layers['Q0'] = (props['runoff_mean'] - props['hbv_hl1']) * props['hbv_ck0']
                props['runoff_mean'] -= soil_layers['Q0']
            else:
                soil_layers['Q0'] = 0.0

            if props['runoff_mean'] > 0.0:
                soil_layers['Q1'] = props['runoff_mean'] * props['hbv_ck1']
                props['runoff_mean'] -= soil_layers['Q1']
            else:
                soil_layers['Q1'] = 0.0

            if props['runoff_mean'] > props['hbv_perc']:
                props['runoff_mean'] -= props['hbv_perc']
                props['lower_reservoir'] += props['hbv_perc']
            else:
                props['lower_reservoir'] += props['runoff_mean']
                props['runoff_mean'] = 0.0

            if props['lower_reservoir'] > 0.0:
                soil_layers['Q2'] = props['lower_reservoir'] * props['hbv_ck2']
                props['lower_reservoir'] -= soil_layers['Q2']
            else:
                soil_layers['Q2'] = 0.0

            soil_layers['Qall'] = soil_layers['Q0'] + soil_layers['Q1'] + soil_layers['Q2']

        self.soils = dict(zip(self.stw1[i]['properties']['WSHD_ID']), soil_layers)


    def runoff(self):
        def u(i):
            return np.array(
                - ((self.p_base - 2 * i + 2) * np.abs(self.p_base - 2 * i + 2) + (2 * i - self.p_base) * np.abs(
                    2 * i - self.p_base) - 4 * self.p_base) / (2 * self.p_base ^ 2)).clip(min=0)

        v = np.array([0, 0, 0, 0, 2, 3, 0, 2.3, 1.4, 2, 0, 0, 0, 0, 0, 0.2, 1.2, 2.2, 1.5, 0, 0, 0, 0])
        Q = np.convolve(v, u(np.arange(1, 100)))
        print np.trapz(v)
        print np.trapz(Q)
        print 'area of triangle', np.trapz(u(np.arange(0, 12)))
        print self.p_base, u(np.arange(0, 12))

        return Q
