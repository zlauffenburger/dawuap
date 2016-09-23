from __future__ import division
from abc import ABCMeta, abstractmethod
import cPickle as pickle
import numpy as np
import rasterstats as rst
import rasterio as rio


class Soil(object):
    def __init__(self, ur=0, lr=0, q0=0, q1=0, q2=0, qall=0, base=0):

        self.upper_reservoir = ur
        self.lower_reservoir = lr

        self.Q0 = q0
        self.Q1 = q1
        self.Q2 = q2
        self.Qall = qall

        self.base = base


    @property
    def runoff(self, value):
        self._runoff = value

    @property.getter
    def runoff(self):
        return self._runoff[0]

    def __hash__(self):
        return hash(self.upper_reservoir,
                    self.lower_reservoir,
                    self.Q0,
                    self.Q1,
                    self.Q2,
                    self.Qall)

    def __eq__(self, other):
        return (self.upper_reservoir,
                    self.lower_reservoir,
                    self.Q0,
                    self.Q1,
                    self.Q2,
                    self.Qall) == (other.upper_reservoir,
                                   other.lower_reservoir,
                                   other.Q0,
                                   other.Q1,
                                   other.Q2,
                                   other.Qall)

    def __ne__(self, other):
        return not (self == other)


class RRmodel(object):
    """
    This is a base class that defines the basic interface of rainfall runoff models
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def runoff(self):
        pass


class HBV(RRmodel):
    """Implementation of the classic HVB rainfall-runoff model

    """

    def __init__(self, swe_o, pond_o, sm_o, soils_o={}, **params):

        # Geomtry information
        #self.pixel_area = params['image_res']*params['image_res']
        #self.catchment_area = params['catch_area']
        #self.t_step = params['time_step']

        # snow paramters
        self.t_thres = params['pp_temp_thres']
        self.ddf = params['ddf']

        # distributed soil paramters
        self.fcap = params['soil_max_wat']
        self.beta = params['soil_beta']
        self.lp = params['aet_lp_param']

        # self.hl1 = params['storage_parameter_1']
        # self.hl2 = params['storage_parameter_2']
        #
        # self.ck0 = params['surface_conductance']
        # self.ck1 = params['mid_layer_conductance']
        # self.ck2 = params['lower_layer_conductance']

        self.soils = soils_o # Soil reinitialization file. Dictionary with geo_json objects

        # snow state and flux variables
        self.swe = swe_o
        self.pond = pond_o
        self.melt = np.zeros_like(self.swe)

        # soil state and flux variables
        self.sm = sm_o
        self.delta_sm = np.zeros_like(self.swe)
        self.ovlnd_flow = np.zeros_like(self.swe)

        self.aet = np.zeros_like(self.swe)

    def snow_pack(self, incid_precip, t_max, t_min):

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
        stw1 = rst.zonal_stats(shp_wtshds, self.ovlnd_flow, nodata=-32768, affine=affine, geojson_out=True,
                                    prefix='runoff_', stats=stats)

        soil_layers = {}

        for i in range(len(stw1)):
            if not self.soils: # if there is no dictionary from a previous time step, create a new soil object
                soil_layers = Soil()
            else:
                soil_layers = self.soils[i][1]

            props = stw1[i]['properties']
            soil_layers.upper_reservoir += props['runoff_mean']
            if soil_layers.upper_reservoir > props['hbv_hl1']:
                soil_layers.Q0 = (soil_layers.upper_reservoir - props['hbv_hl1']) * props['hbv_ck0']
                soil_layers.upper_reservoir -= soil_layers.Q0
            else:
                soil_layers.Q0 = 0.0

            if soil_layers.upper_reservoir > 0.0:
                soil_layers.Q1 = soil_layers.upper_reservoir * props['hbv_ck1']
                soil_layers.upper_reservoir -= soil_layers.Q1
            else:
                soil_layers.Q1 = 0.0

            if soil_layers.upper_reservoir > props['hbv_perc']:
                soil_layers.upper_reservoir -= props['hbv_perc']
                soil_layers.lower_reservoir += props['hbv_perc']
            else:
                soil_layers.lower_reservoir += soil_layers.upper_reservoir
                soil_layers.upper_reservoir = 0.0

            if soil_layers.lower_reservoir > 0.0:
                soil_layers.Q2 = soil_layers.lower_reservoir * props['hbv_ck2']
                soil_layers.lower_reservoir -= soil_layers.Q2
            else:
                soil_layers.Q2 = 0.0

            soil_layers.Qall = soil_layers.Q0 + soil_layers.Q1 + soil_layers.Q2

        self.soils = dict(zip(stw1[i]['properties']['WSHD_ID']), soil_layers)
        pickle.dump(self.soils, open("soils.pickled", "wb"))

    def calculate_runoff(self, i):
        """

        :param i:
        :return:
        """
        def u(j, p_base):
            """

            :type p_base: scalar with base time of unit hydrograph
            """
            return np.array(
                - ((p_base - 2 * j + 2) * np.abs(p_base - 2 * j + 2) + (2 * j - p_base) * np.abs(
                    2 * j - p_base) - 4 * p_base) / (2 * p_base ^ 2)).clip(min=0)

        base = self.soils[i][1].uh_base
        q = self.soils[i][1].Qall
        delta_runoff = [q*u(k, base) for k in range(base)]
        self.soils[i][1].runoff += delta_runoff

    def runoff(self, i):
        return self.soils[i][1].runoff
