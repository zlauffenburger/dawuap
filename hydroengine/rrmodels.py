from __future__ import division
from .RRmodel import RRmodel
import cPickle as pickle
import numpy as np
import rasterstats as rst
import rasterio as rio

class Soil(object):
    def __init__(self, ur=0., lr=0., q0=0., q1=0., q2=0., qall=0., base=25):

        self.upper_reservoir = ur
        self.lower_reservoir = lr

        self.Q0 = q0
        self.Q1 = q1
        self.Q2 = q2
        self.Qall = qall

        self.uh_base = base
        self._runoff = np.zeros(base)


    @property
    def runoff(self, value):
        self._runoff = value

    @runoff.getter
    def runoff(self):
        return self._runoff

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


class HBV(RRmodel):
    """Implementation of the classic HVB rainfall-runoff model

    """

    def __init__(self, dt, swe_o, pond_o, sm_o, soils_o=[], **params):

        # Geomtry information
        #self.pixel_area = params['image_res']*params['image_res']
        #self.catchment_area = params['catch_area']
        super(HBV, self).__init__(dt)

        # snow paramters
        self.t_thres = params['pp_temp_thres']
        self.ddf = params['ddf'] * self.dt / 86400

        # distributed soil paramters
        self.fcap = params['soil_max_wat']
        self.beta = params['soil_beta']
        self.lp = params['aet_lp_param']

        # base of the unit hydrograph
        #self.base = params['p_base']

        # self.hl1 = params['storage_parameter_1']
        # self.hl2 = params['storage_parameter_2']
        #
        # self.ck0 = params['surface_conductance']
        # self.ck1 = params['mid_layer_conductance']
        # self.ck2 = params['lower_layer_conductance']

        self.soils = soils_o  # Soil reinitialization file. Dictionary with geo_json objects

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
        self.pond = np.zeros_like(self.swe)
        self.melt = np.zeros_like(self.swe)

        ind_allswe = np.less_equal(t_max, self.t_thres)
        ind_allrain = np.greater(t_min, self.t_thres)
        ind_mixed = np.logical_not(np.logical_or(ind_allswe, ind_allrain))

        swe[ind_allswe] = incid_precip[ind_allswe]
        rain[ind_allrain] = incid_precip[ind_allrain]
        swe[ind_mixed] = (incid_precip * (np.array((self.t_thres - t_min) / (t_max - t_min).clip(min=0.01)).clip(min=0)))[ind_mixed]
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
        self.aet = (pot_et.clip(min=0) * (self.sm/(self.fcap*self.lp)).clip(min=0.0, max=1.0))

        # For cell where soil moisture is already at capacity, runoff is all the ponded water plus excess water in soil
        ind_sm_geq_fcap = np.greater_equal(self.sm, self.fcap)
        self.ovlnd_flow[ind_sm_geq_fcap] = self.pond[ind_sm_geq_fcap] + (self.sm - self.fcap)[ind_sm_geq_fcap]
        self.sm[ind_sm_geq_fcap] = self.fcap[ind_sm_geq_fcap]

        # in all other cells calculate the portion of ponded water that goes into the soil storage
        ind_sm_lt_fcap = np.logical_not(ind_sm_geq_fcap)
        self.delta_sm[ind_sm_lt_fcap] = (self.pond * (1 - np.power((self.sm/self.fcap), self.beta).clip(max=1.0, min=0)))[ind_sm_lt_fcap]
        self.sm[ind_sm_lt_fcap] += self.delta_sm[ind_sm_lt_fcap]
        self.ovlnd_flow[ind_sm_lt_fcap] = (self.pond - self.delta_sm)[[ind_sm_lt_fcap]]

        # Check if cell exceed storage capacity after adding delta_sm
        ind_sm_geq_fcap = np.greater_equal(self.sm, self.fcap)
        self.ovlnd_flow[ind_sm_geq_fcap] += (self.sm - self.fcap)[ind_sm_geq_fcap]
        self.sm[ind_sm_geq_fcap] = self.fcap[ind_sm_geq_fcap]

        # if there is sufficient soil moisture to satisfy aet, reduce sm
        ind_sm_gt_aet = np.greater(self.sm, self.aet)
        self.sm[ind_sm_gt_aet] -= self.aet[ind_sm_gt_aet]
        # otherwise take all water storage and limit aet
        ind_sm_leq_aet = np.logical_not(ind_sm_gt_aet)
        self.aet[ind_sm_leq_aet] = self.sm[ind_sm_leq_aet]
        self.sm[ind_sm_leq_aet] = 0.0

    def precipitation_excess(self, shp_wtshds, affine=None, stats=['mean'], **kwargs):

        affine = affine
        # builds a geojson object with required statistics for each catchment
        #print "Calculating zonal statistics for each HRU..."
        nodata = kwargs['nodata'] if kwargs.has_key('nodata') else -32767
        stw1 = rst.zonal_stats(shp_wtshds, self.ovlnd_flow, nodata=nodata, affine=affine, geojson_out=True,
                                    prefix='runoff_', stats=stats)

        soil_layers = {}
        soils = []

        for i in range(len(stw1)):
            #print "processing catchment ", i
            props = stw1[i]['properties']

            if not self.soils: # if there is no dictionary from a previous time step, create a new soil object
                soil_layers = Soil(base=props['hbv_pbase'])
            else:
                soil_layers = self.soils[i][1]

            #props = stw1[i]['properties']

            # conversion stuff
            hbv_ck0 = self.dt / props['hbv_ck0'] / 86400
            hbv_ck1 = self.dt / props['hbv_ck1'] / 86400
            hbv_ck2 = self.dt / props['hbv_ck2'] / 86400
            # If a catchment is not covered by the climate layer, zonal_stats results
            # in None and raises an exception
            if props['runoff_mean'] is None:
                props['runoff_mean'] = 0
                # print("WARNING: Catchment " + str(props['Subbasin']) + "  is probably"
                #       " outside the region covered by the climate grids")
            soil_layers.upper_reservoir += props['runoff_mean']
            if soil_layers.upper_reservoir > props['hbv_hl1']:
                soil_layers.Q0 = (soil_layers.upper_reservoir - props['hbv_hl1']) * hbv_ck0
                soil_layers.upper_reservoir -= soil_layers.Q0
            else:
                soil_layers.Q0 = 0.0

            if soil_layers.upper_reservoir > 0.0:
                soil_layers.Q1 = soil_layers.upper_reservoir * hbv_ck1
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
                soil_layers.Q2 = soil_layers.lower_reservoir * hbv_ck2
                soil_layers.lower_reservoir -= soil_layers.Q2
            else:
                soil_layers.Q2 = 0.0

            soil_layers.Qall = soil_layers.Q0 + soil_layers.Q1 + soil_layers.Q2
            self._calculate_runoff(soil_layers)
            soils.append((stw1[i]['id'], soil_layers))


        self.soils = soils
        pickle.dump(self.soils, open("soils.pickled", "wb"))

    def _calculate_runoff(self, soil_layer):
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

        runoff = np.roll(soil_layer._runoff, -1, axis=0)
        runoff[-1] = 0
        base = soil_layer.uh_base
        q = soil_layer.Qall
        delta_runoff = [q*u(k+1, base) for k in range(base)]
        soil_layer._runoff = runoff + delta_runoff

    def pickle_current_states(self):

        # save surface states
        pickle.dump(self.sm, open("sm.pickled", "wb"))
        pickle.dump(self.swe, open("swe.pickled", "wb"))
        pickle.dump(self.pond, open("pond.pickled", "wb"))

        # save soil states
        pickle.dump(self.soils, open("soils.pickled", "wb"))

    def unpickle_current_states(self):
        pass

    def write_current_states(self, current_ts, ext, callback):
        """Saves state and diagnostic variables to disk. Filename of written state includes current_ts

        The object writing operation is done by a callback function with prototype
        type(string, object) -> None

            ```callback(string, object)```

        Function raises IOError if fails to write

        :param current_ts: string with current time stamp
        :param ext: extension to be appended at the end of the output filenames
        :param callback: callback function to handle writing to disk


        """
        try:
            callback("swe_" + str(current_ts) + "." + str(ext).strip('.'), self.swe)
            callback("pond_" + str(current_ts) + "." + str(ext).strip('.'), self.pond)
            callback("melt_" + str(current_ts) + "." + str(ext).strip('.'), self.melt)
            callback("aet_" + str(current_ts) + "." + str(ext).strip('.'), self.aet)
        except IOError as e:
            print e.message
            raise e


    @property
    def runoff(self):
        return [x[1].runoff[0] for x in self.soils]

    def run_time_step(self, incid_precip, t_max, t_min, pot_et, shp_wtshds, affine=None, stats=['mean'], **kwargs):

        self.snow_pack(incid_precip, t_max, t_min)
        self.soil_processes(pot_et)
        self.precipitation_excess(shp_wtshds, affine, stats, **kwargs)

        return self.runoff

