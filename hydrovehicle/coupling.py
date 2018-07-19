import numpy as np
import numpy.ma as ma
import econengine as econ
from utils.crop_coefficient import retrieve_crop_coefficient


__all__ = ['HydroEconCoupling']


class HydroEconCoupling(object):
    """Couples the hydrologic and economic models"""

    def __init__(self, routing_obj, water_users_lst):

        self.nodes = routing_obj
        self.water_users = water_users_lst

        self.ma_farms_table = self._build_water_user_matrix()
        self.farm_idx = np.where(self.ma_farms_table[:, 1:])

        self.applied_water_factor = self.ma_farms_table.copy()

    @staticmethod
    def apply_to_all_members(sequence, attrib, *args, **kwargs):
        lst = []
        for obj in sequence:
            try:
                lst.append(getattr(obj, attrib)(*args, **kwargs))
            except TypeError:
                lst.append(getattr(obj, attrib))

        return lst

    def _build_water_user_matrix(self):
        """ Loop through nodes in the network, find farms diverting from it and construct matrix of
         farms associated to each node."""
        nodes = []

        for ids in self.nodes.conn.index:
            li = [ids]
            for farm in self.water_users:
                if farm.get('source_id') == ids:
                    li.append(econ.Farm(**farm))
                else:
                    li.append(None)
            nodes.append(li)
        arr_nodes = ma.array(nodes, mask=[np.array(nodes) == None])

        return arr_nodes

    def simulate_all_users(self, lst_observations):
        # type: (list) -> None

        for obs in lst_observations:
            for farm in self.ma_farms_table[:,1:][self.farm_idx]:
                if obs.get("farm_id") == farm.source_id:
                    farm.simulate(**obs)

        self._calculate_applied_water_factor()

    def total_water_use_per_node(self):
        """Returns an array with total water used from
         each node from all users diverting from it.

         Note that this represents water used, not water diverted.

        :param node: NxM numpy array of WaterUsers. each row is one of hte N nodes
         in the network. The first column are node IDs, which if the subsequent columns
         represents one of the M-1 water users in the system.

        :returns: an array of length N with the total water used from each node.

        """

        wu = self.ma_farms_table.copy()
        wu[:,1:][self.farm_idx] = self.apply_to_all_members(self.ma_farms_table[:,1:][self.farm_idx], "watersim")
        wu[:, 1:][self.farm_idx] = self.apply_to_all_members(wu[:, 1:][self.farm_idx], "sum")
        print np.vstack((wu[:,0], wu[:, 1:].sum(axis=1))).T


    def _calculate_applied_water_factor(self):
        """Calculates a matrix of arrays with the water diversion adjustment factors per crop and farm.

         The factor takes into account the irrigation efficient as well as the length of the crop period
          expressed as the accumulation of crop coefficients. The factor is defined as follows:

          ::

          f:= Sum_t(Kc_t) * Ieff

        The actual daily water diverted (D) to supply water for each crop can then be calculated as:

         ::

         D = Wtot_t * Kc_t / f

        Factor f and the subsequent calculation of D is per crop, so thhis function yields a vector per farm, with one
        f per crop.


        """
        from dateutil import parser
        import datetime

        Kcs = np.vectorize(retrieve_crop_coefficient)

        for i, farm in enumerate(self.ma_farms_table[:,1:][self.farm_idx]):
            try:
                dates = zip(farm.crop_start_date,
                            farm.crop_cover_date,
                            farm.crop_end_date,
                            farm.crop_id)
            except TypeError, e:
                print "Water User %s does not have information on crop planting dates. Did you forget to " \
                      "simulate a scenario?" %farm.name
                exit(-1)

            lst_kc = []
            for s, c, e, cropid, in dates:

                date_array = [(parser.parse(s) + datetime.timedelta(days=x)).strftime("%m/%d/%Y")
                              for x in range(0, (parser.parse(e) - parser.parse(s)).days + 1)]
                lst_kc.append(
                      Kcs(date_array, s, c, e, cropid).sum() * farm.irr_eff * farm.irr
                )
                self.applied_water_factor[:, 1:][self.farm_idx][i] = np.array(lst_kc)

