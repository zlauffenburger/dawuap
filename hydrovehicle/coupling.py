import numpy as np
import numpy.ma as ma
import econengine as econ
from utils.crop_coefficient import retrieve_crop_coefficient

__all__ = ['HydroEconCoupling']


class HydroEconCoupling(object):
    """Couples the hydrologic and economic models"""

    def __init__(self, routing_obj, water_users):

        self.nodes = routing_obj
        self.water_users = water_users

        self.ma_farms_table = self._build_water_user_matrix()
        self.farm_idx = np.where(self.ma_farms_table[:, 1:])

    @staticmethod
    def apply_to_all_members(sequence, attrib, *args, **kwargs):
        lst = []
        for obj in sequence:
            lst.append(getattr(obj, attrib)(*args, **kwargs))
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
                if obs.get("farm_id") == farm.id:
                    farm.simulate(**obs)

    def node_total_water_use(node):
        """Returns an array with total water used from
         each node from all users diverting from it.

         Note that this represents water used, not water diverted.

        :param node: NxM numpy array of WaterUsers. each row is one of hte N nodes
         in the network. The first column are node IDs, which if the subsequent columns
         represents one of the M-1 water users in the system.

        :returns: an array of length N with the total water used from each node.

        """
        wu = map(lambda x: x.watersim.sum() if isinstance(x, econfuncs.WaterUser)
                 else 0., node[1:])
        return np.append(node[0], sum(wu))

    def calculate_applied_water_per_crop(self):
        """Returns a matrix of arrays of water used per crop for each water user and
         water diversion node.

         Note that the values represents water used, not water diverted.

        :param node: NxM numpy array of WaterUsers. each row is one of hte N nodes
         in the network. The first column are node IDs, the subsequent columns
         are each of the M-1 WaterUser objects represented in the system.

        :returns: an array of length N with arrays of water used per crop
        grown by the water user from each node.

        """

        pass

        # aw = self.ma_farms_table.copy()
        #
        # f = []
        # for farm, idx in self.ma_farms_table[self.farm_idx]:
        #     tot_kc = []
        #     farm.watersim *
        # wu = map(lambda x: x.watersim if isinstance(x, econ.WaterUser)
        #          else 0., node[1:])
        # retrieve_crop_coefficient()
        # return np.append(node[0], sum(wu))
