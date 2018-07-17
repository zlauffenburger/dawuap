import numpy as np
from econengine import econfuncs
from utils.crop_coefficient import retrieve_crop_coefficient

__all__ = ['node_total_water_use']


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


def calculate_applied_water(node):
    """Returns a matrix of arrays of water used per crop for each water user and
     water diversion node.

     Note that the values represents water used, not water diverted.

    :param node: NxM numpy array of WaterUsers. each row is one of hte N nodes
     in the network. The first column are node IDs, the subsequent columns
     are each of the M-1 WaterUser objects represented in the system.

    :returns: an array of length N with arrays of water used per crop
    grown by the water user from each node.

    """
    wu = map(lambda x: x.watersim if isinstance(x, econfuncs.WaterUser)
             else 0., node[1:])
    retrieve_crop_coefficient()
    return np.append(node[0], sum(wu))
