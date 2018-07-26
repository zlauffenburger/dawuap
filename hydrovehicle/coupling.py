import numpy as np
import econengine as econ
import utils
from utils.crop_coefficient import retrieve_crop_coefficient
from dateutil import parser
import datetime
import rasterio as rio
from rasterio.features import rasterize

__all__ = ['HydroEconCoupling']


class HydroEconCoupling(object):
    """Couples the hydrologic and economic models"""

    def __init__(self, routing_obj, water_users_lst, precip_arr, transform):

        self.nodes = routing_obj
        self.water_users = water_users_lst

        self.farms_table = self._build_water_user_matrix()
        self.farm_idx = np.where(self.farms_table[:, 1:])

        self.applied_water_factor = np.zeros_like(self.farms_table)

        self.array_supplemental_irrigation = np.zeros_like(precip_arr)

        self.transform = transform

    @staticmethod
    def apply_to_all_members(sequence, attrib, *args, **kwargs):
        """Returns a list with the results of applying attrib to a sequence of objects. Parameter attrib is
        a string with the method name and can be an attribute obj.attrib or a method obj.attrib(*args, **kwargs)"""

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
                    li.append(0)
            nodes.append(li)
        arr_nodes = np.array(nodes)

        return arr_nodes

    def simulate_all_users(self, lst_scenarios):
        # type: (list) -> None

        for obs in lst_scenarios:
            for farm in self.farms_table[:, 1:][self.farm_idx]:
                if obs.get("farm_id") == farm.source_id:
                    farm.simulate(**obs)

        self._calculate_applied_water_factor()

    def calculate_water_diversion_per_node(self, date):
        """Returns a """

        # obtain vector of crop coefficients
        vect_retrieve_kcs = np.vectorize(retrieve_crop_coefficient, excluded=['current_date'])

        # Obtains water simulated per crop
        s = self.apply_to_all_members(self.farms_table[:, 1:][self.farm_idx], "crop_start_date")
        c = self.apply_to_all_members(self.farms_table[:, 1:][self.farm_idx], "crop_cover_date")
        e = self.apply_to_all_members(self.farms_table[:, 1:][self.farm_idx], "crop_end_date")
        cropid = self.apply_to_all_members(self.farms_table[:, 1:][self.farm_idx], "crop_id")

        current_kcs = vect_retrieve_kcs(date, s, c, e, cropid)

        # Obtain water simulated per crop
        Xw = np.vstack(
            self.apply_to_all_members(
                self.farms_table[:, 1:][self.farm_idx], "watersim"
            )
            )

        # Obtain applied water factor
        f = np.vstack(self.applied_water_factor[:, 1:][self.farm_idx])

        # Calculated water diverted for each crop and farm

        # diversions per per farm, crop and node
        d = Xw * np.divide(current_kcs, f, where=f!=0)
        D = self.farms_table.copy()
        D[:, 1:][self.farm_idx] = tuple(d)

        # Total diversions per node
        # First sum all waterdiverted per crop in each farm
        dtot = [fm.sum() for fm in d]
        Dtot = self.farms_table.copy()
        Dtot[:, 1:][self.farm_idx] = dtot
        Dtot = np.vstack((Dtot[:, 0], Dtot[:, 1:].sum(axis=1)))
        return Dtot, D

    def _calculate_applied_water_factor(self):
        """Sets member variable ``applied_water_factor``, a masked matrix of arrays with the water diversion
        adjustment factors per crop and farm.

         The factor takes into account the irrigation efficient as well as the length of the crop period
          expressed as the accumulation of crop coefficients. The factor is defined as follows:

          ::

          f:= Sum_t(Kc_t) * Ieff

        The actual daily water diverted (D) to supply water for each crop can then be calculated as:

         ::

         D = Wtot_t * Kc_t / f

        Factor f and the subsequent calculation of D is calculated per crop. Thus, th function yields a
        vector per farm, with one f per crop.


        """

        Kcs = np.vectorize(retrieve_crop_coefficient)

        lst_kc = []
        for i, farm in enumerate(self.farms_table[:, 1:][self.farm_idx]):
            try:
                dates = zip(farm.crop_start_date,
                            farm.crop_cover_date,
                            farm.crop_end_date,
                            farm.crop_id,
                            farm.irr_eff,
                            farm.irr)
            except TypeError, e:
                print "Water User %s does not have information on crop planting dates. Did you forget to " \
                      "simulate a scenario?" %farm.name
                exit(-1)

            lst = []
            for s, c, e, cropid, i_eff, i_mask, in dates:

                date_array = [(parser.parse(s) + datetime.timedelta(days=x)).strftime("%m/%d/%Y")
                              for x in range(0, (parser.parse(e) - parser.parse(s)).days + 1)]
                lst.append(
                      Kcs(date_array, s, c, e, cropid).sum() * i_eff * i_mask
                )
            lst_kc.append(np.array(lst))

        self.applied_water_factor[:, 1:][self.farm_idx] = lst_kc

    def _rasterize_water_user_polygons(self, fn_water_user_shapes, fill):

        shapes = utils.VectorParameterIO(fn_water_user_shapes).read_features()

        t = self.array_supplemental_irrigation = \
           rasterize(shapes,
                     self.array_supplemental_irrigation.shape,
                     fill=fill,
                     transform=self.transform)
        return t


