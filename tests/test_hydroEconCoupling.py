import numpy as np
import utils
from hydrovehicle.coupling import HydroEconCoupling
import hydroengine as hyd
import econengine as econ
import json
import nose


class TestHydroEconCoupling(object):

    @classmethod
    def setup_class(self):

        print "SETUP!"
        network_file_name = "./test_data/mt_network.shp"
        graph = utils.ParseNetwork(network_file_name)
        self.adj_net = graph.conn_matrix
        self.dt = 86400

        self.mc = hyd.Routing(self.adj_net, self.dt)

        # Open water user object
        with open('./test_data/Farms.json') as json_farms:
            farms = json.load(json_farms)

        # retrieve the list of farms in the json input
        self.lst_farms = farms['farms']

        # Open economic scenario object
        with open('./test_data/Scenario.json') as json_scenario:
            scenario = json.load(json_scenario)

        self.scenario = scenario

        precip = utils.RasterParameterIO('./test_data/precip.nc')
        self.affine = precip.transform
        self.pp_data = precip.array.clip(min=0)[0, :, :]

        self.coupling = HydroEconCoupling(self.mc, self.lst_farms,
                                          self.pp_data, self.affine)

        self.coupling.simulate_all_users(self.scenario)

    @classmethod
    def teardown_class(cls):
        print "TEAR DOWN!"

    def setup(self):
        pass

    def teardown(self):
        pass

    def test__init__(self):
        a = HydroEconCoupling(self.mc, self.lst_farmsm, self.pp_data, self.affine)
        nose.tools.assert_is_instance(a, HydroEconCoupling)
        nose.tools.assert_is_instance(a.nodes, hyd.Routing)
        nose.tools.assert_is_instance(a.water_users, list)

    def test__build_water_user_matrix(self):
        nose.tools.assert_is_instance(self.coupling.farms_table, np.ndarray)
        nose.tools.assert_equal(2, np.count_nonzero(self.coupling.farms_table[:, 1:]))
        np.testing.assert_array_equal((np.array([106, 197]), np.array([1, 0])), self.coupling.farm_idx)

    def test_simulate_all_users(self):

        self.coupling.simulate_all_users(self.scenario)
        for i, users in enumerate(self.coupling.farms_table[:, 1:][self.coupling.farm_idx]):
            np.testing.assert_array_equal(self.scenario[i].get('crop_start_date'),
                                          users.crop_start_date)

    def test__calculate_applied_water_factor(self):
        self.coupling._calculate_applied_water_factor()
        nose.tools.assert_is_instance(self.coupling.applied_water_factor[:, 1:][self.coupling.farm_idx][0], np.ndarray)

    def test_calculate_water_diversion_per_node(self):
        cur_date = "7/02/2014"
        self.coupling.calculate_water_diversion_per_node(cur_date)

    def test__rasterize_water_user_polygons(self):
        self.coupling._rasterize_water_user_polygons('./test_data/Counties.shp', 0, )






