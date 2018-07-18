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

        self.observs = [{
            'farm_id': 198,
            'evapotranspiration': [5., 5.],
            'prices': [5.82, 125],
            'costs': [111.56, 193.95],
            'land_constraint': 100,
            'water_constraint': 100,
            'crop_start_date': ["5/15/2014", "5/15/2014", "5/15/2014", "5/15/2014", "5/15/2014",
                                "5/15/2014", "5/15/2014", "5/15/2014"],
            'crop_cover_date': ["7/02/2014", "7/02/2014", "7/02/2014", "7/02/2014", "7/02/2014",
                                "7/02/2014", "7/02/2014", "7/02/2014"],
            'crop_end_date': ["8/25/2014", "8/25/2014", "8/25/2014", "8/25/2014", "8/25/2014",
                              "8/25/2014", "8/25/2014", "8/25/2014"],
        },
        {
            'farm_id': 107,
            'evapotranspiration': [5., 5.],
            'prices': [5.82, 125],
            'costs': [111.56, 193.95],
            'land_constraint': 100,
            'water_constraint': 100,
            'crop_start_date': ["5/15/2014", "5/15/2014", "5/15/2014", "5/15/2014", "5/15/2014",
                                "5/15/2014", "5/15/2014", "5/15/2014"],
            'crop_cover_date': ["7/02/2014", "7/02/2014", "7/02/2014", "7/02/2014", "7/02/2014",
                                "7/02/2014", "7/02/2014", "7/02/2014"],
            'crop_end_date': ["8/25/2014", "8/25/2014", "8/25/2014", "8/25/2014", "8/25/2014",
                              "8/25/2014", "8/25/2014", "8/25/2014"],
        }]

        self.coupling = HydroEconCoupling(self.mc, self.lst_farms)

        #self.mat_farms = _build_water_user_matrix()
        #self.farm_idx = np.where(isinstance(econ.WaterUser, self.mat_farms))

    @classmethod
    def teardown_class(cls):
        print "TEAR DOWN!"

    def setup(self):
        pass

    def teardown(self):
        pass

    def test__init__(self):
        a = HydroEconCoupling(self.mc, self.lst_farms)
        nose.tools.assert_is_instance(a, HydroEconCoupling)
        nose.tools.assert_is_instance(a.nodes, hyd.Routing)
        nose.tools.assert_is_instance(a.water_users, list)

    def test__build_water_user_matrix(self):
        nose.tools.assert_is_instance(self.coupling.ma_farms_table, np.ndarray)
        nose.tools.assert_equal(2, np.count_nonzero(self.coupling.ma_farms_table[:, 1:]))
        np.testing.assert_array_equal((np.array([106, 197]), np.array([1, 0])), self.coupling.farm_idx)

    def test_simulate_all_users(self):

        self.coupling.simulate_all_users(self.observs)
        for users in self.coupling.water_users:
            print users
