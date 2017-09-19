from context import utils
import filecmp


class TestWriteOutputTimeSeries(object):
    @classmethod
    def setup_class(cls):
        print("setting up class " + cls.__name__)
        TestWriteOutputTimeSeries.network_geojson = '../tests/test_data/mt_network.geojson'
        TestWriteOutputTimeSeries.network_shp = '../tests/test_data/mt_network.shp'
        TestWriteOutputTimeSeries.graph = utils.ParseNetwork(TestWriteOutputTimeSeries.network_geojson)
        TestWriteOutputTimeSeries.adj_net = TestWriteOutputTimeSeries.graph.conn_matrix
        TestWriteOutputTimeSeries.num_links = len(TestWriteOutputTimeSeries.adj_net.index)

    @classmethod
    def teardown_class(cls):
        print("tearing down class " + cls.__name__)

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_write_json(self):

        data = [list(TestWriteOutputTimeSeries.adj_net.index + ts) for ts in range(10)]
        a = utils.WriteOutputTimeSeries(TestWriteOutputTimeSeries.adj_net, "2009-09-30")
        a.write_json(data)
        assert filecmp.cmp('streamflows.json', '../tests/test_data/streamflows.json', shallow=False)
