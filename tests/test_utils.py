import utils
import nose


class TestParseNetwork(object):

    @classmethod
    def setup_class(klass):
        print "SETUP!"

    @classmethod
    def teardown_class(klass):
        print "TEAR DOWN!"

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_init(self):
        a = utils.parse_network('test_data/mt_network.geojson')
        nose.tools.assert_equal('test_data/mt_network.geojson', a.fn_vector)

    def test_parse_network(self):
        pass
