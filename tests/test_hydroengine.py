from .context import hydroengine
import numpy as np
import nose


conn = np.zeros((7,7))
conn[0,2] = 1
conn[1,2] = 1
conn[2,3] = 1
conn[3,4] = 1
conn[3,5] = 1
conn[4,5] = 1
conn[5,6] = 1

class TestInitialization(object):

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
        a = hydroengine.routing(conn, 0.5)
        np.testing.assert_allclose(conn, a.conn)
        nose.tools.assert_equal(a.n, np.shape(conn)[0])


class TestRouting(object):
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

    def test_muskingum_routing(self):
        a = hydroengine.routing(conn, 0.5)
        K =  np.ones(a.n)*0.7
        K[2] = 0.6
        e = np.ones_like(K)*0.42
        Qk = np.ones_like(K)*0.2
        #a.muskingum_routing(Qk,K,e, 0)
        Qk1 = a.muskingum_routing(Qk, K, e, 0, 0)
        print Qk1






