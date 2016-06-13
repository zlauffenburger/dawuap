from .context import hydroengine
import numpy as np
import nose
import matplotlib.pyplot as plt


class Test_hbv(object):
    @classmethod
    def setup_class(klass):
        print "SETUP!"

    @classmethod
    def teardown_class(klass):
        print "TEAR DOWN!"
        plt.close()

    def TestInitialization(self):
        Q = hydroengine.hbv().excess_precip_to_runoff()
        plt.plot(Q)
        plt.show()

    def Test

