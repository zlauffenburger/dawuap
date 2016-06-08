from .context import hydroengine
import numpy as np
import nose
import matplotlib.pyplot as plt


class TestInitialization(object):
    def TestInitialization(self):
        Q = hydroengine.hbv().excess_precip_to_runoff()
        plt.plot(Q)
        plt.show()
