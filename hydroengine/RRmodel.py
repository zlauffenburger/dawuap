from abc import ABCMeta, abstractmethod


class RRmodel(object):
    """
    This is a base class that defines the basic interface of rainfall runoff models
    """
    __metaclass__ = ABCMeta

    def __init__(self, dt):
        self.dt = dt
        pass

    @abstractmethod
    def runoff(self):
        pass

    @abstractmethod
    def run_time_step(self, *args):
        pass
