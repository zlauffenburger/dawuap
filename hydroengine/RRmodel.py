from abc import ABCMeta, abstractmethod


class RRmodel(object):
    """
    Base class defining the basic interface of rainfall runoff models.

    The actual rainfall runoff models are implemented in a derived class. This class
    must override the runoff() and run_time_step() methods
    """
    __metaclass__ = ABCMeta

    def __init__(self, dt):
        self.dt = dt
        pass

    @abstractmethod
    def runoff(self):
        """
        Virtual method to be overriden by a derived class.

        :return: Derived class must return a list of streamflows, with one element per node in the domain.

        """
        pass

    @abstractmethod
    def run_time_step(self, *args):
        """
        Virtual method to be overriden by a derived class.
        Derived class must implement the run_time_step method(). This method must run the
        model forward one time step.

        No specifications for input arguments and return values
        """
        pass
