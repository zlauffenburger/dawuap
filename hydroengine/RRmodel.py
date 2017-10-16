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

    @abstractmethod
    def pickle_current_states(self):
        """Virtual method to be overriden by a derived class.

        Derived class must implement a function to pickle all state variables and permit to restart model from
        pickled states

        Function must take no arguments. No specification for return values
        """
        pass

    @abstractmethod
    def unpickle_current_states(self):
        """Virtual method to be overriden by a derived class.

        Derived class must implement a function to read all pickled state variables written by
        pickle_current_states() and permit to restart model these states

        Function must take no arguments. No specification for return values
        """
        pass

    @abstractmethod
    def write_current_states(self, current_ts, ext, callback):
        """Virtual method to be overriden by a derived class.

        No specification for return values.
        Raises IOError exception if writing fails

        This function writes to disk the model state for any time step in a format that can be used for
        plotting, visualization.

        :param current_ts: Current timestep, string
        :param ext: extension to be attached at the end of the output files, str
        :param callback: Callback function to write states to disk

        """
        pass



