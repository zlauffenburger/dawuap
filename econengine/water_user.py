from abc import ABCMeta, abstractmethod


class WaterUser(object):
    """Abstract base class for water users

    The actual water users are implemented in a derived class. This class
        must override the simulate(calibrate) and calibrate() methods
    """

    __metaclass__ = ABCMeta

    def __init__(self, identifier, source_id, name):
        self.id = identifier
        self.name = name
        self.source_id = source_id

    @abstractmethod
    def simulate(self):
        """
        Virtual method to be overriden by a derived class.

        :return: Derived class must return a list of water demands

        """
        pass