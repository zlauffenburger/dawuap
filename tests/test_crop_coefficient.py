from .context import utils
from nose.tools import with_setup


class TestCropCoefficient(object):

    @classmethod
    def setup_class(cls):
        print
        "SETUP!"
        cls.start = "06/15/2014"
        cls.cover = "07/12/2014"
        cls.end = "09/12/2014"

    @classmethod
    def teardown_class(klass):
        print
        "TEAR DOWN!"

    def setup_function(self):
        pass

    def teardown_function(self):
        pass

    #@with_setup(setup_function, teardown_function)
    def rest_retrieve_crop_coefficient(self):
        cur_date = "25/06/2014"
        cropid = 8
        utils.retrieve_crop_coefficient(cur_date, self.start, self.cover, self.end, cropid)

