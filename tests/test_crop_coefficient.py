from __future__ import print_function
from context import utils


class TestCropCoefficient(object):

    def setup(self):
        start = "15/06/2017"
        cover = "10/07/2017"
        end = "10/09/2017"

    def teardown(self):
        pass

    def rest_retrieve_crop_coefficient_rapid_growth(self):
        cur_date = "25/06/217"
        cropid = 8
        utils.retrieve_crop_coefficient(cur_date, self.start, self.cover, self.end, cropid)
