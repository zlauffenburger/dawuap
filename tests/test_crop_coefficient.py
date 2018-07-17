from .context import utils
import nose

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

    def test_retrieve_crop_coefficient_fast_growth(self):
        cur_date = "06/25/2014"
        cropid = 15 # Peas
        res = utils.retrieve_crop_coefficient(cur_date, self.start, self.cover, self.end, cropid)
        nose.tools.assert_true(0.29 <= res <= 0.37)

    def test_retrieve_crop_coefficient_early(self):
        cur_date = "04/25/2014"
        cropid = 15  # Peas
        res = utils.retrieve_crop_coefficient(cur_date, self.start, self.cover, self.end, cropid)
        nose.tools.assert_equal(res, 0.0)

    def test_retrieve_crop_coefficient_late(self):
        cur_date = "10/25/2014"
        cropid = 15  # Peas
        res = utils.retrieve_crop_coefficient(cur_date, self.start, self.cover, self.end, cropid)
        nose.tools.assert_equal(res, 0.0)

    def test_retrieve_crop_coefficient_peak(self):
        cur_date = "9/11/2014"
        cropid = 15  # Peas
        res = utils.retrieve_crop_coefficient(cur_date, self.start, self.cover, self.end, cropid)
        nose.tools.assert_true(0.17 <= res <= 0.22)

