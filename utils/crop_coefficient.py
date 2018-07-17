import pandas as pd
import numpy as np
from datetime import date
import pkgutil
data = pkgutil.get_data('', 'utils/crop_coefficients.txt')


def retrieve_crop_coefficient(current_date, start_date, cover_date, end_date, crop_id, kc_table="crop_coefficients.txt"):
    "Returns crop coefficient for current_date interpolated from agMet lookup table"

    df_kc = pd.read_table(data, index_col="crop_id")
    current_date = date(current_date)
    start_date = date(start_date)
    cover_date = date(cover_date)
    end_date = date(end_date)
    crop_id = int(crop_id)

    if (current_date > start_date) & (current_date < cover_date):
        frac_growing_season = (start_date - current_date).days / (cover_date - start_date).days
    elif (current_date >= cover_date) & (current_date <= end_date):
        frac_growing_season = (cover_date - current_date).days / (end_date - cover_date).days
    else:
        return 0.0

    fl = df_kc[crop_id, frac_growing_season.floor]
    ceil = df_kc[crop_id, frac_growing_season.ceiling]

    return fl + (ceil - fl) * 0.1 * frac_growing_season


    print df_kc[crop_id]
