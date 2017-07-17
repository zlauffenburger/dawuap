from __future__ import division
import numpy as np


def hamon_pe(t_avg, latitude, doy):
    """
    Computes potential evapotranspiration using the Hamon method, which depends on average temperature,
    daytime length, and saturation vapor pressure

    :param t_avg: Average temperature, C
    :param latitude: latitude, deg
    :param doy: day of year, starting on jan 1st
    :return: evapotranspiration, mm/day
    """

    p = np.arcsin(0.39795*np.cos(0.2163108 + 2.0 * np.arctan(0.9671396*np.tan(0.00860*(doy-186)))))
    daylength = 24.0 - (24.0/np.pi) * (np.arccos((np.sin(0.8333*np.pi/180.0)
                                                  + np.sin(latitude*np.pi/180.0)
                                                  * np.sin(p))/(np.cos(latitude*np.pi/180.0)*np.cos(p))))
    e_sat = 0.6108*np.exp((17.27 * t_avg) / (237.3 + t_avg))
    pe = (715.5 * daylength * e_sat/24.0)/(t_avg + 273.2)

    return pe
