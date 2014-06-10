#!/bin/env python

"""
Second demonstration script for Scipy2014 build off
iPython notebook
=====================

Takes four processed radar files in CF Radial,
does initial QC and then mesh maps the 3 X-Band systems
and co-maps the C-Band system.
"""

from matplotlib import use
use('agg')
import numpy as np
import pyart
import sys
from matplotlib import pyplot as plt
from matplotlib import rc
from netCDF4 import num2date

if __name__=="__main__":
  print pyart.__version__

  #read data

  base_name = sys.argv[1]
  xnw_name = base_name + "corrected_xnwsapr.nc"
  xsw_name = base_name + "corrected_xswsapr.nc"
  xse_name = base_name + "corrected_xsesapr.nc"
  c_name = base_name + "corrected_csapr.nc"
  xnw_radar = pyart.io.read(xnw_name)
  xsw_radar = pyart.io.read(xsw_name)
  xse_radar = pyart.io.read(xse_name)
  c_radar = pyart.io.read(c_name)

  #QC

  bad_rain = np.where(xnw_radar.fields['rain_rate_A']['data'] < 0.0)
  xnw_radar.fields['rain_rate_A']['data'][bad_rain] = np.ma.masked
  bad_rain = np.where(xsw_radar.fields['rain_rate_A']['data'] < 0.0)
  xsw_radar.fields['rain_rate_A']['data'][bad_rain] = np.ma.masked
  bad_rain = np.where(xse_radar.fields['rain_rate_A']['data'] < 0.0)
  xse_radar.fields['rain_rate_A']['data'][bad_rain] = np.ma.masked
  bad_rain = np.where(xsw_radar.fields['corrected_reflectivity']['data'] < -80.0)
  xsw_radar.fields['corrected_reflectivity']['data'][bad_rain] = np.ma.masked
  bad_rain = np.where(xse_radar.fields['corrected_reflectivity']['data'] < -80.0)
  xse_radar.fields['corrected_reflectivity']['data'][bad_rain] = np.ma.masked
  bad_rain = np.where(xnw_radar.fields['corrected_reflectivity']['data'] < -80.0)
  xnw_radar.fields['corrected_reflectivity']['data'][bad_rain] = np.ma.masked

  #plot a figure

  f, plts = plt.subplots(figsize = [15,4], ncols = 4)
  c_display = pyart.graph.RadarMapDisplay(c_radar)
  plt.sca(plts[0])
  c_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97, title = "C band")
  plt.sca(plts[1])
  xnw_display = pyart.graph.RadarMapDisplay(xnw_radar)
  xnw_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97, title = "NW X band")
  plt.sca(plts[2])
  xsw_display = pyart.graph.RadarMapDisplay(xsw_radar)
  xsw_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97, title = "SW X band")
  plt.sca(plts[3])
  xse_display = pyart.graph.RadarMapDisplay(xse_radar)
  xse_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97, title = "SE X band")
  plt.savefig(base_name + 'four_ppis.png')
  plt.close(f)
