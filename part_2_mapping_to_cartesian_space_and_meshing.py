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

  print "Correcting"

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
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97,
                         title = "C band")
  plt.sca(plts[1])
  xnw_display = pyart.graph.RadarMapDisplay(xnw_radar)
  xnw_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97,
                         title = "NW X band")
  plt.sca(plts[2])
  xsw_display = pyart.graph.RadarMapDisplay(xsw_radar)
  xsw_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97,
                         title = "SW X band")
  plt.sca(plts[3])
  xse_display = pyart.graph.RadarMapDisplay(xse_radar)
  xse_display.plot_ppi_map('rain_rate_A', 0, vmin=0, vmax=150, resolution = 'l',
                         min_lat=36, max_lat=37, min_lon=-98, max_lon=-97,
                         title = "SE X band")
  plt.savefig(base_name + 'four_ppis.png')
  plt.close(f)

  #Map X-Bands to a single cartesian grid

  print "Mapping"

  mesh_mapped_x = pyart.map.grid_from_radars((xnw_radar,xsw_radar,xse_radar),
                                        grid_shape=(35, 201, 201),
                                        grid_limits=((0, 17000), (-50000, 40000),
                                        (-60000, 40000)),
                                        grid_origin = (36.57861, -97.363611),
                                        fields=xnw_radar.fields.keys(),
                                        refl_field='corrected_reflectivity',
                                        max_refl=100., copy_field_data=False)

  #Save out

  print "Saving"

  pyart.io.write_grid('/data/scipy2014/x_mesh2.nc', mesh_mapped_x)

  #Make nice graphics

  display = pyart.graph.GridMapDisplay(mesh_mapped_x)
  # create the figure
  font = {'size': 16}
  rc('font', **font)
  fig = plt.figure(figsize=[15, 8])

  # panel sizes
  map_panel_axes = [0.05, 0.05, .4, .80]
  x_cut_panel_axes = [0.55, 0.10, .4, .30]
  y_cut_panel_axes = [0.55, 0.50, .4, .30]
  colorbar_panel_axes = [0.05, 0.90, .4, .03]

  # parameters
  level = 3
  vmin = 0
  vmax = 64
  lat = 36.57861
  lon = -97.363611

  # panel 1, basemap, radar reflectivity a
  ax1 = fig.add_axes(map_panel_axes)
  display.plot_basemap()
  display.plot_grid('corrected_reflectivity', level=level, vmin=vmin, vmax=vmax)
  display.plot_crosshairs(lon=lon, lat=lat)

  # colorbar
  cbax = fig.add_axes(colorbar_panel_axes)
  display.plot_colorbar(cax=cbax)

  # panel 2, longitude slice.
  ax2 = fig.add_axes(x_cut_panel_axes)
  display.plot_longitude_slice('corrected_reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax)
  ax2.set_xlabel('Distance from SGP CF (km)')

  # panel 3, latitude slice
  ax3 = fig.add_axes(y_cut_panel_axes)
  display.plot_latitude_slice('corrected_reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax)

  # add a title
  slc_height = mesh_mapped_x.axes['z_disp']['data'][level]
  dts = num2date(mesh_mapped_x.axes['time']['data'], mesh_mapped_x.axes['time']['units'])
  datestr = dts[0].strftime('%H:%M Z on %Y-%m-%d')
  title = 'Sliced at ' + str(slc_height) + ' meters at ' + datestr
  fig.text(0.5, 0.9, title)
  plt.savefig(base_name+"reflectivity_x_three_panel.png")

  display = pyart.graph.GridMapDisplay(mesh_mapped_x)
  # create the figure
  font = {'size': 16}
  rc('font', **font)
  fig = plt.figure(figsize=[15, 8])

  # panel sizes
  map_panel_axes = [0.05, 0.05, .4, .80]
  x_cut_panel_axes = [0.55, 0.10, .4, .30]
  y_cut_panel_axes = [0.55, 0.50, .4, .30]
  colorbar_panel_axes = [0.05, 0.90, .4, .03]

  # parameters
  level = 2
  vmin = -80
  vmax = 80
  lat = 36.57861
  lon = -97.763611

  # panel 1, basemap, radar reflectivity a
  ax1 = fig.add_axes(map_panel_axes)
  display.plot_basemap()
  display.plot_grid('rain_rate_A', level=level, vmin=vmin, vmax=vmax)
  display.plot_crosshairs(lon=lon, lat=lat)

  # colorbar
  cbax = fig.add_axes(colorbar_panel_axes)
  display.plot_colorbar(cax=cbax)

  # panel 2, longitude slice.
  ax2 = fig.add_axes(x_cut_panel_axes)
  display.plot_longitude_slice('rain_rate_A', lon=lon, lat=lat, vmin=vmin, vmax=vmax)
  ax2.set_xlabel('Distance from SGP CF (km)')

  # panel 3, latitude slice
  ax3 = fig.add_axes(y_cut_panel_axes)
  display.plot_latitude_slice('rain_rate_A', lon=lon, lat=lat, vmin=vmin, vmax=vmax)

  # add a title
  slc_height = mesh_mapped_x.axes['z_disp']['data'][level]
  dts = num2date(mesh_mapped_x.axes['time']['data'], mesh_mapped_x.axes['time']['units'])
  datestr = dts[0].strftime('%H:%M Z on %Y-%m-%d')
  title = 'Sliced at ' + str(slc_height) + ' meters at ' + datestr
  fig.text(0.5, 0.9, title)
  plt.savefig(base_name+"rain_rate_x_three_panel.png")
