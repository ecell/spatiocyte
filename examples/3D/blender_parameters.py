#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#        This file is part of Spatiocyte particle simulator package
#
#             Copyright (C) 2009-2014 RIKEN, Keio University
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# Spatiocyte is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# Spatiocyte is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public
# License along with Spatiocyte -- see the file COPYING.
# If not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
# 
#END_HEADER
#
# written by Satya Arjunan <satya.arjunan@gmail.com>
# Spatiocyte, RIKEN Quantitative Biology Center, Osaka, Japan
#

import math

start_frame = 1
end_frame = 350
resolution_x = 1920
resolution_y = 1080
resolution_percentage = 50
render_samples = 20
lamp_shadow_size = 0.1
lamp_strength = 3
plane_scale = 5
background_strength = 0.3
species_radius_delta = 0.01
visible_planes = [1, 1, 1, 0, 0, 0]
camera_rotation = (63.18*math.pi/180.0,0*math.pi/180.0,136.6*math.pi/180.0)
camera_location = (65.35, 69.33, 50.31)
time_location = (46.43, 15.87, 37.37)
lamp_location = (7.88, 37.38, 73.27)
lamp_rotation = (-7.82*math.pi/180.0,0.69*math.pi/180.0,92.31*math.pi/180.0)
plane_disp = [1.0, 1.25, 1.5]
#plane_disp = [0.5, 0, 0.5]
#set True if using GPU device to render
GPU_device = False
#for GPU device, set tile_x = 512, tile_y = 512
tile_x = 256
tile_y = 256
plane_material_name = 'White'
time_material_name = 'White'
filename = 'CoordinateLog.csv'
species_material_names = ['DarkRed_glossy','DarkBlue_glossy','Green_glossy', 'BrightGreen_glossy','Magenta_glossy','DarkOrange_glossy','Cyan_glossy']

#start_frame = 5
#end_frame = 350
#resolution_x = 1920
#resolution_y = 1920
#resolution_percentage = 100
#render_samples = 1000
#lamp_shadow_size = 0.1
#lamp_strength = 3
#plane_scale = 5
#background_strength = 0.3
#species_radius_delta = 0.01
#visible_planes = [1, 1, 1, 0, 0, 0]
#camera_rotation = (63.18*math.pi/180.0,0*math.pi/180.0,136.6*math.pi/180.0)
#camera_location = (65.35, 69.33, 50.31)
#time_location = (46.43, 15.87, 37.37)
#lamp_location = (7.88, 37.38, 73.27)
#lamp_rotation = (-7.82*math.pi/180.0,0.69*math.pi/180.0,92.31*math.pi/180.0)
#plane_disp = [1.0, 1.25, 1.5]
##plane_disp = [0.5, 0, 0.5]
##set True if using GPU device to render
#GPU_device = True
##for GPU device, set tile_x = 512, tile_y = 512
#tile_x = 256
#tile_y = 256
#plane_material_name = 'White'
#time_material_name = 'Black'
#filename = 'CoordinateLog.csv'
#species_material_names = ['Yellow_glossy','LightBlue_glossy','DarkOrange_glossy','Cyan_glossy']
