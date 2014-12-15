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

import sys
import os

path, file = os.path.split(os.path.abspath(__file__))
sys.path.append(path)

from blender import *

if __name__ == "__main__": 
  f, world_vec, species_size = init_coord_file(filename)
  for i in range(start_frame):
    for j in range(species_size):
      time, c = load_coords(f)
  delete_home_scenes()
  spheres = set_new_scene(world_vec, species_size)
  i = 0
  for j in range(species_size):
    time, c = load_coords(f)
    if len(c):
      for k in range(0, int(len(c)/3)):
        print_sphere((c[k*3],c[k*3+1],c[k*3+2]), spheres[j])
  update_time(time)
  save('frame.blend')
  print('Saved frame.blend file')

