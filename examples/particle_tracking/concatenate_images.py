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


import os
import glob

left = sorted(glob.glob('image????.png'))
right = sorted(glob.glob('image???????.png'))

left_interval = 0.01
right_interval = 0.5
right_frame_cnt = int(right_interval/left_interval)

start_frame = int(5000/right_frame_cnt)
end_frame = int(len(left)/right_frame_cnt)

for i in range(start_frame, end_frame):
#for i in range(int(1000/right_frame_cnt)):
  for j in range(right_frame_cnt):
    idx = i*right_frame_cnt+j
    string = 'convert %s %s +append concatenated%04d.png' %(left[idx],right[i],idx)
    os.system(string)



