import os

from paraview.simple import*
LoadPlugin('/home/chin/root/lib/paraview-4.1/libGMVReader.so',False,globals())

#for i in range(0,9):
#  for j in range(0,9):
#    for k in range(0,9):
camera = GetActiveCamera()
#camera.Elevation(280)
view = GetActiveView()
view.Background = [1,1,1]
view.ViewSize = [900,400]
#camera.Roll(45)
for i in range(323):
  view.CameraViewAngle = 25
  filename = 'migration.%03d'%(i,)
  print filename
  os.system("xgmvsw " + filename)
#  xgmvsw filename
  reader = OpenDataFile("/home/chin/wrk/temp9/spatiocyte/samples/migration/rif.gmv")
  #view=GetRenderView()
  #SetActiveView(view)
  #view.Elevation(45)
  dp = GetDisplayProperties(reader)
  dp.Representation = 'Surface With Edges'
  reader.GetProperty('ImportPolygons').SetData(0)
  reader.GetProperty('IgnorestoredtimestepvaluesPROBTIMEkeyword').SetData(0)
  a = reader.PointData[4] #thetaN
  #4,5,6,7,9-PI3Km-PTENm
  #a = reader.PointData[3] #pip2
  print a.GetName()
  r = a.GetRange()
  dp.LookupTable = AssignLookupTable(a,'Cool to Warm')
  dp.ColorArrayName = a.GetName()
  dp.LookupTable.ColorSpace = 'Diverging'
  #view.CameraViewAngle = 25
  #view.CameraFocalPoint = [0,0,0]
  #CameraViewAngle = 90
  #ResetCamera()
  Show()
  Render()
  #camera.Zoom(2)
  WriteImage('mechanocyte%03d.png'%(i,))
  #ResetCamera()
  Delete()


#LoadPlugin('/home/chin/root/lib/paraview-4.1/libGMVReader.so',False,globals())
#MakeCoolToWarm(r[0], r[1])
#dp.ColorAttributeType = 'POINT_DATA'
#dpr = dp.GetRange()
#dp.LookupTable = MakeBlueToRedLT(0.0, 0.5)
#dp.ColorAttributeType = 'POINT_DATA'
#dp.ColorArrayName = 'THETAEQ'
#dir(reader)
#reader.ListProperties()
#dir(reader.GetProperty('ImportPolygons'))
#reader.GetProperty('ImportPolygons').GetData()
#reader.GetProperty('ImportPolygons').SetData(0)
#reader.GetProperty('IgnorestoredtimestepvaluesPROBTIMEkeyword').SetData(0)
#Sphere()
#from paraview.servermanager import*
#GetDisplayProperties(reader).Representation = 'Surface with Edges'
#Show()
#Render()


