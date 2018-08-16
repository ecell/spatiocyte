import os

if __name__ == '__main__':
  #print("egfrd");
  #os.chdir("egfrd")
  #os.system("python run_all.py")
  #os.system("python run_all_dense.py")
  #print("fastbd");
  #os.chdir("../fastbd")
  #os.system("python run_all.py")
  #os.system("python run_all_dillute.py")
  print("smoldyn");
  os.chdir("../smoldyn")
  os.system("python run_all.py")
  os.system("python run_all_dillute.py")
  os.system("python run_all_excluded_volume.py")
  os.system("python run_all_excluded_volume_dillute.py")
  print("spatiocyte");
  os.chdir("../spatiocyte")
  #os.system("python run_all.py")
  #os.system("python run_all_dillute.py")
  #os.system("python run_all_point.py")
  os.system("python run_all_point_dillute.py")

