//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//

#include <libecs/MigrationProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(MigrationProcess, Process); 

void MigrationProcess::initializeThird()
{
  ipf = new char[20];
  for(int i(0);i<FileName.size();i++)
    {
      ipf[i] = FileName[i];
    }
  ipf[FileName.size()] = '\0';
  id1 = 11;
  id2 = 12;
  id3 = 1e-2;
  id4 = 1;
  cmdt = new double[12];
  lag = new char[10];
  lag = "lagrangian";
  nw = new char [2];
  nw = "nw"; 
  getCompartmentLength();

  std::ifstream checkFile("migration.dat");
  if(checkFile)
    {
      remove("migration.dat");
    }

  idot1=strchr(ipf,'.'); 
  idot=(idot1-ipf)+1;
  if(idot<=0 || idot>20)
  std::cout<<"invalid inputfile name"<<std::endl;
  extractnumeric(ipf,idot);
  openfile(ipf,id1,idot);
  delt=dscp[idfrm-1][2];
  timedump(delt);
  setScalingFactor();
  setCenterPoint();
  initForces();
  updateComp();
}

void MigrationProcess::initializeFifth()
{
  theInterval = tstp;
  theTime= time_; 
  thePriorityQueue->move(theQueueID);
}

void MigrationProcess::fire()
{
  for (int i(0);i<5;i++)
    {
      theSpecies[i+2]->clearMolecules();
    }
  updateComp();
  theInterval = tstp;
  theTime= time_; 
  thePriorityQueue->moveTop();
}

void MigrationProcess::setScalingFactor()
{
  double mechanocyteLengthX(maxhvecX-minhvecX);
  double mechanocyteLengthY(maxhvecY-minhvecY);
  double mechanocyteLengthZ(maxhvecZ-minhvecZ);
  double lengthXunnormalized(lengthX*(2*voxelRadius));
  double lengthYunnormalized(lengthY*(2*voxelRadius));
  double lengthZunnormalized(lengthZ*(2*voxelRadius));
  double scale(8*voxelRadius);
  double ratioX((lengthXunnormalized-scale)/mechanocyteLengthX);
  double ratioY((lengthYunnormalized-scale)/mechanocyteLengthY);
  double ratioZ((lengthZunnormalized-scale)/mechanocyteLengthZ);
  scalingFactor = std::min(ratioX,(std::min(ratioY,ratioZ)));
  translate = 4*voxelRadius;
}


void MigrationProcess::constructComp()
{
  for (int i(0);i<nq;i++)
    {
      for (int j(0);j<3;j+=2)//level
        {
          std::vector<Point> aQuad(getQuad(i,j)); 
          Point bottomLeft;
          Point topRight;
          getBox(aQuad,bottomLeft,topRight);
          setQuadVoxels(aQuad,bottomLeft,topRight);
        }
    }
  for (int i(0);i<nl;i++)
    {
      for (int j(0);j<2;j++)//level
        {
          std::vector<Point> aQuad(getEdgeQuad(i,j));
          Point bottomLeft;
          Point topRight;
          getBox(aQuad,bottomLeft,topRight);
          setQuadVoxels(aQuad,bottomLeft,topRight);
        }
    }
}

std::vector<Point> MigrationProcess::getQuad(int i, int j)
{
  std::vector<Point>aQuad;
  for (int k(0);k<4;k++)//no. of nodes on 1 aQuad
    {
      aQuad.push_back(Point(((hvec[isoq[i][k]-1][j][0]+(-minhvecX))*
                             scalingFactor+translate)/(2*voxelRadius),
                            ((hvec[isoq[i][k]-1][j][1]+(-minhvecY))*
                             scalingFactor+translate)/(2*voxelRadius),
                            ((hvec[isoq[i][k]-1][j][2]+(-minhvecZ))*
                             scalingFactor+translate)/(2*voxelRadius)));
    } 
  return aQuad;
}

std::vector<Point> MigrationProcess::getEdgeQuad(int i, int j)
{
  std::vector<Point>aQuad;
  aQuad.push_back(Point(((hvec[isol[i][0]-1][j+1][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][0]-1][j+1][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][0]-1][j+1][2]+(-minhvecZ))*
                         scalingFactor+translate)/(2*voxelRadius)));
  aQuad.push_back(Point(((hvec[isol[i][1]-1][j+1][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j+1][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j+1][2]+(-minhvecZ))*
                         scalingFactor+translate)/(2*voxelRadius)));
  aQuad.push_back(Point(((hvec[isol[i][1]-1][j][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j][2]+(-minhvecZ))*
                         scalingFactor+translate)/(2*voxelRadius)));
  aQuad.push_back(Point(((hvec[isol[i][0]-1][j][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][0]-1][j][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][0]-1][j][2]+(-minhvecZ))*
                         scalingFactor+translate)/(2*voxelRadius)));
    return aQuad;
}

void MigrationProcess::getBox(std::vector<Point>& aQuad,Point& bottomLeft, 
                              Point& topRight)
{
  double minimumX=aQuad[0].x;
  double minimumY=aQuad[0].y;
  double minimumZ=aQuad[0].z;
  double maximumX=aQuad[0].x;
  double maximumY=aQuad[0].y;
  double maximumZ=aQuad[0].z;

  for(int a=1;a<4;a++)
    {
      minimumX=std::min(aQuad[a].x,minimumX);
      minimumY=std::min(aQuad[a].y,minimumY);
      minimumZ=std::min(aQuad[a].z,minimumZ);
      maximumX=std::max(aQuad[a].x,maximumX);
      maximumY=std::max(aQuad[a].y,maximumY);
      maximumZ=std::max(aQuad[a].z,maximumZ);   	
    }

  bottomLeft.x = minimumX-(normVoxelRadius+theSpatiocyteStepper
                           ->getColLength());
  bottomLeft.y = minimumY-(normVoxelRadius+theSpatiocyteStepper
                           ->getLayerLength());
  bottomLeft.z = minimumZ-(normVoxelRadius+theSpatiocyteStepper
                           ->getRowLength());
  topRight.x = maximumX+(normVoxelRadius+theSpatiocyteStepper
                         ->getColLength());
  topRight.y = maximumY+(normVoxelRadius+theSpatiocyteStepper 
                         ->getLayerLength());
  topRight.z = maximumZ+(normVoxelRadius+theSpatiocyteStepper
                         ->getRowLength());
  topRight.x = std::min(topRight.x, lengthX);
  topRight.y = std::min(topRight.y, lengthY);
  topRight.z = std::min(topRight.z, lengthZ);
  bottomLeft.x = std::max(bottomLeft.x, 0e0);
  bottomLeft.y = std::max(bottomLeft.y, 0e0);
  bottomLeft.z = std::max(bottomLeft.z, 0e0);
}

void MigrationProcess::setQuadVoxels(std::vector<Point>& aQuad,
                                     Point& bottomLeft, Point& topRight)
{
  Point fixsurfaceNormal;
  double fixsurfaceDisplace(0);
  unsigned blRow(0);
  unsigned blLayer(0);
  unsigned blCol(0);
  theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
  unsigned trRow(0);
  unsigned trLayer(0);
  unsigned trCol(0);
  theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);

  for(unsigned d(blRow); d <= trRow; ++d)
    {
      for(unsigned e(blLayer); e <= trLayer; ++e)
        {
          for(unsigned f(blCol); f <= trCol; ++f)
            {
              unsigned m(theSpatiocyteStepper->global2coord(d, e, f));
              Point n(theSpatiocyteStepper->coord2point(m));
              calculateSurfaceNormal(aQuad[0],aQuad[1],aQuad
                                     [2],fixsurfaceNormal, fixsurfaceDisplace);
              if(isOnAboveSurface(n,fixsurfaceNormal,fixsurfaceDisplace) 
                 && isOnBelowSideSurface(aQuad[0],aQuad[1],n,fixsurfaceNormal)
                 && isOnBelowSideSurface(aQuad[1],aQuad[2],n,fixsurfaceNormal)
                 && isOnBelowSideSurface(aQuad[2],aQuad[3],n,fixsurfaceNormal)
                 && isOnBelowSideSurface(aQuad[3],aQuad[0],n,fixsurfaceNormal))
                { 	
              	  Voxel& aVoxel((*theLattice)[m]);	
                  for (unsigned i(0); i!=theAdjoiningCoordSize; i++)
                    {
                      unsigned coord(aVoxel.adjoiningCoords[i]);           
                      Point o(theSpatiocyteStepper->coord2point(coord));
                      if(isOnAboveSurface(o,fixsurfaceNormal,
                         fixsurfaceDisplace)==false)
                        {
                          if(getID((*theLattice)[m]) != theVacantSpecies
                             ->getID())
                            {  
                              theVacantSpecies->addCompVoxel(m);
                            }
                          break;
                        }
                    }
                }
            }
        }
    }    
}

void MigrationProcess::populateMolecules()
{
  for(int i(0); i !=theVacantCompSpecies.size(); ++i)
    {
      Species* aSpecies(theVacantCompSpecies[i]);
      for(int j(0); j!=aSpecies->size(); ++j)
        {
          Voxel* aVoxel(aSpecies->getMolecule(j));
          if(getID(aVoxel)== theVacantSpecies->getID())
            {
              //reassign molecule to aVoxel
              aVoxel->idx = j+theStride*aSpecies->getID();
            }
        }
    }

 for (int i(0);i<theVacantCompSpecies.size();++i)
    {
      Species* aSpecies(theVacantCompSpecies[i]);
      for (int j(0);j<aSpecies->size();++j)
        {
          Voxel* aVoxel(aSpecies->getMolecule(j));
          if(getID(aVoxel) != aSpecies->getID())
            {
              advectSurfaceMolecule(aSpecies, j);
            }
        }
    }
}


void MigrationProcess::advectSurfaceMolecule(Species* aSpecies,
                                             unsigned lastIndex)
{     
  std::vector<unsigned> stageList;
  std::vector<unsigned> stageIndex;
  stageIndex.push_back(0);
  stageList.push_back(aSpecies->getCoord(lastIndex));
  stageIndex.push_back(stageList.size());
  while (1)
    {
      for (int j(stageIndex[stageIndex.size()-2]); j != stageIndex.back(); ++j)
        {
          Voxel* currentMol(&(*theLattice)[stageList[j]]);
          for (int k(0);k<theAdjoiningCoordSize;++k)
            {
              const unsigned adjCoord(currentMol->adjoiningCoords[k]);
              Voxel* adj(&(*theLattice)[adjCoord]);
              //check if adj is vacant
              if(getSpecies(adj) == theVacantSpecies)
                {
                  replaceMolecules(stageList, stageIndex, adj, currentMol,
                                   aSpecies, lastIndex);
                  return;
                }
              //check if adj is on surface
              if(getSpecies(adj)->getVacantSpecies() == 
                 theVacantSpecies)
                {
                  if(std::find(stageList.begin(),stageList.end(),
                                    adjCoord)==stageList.end())
                    {
                      stageList.push_back(adjCoord);
                    }
                }
            }
        }
      if (stageList.size()==1)
        {
          stageList[0] = getSurfaceAdjCoord(stageList[0]);
        }
      else
        {
          stageIndex.push_back(stageList.size());
        }
    }
}

unsigned MigrationProcess::getSurfaceAdjCoord(unsigned aCoord)
{
  Voxel* aVoxel(&(*theLattice)[aCoord]);
  for (int i(0);i<theAdjoiningCoordSize;++i)
    {
      const unsigned adjCoord(aVoxel->adjoiningCoords[i]);
      Voxel* adj(&(*theLattice)[adjCoord]);
      for (int j(0);j<theAdjoiningCoordSize;++j)
        {
          const unsigned adjadjCoord(adj->adjoiningCoords[j]);
          Voxel* adjadj(&(*theLattice)[adjadjCoord]);                                     if (getSpecies(adjadj)==theVacantSpecies || 
              getSpecies(adjadj)->getVacantSpecies()==theVacantSpecies)
            {
              return adjadjCoord;
            }                    
        }
    }
  std::cout<<"Error: Couldnt get surface adjacent coord."<<std::endl;
}

void MigrationProcess::replaceMolecules(std::vector<unsigned> stageList, 
                                        std::vector<unsigned> stageIndex, 
                                        Voxel* adj, Voxel* currentMol, 
                                        Species* currSpecies, 
                                        unsigned lastIndex)
{
  Species* aSpecies(getSpecies(currentMol));
  if(stageIndex.size() == 2)
    {
      theSpecies[2]->softAddMolecule(currSpecies->getMolecule(lastIndex));
      currSpecies->softReplaceMolecule(lastIndex, adj);
      theSpecies[3]->softAddMolecule(adj);
      return;
    }
  unsigned index(aSpecies->getIndex(currentMol));
  aSpecies->softReplaceMolecule(index, adj);
  adj = currentMol;
  for (int i(stageIndex.size()-1);i>1;i--)
    {
      bool isAdded(false);
      for (int j(stageIndex[i-2]);
           j != stageIndex[i-1]; j++)
        {
          unsigned coord(stageList[j]);
          for (int k(0);k<theAdjoiningCoordSize;++k)
            {
              if (adj->adjoiningCoords[k] == coord)
                {
                  currentMol = &(*theLattice)[coord];
                  aSpecies = getSpecies(currentMol);
                  if(i == 2)
                    {
                      theSpecies[4]->softAddMolecule(currSpecies->getMolecule
                                                     (lastIndex));
                      currSpecies->softReplaceMolecule(lastIndex, adj);
                      theSpecies[5]->softAddMolecule(adj);
                      return;
                    }
                  index = aSpecies->getIndex(currentMol);
                  aSpecies->softReplaceMolecule(index, adj);
                  adj = currentMol;
                  isAdded = true;
                  break;
                }

            }
          if(isAdded)
            {
              break;
            }
        } 
    }
}



void MigrationProcess::calculateSurfaceNormal(Point& node1, Point& node2, 
                                              Point& node3,
                                              Point& fixsurfaceNormal,
                                              double& fixsurfaceDisplace)
{
  Point AB;
  Point AC;
  AB = sub(node2,node1);
  AC = sub(node3,node1);
  fixsurfaceNormal = cross(AB,AC);
  fixsurfaceNormal = norm(fixsurfaceNormal);
  fixsurfaceDisplace = dot(fixsurfaceNormal, node1);		
}	

bool MigrationProcess::isOnBelowSideSurface(Point& node1, Point& node2,  
Point& aPoint, Point& fixsurfaceNormal)
{
  Point AB;
  Point sideSurfaceNormal;
  double sideSurfaceDisplace(0);
  AB = sub(node2,node1);
  sideSurfaceNormal = cross(AB, fixsurfaceNormal);
  sideSurfaceNormal = norm(sideSurfaceNormal);
  sideSurfaceDisplace = dot(sideSurfaceNormal, node1);
  double disp(point2planeDisp(aPoint, sideSurfaceNormal, 
              sideSurfaceDisplace));
  if(disp-normVoxelRadius <= 0)
    {
      return true;
    }
  return false;
}

bool MigrationProcess::isOnAboveSurface(Point& aPoint, Point& 
fixsurfaceNormal, double& fixsurfaceDisplace)
{
  double disp(point2planeDisp(aPoint, fixsurfaceNormal, 
              fixsurfaceDisplace));
  if(disp >= 0)
    {
      return true;
    }
  return false;
}

void MigrationProcess::setPopulated()
{
  theVacantSpecies->setIsPopulated();
}


void MigrationProcess::initValue()
{
  for (int iq(0);iq<nq;iq++)
    {
      for (int i(0);i<3;i++)
        {
          vfixv[iq][i]=1e0;
        }
    }
  if (time_<=1e-5)
    {
      for (int is(0);is<ns;is++)
        {
          svec[is][0][3]=std::max(svec[is][0][3],1e-4);
        } 
      for (int il(0);il<nl;il++)
        {
          cxprm[il][3]=1e0;
          cxval[il][3]=2e-1;
          int iq = iqol[il];
          double yq = 0e0;
          double xq = 0e0;
          for (int isn(0);isn<4;isn++)
            {
              int is;
              is = isoq[iq-1][isn];
              xq = xq + hvec[is-1][0][0];
              yq = yq + hvec[is-1][0][1];
            }
          double al;
          al = atan(std::abs(yq)/(std::abs(xq)+1e-20))*180.0/3.141592;
          if (xq<0e0 || al>45e0)
            {
              cxprm[il][3]=0e0;
            }
          evec[il][3]=cxprm[il][3];
        }
    }
  else
    {
      for (int il(0);il<nl;il++)
        {
          cxprm[il][3] = evec[il][3];
          cxval[il][3] = 2e-1;
        }
    }
}

void MigrationProcess::updateComp()
{
  initVmaxCnwmin();
  clchm();
  vvec1();
  bc();
  cmstpsiz(cmdt,id3);
  initAvdtTstp();
  avgridmo(lag);
  avfield(id4,nw); 
  clphi();
  clvis();
  clpsi();
  vvec1();
  clsfr();
  clgam(aratio);      
  while(1)
    {
      modriver(icyc,epsl,idebug);
      if(icyc<10)
        {
          break;
        }
      openunit12();
    }
  writeNewFile();
  theVacantSpecies->clearCompVoxels();
  constructComp();
  populateMolecules();
  setPopulated();
}

void MigrationProcess::initForces()
{
  idphi=1; 
  idvis=1; 
  idvfx=1; 
  idgam=1; 
  idpsi=1; 
  idsfr=1; 
  idbfr=0; 
  iddrg=0; 
  idtrc=0; 
  idhyc=0; 

  clphi();
  clvis();
  clpsi();
  for (int i(0);i<2048;i++)
    {
      for (int j(0);j<3;j++)
        {
          cnode[i][j]=1e0;
        }
    }
  govolint(cnode,volint);
  gosurfintn(cnode,surfintv,surfintd,surfinte);  
  tsurf=surfintv+surfintd+surfinte;
  double const1(volint*3.0);
  double const2(4.0*3.141592);
  double const3(1./3.);  
  r0=pow((const1/const2),const3);
  a0=4.0*3.141592*(pow(r0,2));
  aratio=tsurf/a0;
  clgam(aratio);
  initValue();
  vvec1();
  clsfr();
  epsl = 1e-4;
  idebug = 1;
  int kount = 0;
  icyc = 0;
  do
    {
      modriver(icyc,epsl,idebug);
      kount = kount +1;
      openunit12();
    } while (icyc>10);
  idebug = 1;
}

void MigrationProcess::initVmaxCnwmin()
{
  double vmax;
  vmax = sqrt(pow(hvec[0][0][6],2)+pow(hvec[0][0][7],2)+
              pow(hvec[0][0][8],2));
  for (int i(0);i<ns;i++)
    {
      for (int j(0);j<3;j++)
        {
          double loopformax;
          loopformax = sqrt(pow(hvec[i][j][6],2)+pow(hvec[i][j][7],2)+
                            pow(hvec[i][j][8],2)); 
          vmax = std::max(loopformax,vmax);
        }
    }
  double cnwmin;
  cnwmin = svec[0][0][0];
  for (int i(0);i<ns;i++)
    {
      for (int j(0);j<3;j++)
        {
          double loopformin;
          loopformin = svec[i][j][0];
          cnwmin = std::min(loopformin,cnwmin);
        }
    }
  std::cout<<std::setprecision(17)<<"vmax="<<vmax<<"\t"<<"cnwmin="<<cnwmin
    <<std::endl;
  avtstep(avdt);
  if (avdt<1e-9)
    {
      std::cout<<"small avdt!"<<avdt<<std::endl;
      exit(0) ;
    }
}

void MigrationProcess::initAvdtTstp()
{
  std::cout<<std::setprecision(17)<<"chemistry step,advection step"<<"\t"
    <<tstp<<"\t"<<avdt<<std::endl;
  tstp = avdt;
  if (tstp>=tnext-time_)
      {
        tstp = tnext- time_;
        isve = 1;
      }
  else
      {
        isve = 0;
      }
  std::cout<<"advt,tstp"<<"\t"<<avdt<<"\t"<<tstp<<std::endl;
  dfdriver(lag);
  double minval;
  minval = svec[0][0][3];
  for (int i(0);i<ns;i++)
    {
      double loopforminval = svec[i][0][3];
      minval = std::min(loopforminval,minval);
    }
  std::cout<<"min 4"<<"\t"<<minval<<std::endl;
}

void MigrationProcess::writeNewFile()
{
  double alphamin(vbig);
  double ratiomn(vbig);
  double rmean(0e0);
  double rmax(0e0);
  double rmin(vbig);
  double dxy;
  double dz;
  double rxy;
  double alpha;
  double hmax;
  for (int il(0);il<nl;il++)
    {
      int is;
      is = isol[il][0];
      dxy = sqrt(pow((hvec[is-1][1][0]-hvec[is-1][0][0]),2)+
                 pow((hvec[is-1][1][1]-hvec[is-1][0][1]),2));
      dz = hvec[is-1][1][2]-hvec[is-1][0][2];
      rxy = sqrt(pow(hvec[is-1][0][0],2)+pow(hvec[is-1][0][1],2));
      rmean = rmean + rxy;
      rmax = std::max(rmax,rxy);
      rmin = std::min(rmin,rxy);
      alpha = atan(dz/dxy);
      alphamin = std::min(alpha,alphamin);
      ratiomn = std::min((hvec[is-1][1][2]/hvec[is-1][2][2]),ratiomn);
    }
  rmean = rmean/nl;
  double const4(180/3.141592);
  std::cout<<std::setprecision(17)<<"alphamin"<<"\t"<<alphamin*const4
    <<"\t"<<"ratiomn"<<"\t"<<ratiomn<<std::endl;
  time_=time_+tstp;    
  openunit12();
  double maxval;
  maxval = hvec[0][2][2];
  for(int count(0);count<ns;count++)
    {
      maxval = std::max(hvec[count][2][2],maxval);
    }
  hmax = maxval;
  govolint(cnode,volint);
  gosurfintn(cnode,surfintv,surfintd,surfinte);  
  tsurf=surfintv+surfintd+surfinte;
  double const5(volint*3.0);
  double const6(4.0*3.141592);
  double const7(1./3.);  
  r0=pow((const5/const6),const7);
  a0=4.0*3.141592*(pow(r0,2));
  aratio=tsurf/a0;
  std::cout<<"volume, aratio, time"<<"\t"<<" "<<std::setprecision(17)
    <<volint<<" "<<aratio<<" "<<time_<<std::endl;
  wrfile(isve,ipf,idot,rmin,rmax,hmax,volint,aratio);
}

void MigrationProcess::getCompartmentLength()
{
  voxelRadius = theSpatiocyteStepper->getVoxelRadius();
  normVoxelRadius = theSpatiocyteStepper->getNormalizedVoxelRadius();  
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  *theComp = *aComp;
  lengthX = aComp->lengthX;
  lengthY = aComp->lengthY;
  lengthZ = aComp->lengthZ;
}

void MigrationProcess::setCenterPoint()
{
  initminX=hvec[0][0][0];
  initminY=hvec[0][0][1];
  initminZ=hvec[0][0][2];
  initmaxX=hvec[0][0][0];
  initmaxY=hvec[0][0][1];
  initmaxZ=hvec[0][0][2];
  
  for (int i(0);i<ns;i++)
    {
      for (int j(0);j<3;j++)
        {
          initminX=std::min(hvec[i][j][0],initminX);
          initminY=std::min(hvec[i][j][1],initminY);
          initminZ=std::min(hvec[i][j][2],initminZ);
          initmaxX=std::max(hvec[i][j][0],initmaxX);
          initmaxY=std::max(hvec[i][j][1],initmaxY);
          initmaxZ=std::max(hvec[i][j][2],initmaxZ);    
        }
    }

  theComp->lengthX = ((initmaxX-initminX)*scalingFactor)/(2*voxelRadius);
  theComp->lengthY = ((initmaxY-initminY)*scalingFactor)/(2*voxelRadius);
  theComp->lengthZ = ((initmaxZ-initminZ)*scalingFactor)/(2*voxelRadius);

  theComp->centerPoint.x = ((initmaxX-initminX)*scalingFactor/2+translate)/(2*voxelRadius);
  theComp->centerPoint.y = ((initmaxY-initminY)*scalingFactor/2+translate)/(2*voxelRadius);
  theComp->centerPoint.z = ((initmaxZ-initminZ)*scalingFactor/2+translate)/(2*voxelRadius);
}

}

