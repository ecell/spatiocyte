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
  row.resize(4);
  col.resize(4);
  lay.resize(4);
  corn.resize(4);
  for(unsigned i(0);i<4;i++)
    {
      theSpatiocyteStepper->point2global(aQuad[i], row[i], col[i],
                                         lay[i]);
      corn[i]=theSpatiocyteStepper->global2coord(row[i], col[i], lay[i]);
    }

  for(unsigned i(0);i<corn.size();i++)
    {
      theVacantCompSpecies[0]->addMolecule(&(*theLattice)[corn[i]]);
    }
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
  theVacantCompSpecies[0]->setIsPopulated();
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
  for (int i(0);i<corn.size();i++)
    {
      theVacantCompSpecies[0]->clearMolecules();
    }
  while (theVacantSpecies->compVoxelSize() != 0)
    {
      theVacantSpecies->clearCompVoxels();
    }
  constructComp();
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
  Comp* theComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  lengthX = theComp->lengthX;
  lengthY = theComp->lengthY;
  lengthZ = theComp->lengthZ;
}

}


  //To obtain the min and max hvec
  /*std::cout<<"Input for minHvec(Use only for First mesh file): "<<std::endl;
  std::cout<<"minHvecX: "<<minhvecX<<std::endl;
  std::cout<<"minHvecY: "<<minhvecY<<std::endl;
  std::cout<<"minHvecZ: "<<minhvecZ<<std::endl;
  std::cout<<"Input for maxHvec(Use only for Last mesh file): "<<std::endl;
  std::cout<<"maxHvecX: "<<maxhvecX<<std::endl;
  std::cout<<"maxHvecY: "<<maxhvecY<<std::endl;
  std::cout<<"maxHvecZ: "<<maxhvecZ<<std::endl;*/
  //CheckSize
  /*std::cout << "vacant size:" << theVacantSpecies->size() << std::endl;
  std::cout << "first size:" << surfaceCoords.size() << std::endl;
  setPopulated(); 
  std::cout << "vacant size:" << theVacantSpecies->size() << std::endl;
  std::cout << "vacant comp size:" << theVacantCompSpecies[0]->size() 
  << std::endl;*/
/*void MigrationProcess::assignQuad()
{
  quadIndex.resize(nq);
  for (unsigned a=0;a<nq;a++)
    {
      quadIndex[a].resize(4);
    }
  for(int i(0);i<quadIndex.size();i++)	
    {	
      for (int k(0);k<4;k++)
        {
	 	 	    quadIndex[i][k]=isoq[i][k];
				}	 
    }	
}	

void MigrationProcess::assignEdge()
{
	edgeIndex.resize(nl);
	for (unsigned a(0);a<nl;a++)
		{
			edgeIndex[a].resize(2);
		}
	for(int i(0);i<edgeIndex.size();i++)
		{
			for (int k(0);k<2;k++)
				{
					edgeIndex[i][k]=isol[i][k];
				}
		}
}*/

/*void MigrationProcess::assignNeigh()
{
  neigh.resize(nq);
  for (unsigned a=0;a<nq;a++)
    {
      neigh[a].reserve(15);
    }
  for(int i(0);i<nq;i++)	
    {		
      for (int j(0);j<4;j++)
        {
	        for (int k(0);k<nq;k++)	
            {
	            for (int l(0);l<4;l++)
 	              {
		              int nodeID(isoq[k][l]);
		              if(nodeID==isoq[i][j])
		                {
		                  std::vector<int>::iterator it(std::find(neigh[i].
                                                    begin(),neigh[i].end(), k));	
		                  if(k!=i && it==neigh[i].end())
		                    {	
			                    neigh[i].push_back(k);  
			                  }
		                }
		            }
	          }
	      }
    }
}*/
  /*for (int k(0);k<2;k++)//no of nodes on 1 edge   
    {
      if (k==0)//stack 0 
        {
          aQuad[3].x=((hvec[isol[i][k]-1][j][0]+(-minhvecX))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[3].y=((hvec[isol[i][k]-1][j][1]+(-minhvecY))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[3].z=((hvec[isol[i][k]-1][j][2]+(-minhvecZ))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[0].x=((hvec[isol[i][k]-1][j+1][0]+(-minhvecX))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[0].y=((hvec[isol[i][k]-1][j+1][1]+(-minhvecY))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[0].z=((hvec[isol[i][k]-1][j+1][2]+(-minhvecZ))*scalingFactor+translate)/(2*voxelRadius);
        }
      else //stack 1
        {
          aQuad[2].x=((hvec[isol[i][k]-1][j][0]+(-minhvecX))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[2].y=((hvec[isol[i][k]-1][j][1]+(-minhvecY))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[2].z=((hvec[isol[i][k]-1][j][2]+(-minhvecZ))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[1].x=((hvec[isol[i][k]-1][j+1][0]+(-minhvecX))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[1].y=((hvec[isol[i][k]-1][j+1][1]+(-minhvecY))*scalingFactor+translate)/(2*voxelRadius);
          aQuad[1].z=((hvec[isol[i][k]-1][j+1][2]+(-minhvecZ))*scalingFactor+translate)/(2*voxelRadius);
        }
    }*/
/*      aQuad[k].x=((hvec[isoq[i][k]-1][j][0]+(-minhvecX))*scalingFactor+translate)/(2*voxelRadius);
      aQuad[k].y=((hvec[isoq[i][k]-1][j][1]+(-minhvecY))*scalingFactor+translate)/(2*voxelRadius);
      aQuad[k].z=((hvec[isoq[i][k]-1][j][2]+(-minhvecZ))*scalingFactor+translate)/(2*voxelRadius);
      std::cout<<"quad "<<aQuad[k].x<<std::endl;
      std::cout<<"quad "<<aQuad[k].y<<std::endl;
      std::cout<<"quad "<<aQuad[k].z<<std::endl;*/
/*
  for (unsigned count(0);count<surfaceCoords.size();count++)
    {
      theVacantSpecies->addCompVoxel(surfaceCoords[count]);
    }
*/
  /*realHvec.resize(ns);
  for (int i(0);i<ns;i++)
    {
      realHvec[i].resize(3);
    Point bottomLeft;
    Point topRight;
      for (int j(0);j<3;j++)
        {
          realHvec[i][j].resize(3);
        }
    }
  for (int i(0);i<ns;i++)
    {
      for (int j(0);j<3;j++)
        {
	        realHvec[i][j][0] = hvec[i][j][0]+(-minhvecX);
	        realHvec[i][j][0]=(realHvec[i][j][0]*scalingFactor+translate)/
                            (2*voxelRadius);
      	  realHvec[i][j][1] = hvec[i][j][1]+(-minhvecY);
      	  realHvec[i][j][1]=(realHvec[i][j][1]*scalingFactor+translate)/ 
                            (2*voxelRadius);
      	  realHvec[i][j][2] = hvec[i][j][2]+(-minhvecZ);
      	  realHvec[i][j][2]=(realHvec[i][j][2]*scalingFactor+translate)/
                            (2*voxelRadius);
	      }
    }*/
