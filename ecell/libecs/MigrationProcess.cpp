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

void MigrationProcess::prepreinitialize()
{
  SpatiocyteProcess::prepreinitialize();
  theSystem = getSuperSystem();
  theVacantVariable = createVariable("VACANT", theSystem);
  createVariable("PROCESS_COMPARTMENT", theSystem);
  theAddedVariable = createVariable("AddedSurface", theSystem);
  theOverlapVariable = createVariable("OverlapSurface", theSystem);
}

void MigrationProcess::preinitialize()
{
  SpatiocyteProcess::preinitialize();
}

void MigrationProcess::initialize()
{
  if(isInitialized)
    {
      return;
    }
  SpatiocyteProcess::initialize();
  theVacantSpecies = theSpatiocyteStepper->addSpecies(theVacantVariable);
  theAddedSpecies = theSpatiocyteStepper->addSpecies(theAddedVariable);
  theOverlapSpecies = theSpatiocyteStepper->addSpecies(theOverlapVariable);
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      Species* aSpecies(theSpatiocyteStepper->variable2species(
                             (*i).getVariable())); 
      if(aSpecies == NULL)
        {
          theVacantCompVariables.push_back((*i).getVariable());
        }
      else if(!(*i).getCoefficient())
        {
          theVacantCompSpecies.push_back(aSpecies);
        }
    }
  isPriorityQueued = true;
}	

void MigrationProcess::initializeFirst()
{
  SpatiocyteProcess::initializeFirst();
  theComp = theSpatiocyteStepper->system2Comp(theSystem);
  theSuperComp = theSpatiocyteStepper->system2Comp(
      theSystem->getSuperSystem());
  theVacantSpecies->setVacantSpecies(theSuperComp->vacantSpecies);
  theVacantSpecies->setComp(theComp);
  theAddedSpecies->setComp(theComp);
  theOverlapSpecies->setComp(theComp);
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      //to be overwritten by DiffusionProcess in initializeSecond:
      theVacantCompSpecies[i]->setVacantSpecies(theVacantSpecies);
      theVacantCompSpecies[i]->setComp(theComp);
    }
  theAddedSpecies->setVacantSpecies(theVacantSpecies);
  theOverlapSpecies->setVacantSpecies(theVacantSpecies);
}


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
  for (int i(0);i<2048;i++)
    {
      for (int j(0);j<3;j++)
        {
          nodeArea[2048][3] = 0;
        }
    }
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
  getCompartmentLength();
  setScalingFactor();
  setCenterPoint();
  initForces();
  initUpdateComp();
}

void MigrationProcess::initializeFifth()
{
  theInterval = tstp;
  theTime= time_; 
  thePriorityQueue->move(theQueueID);
}

void MigrationProcess::fire()
{
  std::cout<<"time  tstp: "<<time_<<"  "<<tstp<<std::endl;
  updateComp();
  theInterval = tstp*200;
  theTime += theInterval;
  //theTime = libecs::INF; 
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
  translate = 10*voxelRadius;
}


void MigrationProcess::constructComp(bool flagTmpSurface)
{
  totalVacantCompSpeciesCount.clear();
  totalVacantSpeciesCount.clear();
  vacantCompSpeciesCount.clear();
  vacantSpeciesCount.clear();
  totalVacantCompSpeciesCount.resize(theVacantCompSpecies.size());
  int quadCount;
  createList();
  setArea();
  for (int i(0);i<quadList.size();i++)
    {
      quadList[i].voxel.clear();
      quadList[i].voxel.resize(4);
    }
  vacantSpeciesCount.resize(quadList.size());
  for (int i(0);i<vacantSpeciesCount.size();i++)
    {
      vacantSpeciesCount[i].resize(4);      
    }
  vacantCompSpeciesCount.resize(quadList.size());
  for (unsigned i(0);i<vacantCompSpeciesCount.size();i++)
    {
      vacantCompSpeciesCount[i].resize(4);
      for (unsigned j(0);j<vacantCompSpeciesCount[i].size();j++)
        {
          vacantCompSpeciesCount[i][j].resize(theVacantCompSpecies.size());
        }
    }
  for (int i(0);i<nq;i++)
    {
      for (int j(0);j<3;j+=2)//level
        {
          if(!j)
            {
              quadCount = i;
            }
          else
            {
              quadCount = i+nq;
            }
          std::vector<Point> aQuad(getQuad(i,j)); 
          Point bottomLeft;
          Point topRight;
          getBox(aQuad,bottomLeft,topRight);
          setQuadVoxels(aQuad,bottomLeft,topRight,flagTmpSurface,quadCount);
        }
    }
  for (int i(0);i<nl;i++)
    {
      for (int j(0);j<2;j++)//level
        {
          if (!j)
            {
              quadCount = 2*nq+2*i;
            }
          else
            {
              quadCount = (2*nq+2*i)+1;
            }
          std::vector<Point> aQuad(getEdgeQuad(i,j));
          Point bottomLeft;
          Point topRight;
          getBox(aQuad,bottomLeft,topRight);
          setQuadVoxels(aQuad,bottomLeft,topRight,flagTmpSurface,quadCount);
        }
    }
}

void MigrationProcess::setFunctionsForConc(unsigned parentQuadID,unsigned m)
{
  setConc(parentQuadID,m);
  check.clear();
}

void MigrationProcess::createList()
{
  quadList.clear();
  nodeList.clear();
  quadList.resize(2*nq+2*nl);
  nodeList.resize(3*ns);
  setQuadList();
  setEdgeList();
}

void MigrationProcess::setQuadList()
{
  for (int quadCount(0);quadCount<nq;quadCount++)
    {
      for (int level(0);level<3;level+=2)
        { 
          for (int nodeCount(0);nodeCount<4;nodeCount++)
            {
              int stackID(isoq[quadCount][nodeCount]-1);
              double* aVector(hvec[stackID][level]);
              if (!level)
                {
                  Node& aNode(nodeList[stackID]);
                  aNode.point = Point(((aVector[0]+(-minhvecX))*scalingFactor+
                                       translate)/(2*voxelRadius),
                                      ((aVector[1]+(-minhvecY))*scalingFactor+
                                       translate)/(2*voxelRadius),
                                      ((aVector[2]+(-minhvecZ))*scalingFactor
                                      +translate)/(2*voxelRadius));
                  aNode.quadID.push_back(quadCount);
                  quadList[quadCount].nodeID[nodeCount] = stackID; 
                  nodeList[quadList[quadCount].nodeID[nodeCount]].area = 0;
                }
              else 
                {
                  Node& aNode(nodeList[stackID+ns]);
                  aNode.point = Point(((aVector[0]+(-minhvecX))*scalingFactor+
                                       translate)/(2*voxelRadius),
                                      ((aVector[1]+(-minhvecY))*scalingFactor+
                                       translate)/(2*voxelRadius),
                                      ((aVector[2]+(-minhvecZ))*scalingFactor
                                       +translate)/(2*voxelRadius));
                  aNode.quadID.push_back(quadCount+nq);
                  quadList[quadCount+nq].nodeID[nodeCount] = stackID+ns; 
                  nodeList[quadList[quadCount+nq].nodeID[nodeCount]].area = 0;
                }
            } 
        }
    }
}

void MigrationProcess::setEdgeList()
{
  for (int edgeCount(0);edgeCount<nl;edgeCount++)
    {
      for (int level(0);level<2;level++)
        {
          for (int nodeCount(0);nodeCount<2;nodeCount++)
            {
              int stackID(isol[edgeCount][nodeCount]-1);
              double *aVector(hvec[stackID][1]);
              int newQuadID(2*nq+2*edgeCount);
              if(!level)
                {
                  Node& aNode(nodeList[stackID+(2*ns)]);
                  aNode.point = Point(((aVector[0]+(-minhvecX))*
                                       scalingFactor+translate)/(2*voxelRadius),
                                      ((aVector[1]+(-minhvecY))*
                                       scalingFactor+translate)/(2*voxelRadius),
                                      ((aVector[2]+(-minhvecZ))*
                                       scalingFactor+translate)/(2*voxelRadius)); 
                  aNode.quadID.push_back(newQuadID);
                  quadList[newQuadID].nodeID[nodeCount] = stackID;
                  if (nodeCount == 1)
                    {
                      nodeList[stackID].quadID.push_back(newQuadID);
                      quadList[newQuadID].nodeID[nodeCount+1] = 
                        stackID+(2*ns);
                    }
                  else
                    {
                      nodeList[stackID].quadID.push_back(newQuadID);
                      quadList[newQuadID].nodeID[nodeCount+3] = 
                        stackID+(2*ns);
                    }
                }
              else
                { 
                  Node& aNode(nodeList[stackID+(2*ns)]);
                  aNode.quadID.push_back(newQuadID+1);
                  quadList[newQuadID+1].nodeID[nodeCount] = stackID+
                    (2*ns);
                  if (nodeCount == 1)
                    {
                      nodeList[stackID+ns].quadID.push_back(newQuadID+1);
                      quadList[newQuadID+1].nodeID[nodeCount+1] = 
                        stackID+ns;
                    }
                  else
                    {
                      nodeList[stackID+ns].quadID.push_back(newQuadID+1);
                      quadList[newQuadID+1].nodeID[nodeCount+3] = 
                        stackID+ns;
                    }
                }
            }
        }
    }
}

void MigrationProcess::setArea()
{
  Point cenLine[4];
  Point normVector1;
  Point normVector2;
  Point normVector3;
  Point unitNormVector1;
  //cenQuad.resize(2*nq);
  double mHor;
  double mVer;
  double cHor;
  double cVer;
  for (int quadCount(0);quadCount<quadList.size();quadCount++)
    {
      Quad& aQuad(quadList[quadCount]);
      setLineCenter(aQuad);
      Point vector1(sub(nodeList[aQuad.nodeID[1]].point,
                        nodeList[aQuad.nodeID[0]].point));
      Point vector2(sub(nodeList[aQuad.nodeID[3]].point,
                        nodeList[aQuad.nodeID[0]].point)),
      normVector1 = cross(vector1,vector2);
      Point vector3(sub(aQuad.cPoint[1],aQuad.cPoint[3]));
      normVector2 = cross(vector3,normVector1);
      Point vector4(sub(aQuad.cPoint[2],aQuad.cPoint[0]));
      normVector3 = cross(vector4,normVector1);
      Point& node1(nodeList[aQuad.nodeID[0]].point);
      Point& node2(nodeList[aQuad.nodeID[1]].point);
      Point& node3(nodeList[aQuad.nodeID[2]].point);
      Point& node4(nodeList[aQuad.nodeID[3]].point);
      unitNormVector1 = norm(normVector1);
      setCenQuad(nodeList[aQuad.nodeID[0]].point,aQuad.cPoint[3],
                 aQuad.cPoint[0],normVector1,normVector2,normVector3,
                 aQuad.cPoint[4]);
      nodeList[aQuad.nodeID[0]].area += dot(unitNormVector1,
                                           cross(sub(aQuad.cPoint[4],node1),
                                                 sub(aQuad.cPoint[3],
                                                     aQuad.cPoint[0])))/2;
      nodeList[aQuad.nodeID[1]].area += dot(unitNormVector1,
                                           cross(sub(aQuad.cPoint[1],
                                                     aQuad.cPoint[0]),
                                                 sub(aQuad.cPoint[4],
                                                     node2)))/2;
      nodeList[aQuad.nodeID[2]].area += dot(unitNormVector1,
                                           cross(sub(node3,aQuad.cPoint[4]),
                                                 sub(aQuad.cPoint[2],
                                                     aQuad.cPoint[1])))/2;
      nodeList[aQuad.nodeID[3]].area += dot(unitNormVector1,
                                           cross(sub(aQuad.cPoint[2],
                                                     aQuad.cPoint[3]),
                                                 sub(node4,
                                                     aQuad.cPoint[4])))/2;
                                                
    }

}

void MigrationProcess::setConc(unsigned parentQuadID,unsigned m)
{
  unsigned count(0);
  unsigned checkRulesCount(0);
  Point fixsurfaceNormal;
  double fixsurfaceDisplace(0);
  Quad& inQuad(quadList[parentQuadID]);
  std::vector<Point> aQuad(getGaussQuad(nodeList[inQuad.nodeID[0]].point,
                                        inQuad.cPoint[0],inQuad.cPoint[4],
                                        inQuad.cPoint[3]));
  calculateSurfaceNormal(aQuad[2],aQuad[1],aQuad[3],fixsurfaceNormal,
                         fixsurfaceDisplace);
  checkRulesCount = countMol(aQuad,count,m,checkRulesCount,parentQuadID,
                             fixsurfaceNormal,fixsurfaceDisplace);
  count++;
  aQuad = getGaussQuad(nodeList[inQuad.nodeID[1]].point,
                       inQuad.cPoint[1],inQuad.cPoint[4],inQuad.cPoint[0]);
  checkRulesCount = countMol(aQuad,count,m,checkRulesCount,parentQuadID,
                             fixsurfaceNormal,fixsurfaceDisplace);
  count++;
  aQuad = getGaussQuad(nodeList[inQuad.nodeID[2]].point,inQuad.cPoint[2],
                       inQuad.cPoint[4],inQuad.cPoint[1]);
  checkRulesCount = countMol(aQuad,count,m,checkRulesCount,parentQuadID,
                             fixsurfaceNormal,fixsurfaceDisplace);
  count++;
  aQuad = getGaussQuad(nodeList[inQuad.nodeID[3]].point,inQuad.cPoint[3],
                       inQuad.cPoint[4],inQuad.cPoint[2]);
  checkRulesCount = countMol(aQuad,count,m,checkRulesCount,parentQuadID,
                             fixsurfaceNormal,fixsurfaceDisplace);
  if (checkRulesCount != 1)
    {
      Point fixsurfaceNormal1;
      double fixsurfaceDisplace1(0);
      std::cout<<"Every voxel should only satisfy one rule"<<std::endl;
      std::cout<<"checkRulesCount: "<<parentQuadID<<" "<<m<<" "<<
      checkRulesCount<<std::endl;
      Point n(theSpatiocyteStepper->coord2point(m));
      std::cout<<"Point for the voxels that count =/= 1: "<<n.x<<" "<<n.y<<" "<<
        n.z<<std::endl;
      std::cout<<"Points at the left: "<<inQuad.cPoint[3].x<<" "<<
        inQuad.cPoint[3].y<<" "<<inQuad.cPoint[3].z<<std::endl;
      std::cout<<"Points at the right: "<<inQuad.cPoint[1].x<<" "<<
        inQuad.cPoint[1].y<<" "<<inQuad.cPoint[1].z<<std::endl;
      std::cout<<"Points at the top: "<<inQuad.cPoint[2].x<<" "<<
        inQuad.cPoint[2].y<<" "<<inQuad.cPoint[2].z<<std::endl;
      std::cout<<"Points at the bottom: "<<inQuad.cPoint[0].x<<" "<<
        inQuad.cPoint[0].y<<" "<<inQuad.cPoint[0].z<<std::endl;
      calculateSurfaceNormal(inQuad.cPoint[4],inQuad.cPoint[0],inQuad.cPoint[3],
                             fixsurfaceNormal1,fixsurfaceDisplace1);
      if(isOnBelowInternalPlane(inQuad.cPoint[0],inQuad.cPoint[2],n,
                                fixsurfaceNormal1) &&
         isAboveInternalPlane(inQuad.cPoint[1],inQuad.cPoint[3],n,
                              fixsurfaceNormal1))
        {
          std::cout<<"it should be in quad 3"<<std::endl;
        }
      else if(isAboveInternalPlane(inQuad.cPoint[0],inQuad.cPoint[2],n,
                                fixsurfaceNormal1) &&
         isAboveInternalPlane(inQuad.cPoint[1],inQuad.cPoint[3],n,
                              fixsurfaceNormal1))
        {
          std::cout<<"it should be in quad 2"<<std::endl;
        }
      else if(isAboveInternalPlane(inQuad.cPoint[0],inQuad.cPoint[2],n,
                                fixsurfaceNormal1) &&
         isOnBelowInternalPlane(inQuad.cPoint[1],inQuad.cPoint[3],n,
                              fixsurfaceNormal1))
        {
          std::cout<<"it should be in quad 1"<<std::endl;
        }
      else if(isOnBelowInternalPlane(inQuad.cPoint[0],inQuad.cPoint[2],n,
                                fixsurfaceNormal1) &&
         isOnBelowInternalPlane(inQuad.cPoint[1],inQuad.cPoint[3],n,
                              fixsurfaceNormal1))
        {
          std::cout<<"it should be in quad 0"<<std::endl;
        }
      else
        {
          std::cout<<"still not in quad?"<<std::endl;
        }
      exit(0);
    }
}

std::vector<Point> MigrationProcess::getGaussQuad(Point& A, Point& B, 
                                                  Point& C, Point& D)
{
  std::vector<Point> aQuad;
  aQuad.push_back(A);
  aQuad.push_back(B);
  aQuad.push_back(C);
  aQuad.push_back(D);
  return aQuad;
}

void MigrationProcess::setCenQuad(Point& A,Point& B, Point& P, Point& N, 
                                   Point& M, Point& O, 
                                   Point& cenQuad)
{
  cenQuad.y =  (((A.x*M.z*N.x*O.x)/(-M.z*N.x+M.x*N.z))+
               ((A.y*M.z*N.y*O.x)/(-M.z*N.x+M.x*N.z))-
               ((B.x*M.x*N.z*O.x)/(-M.z*N.x+M.x*N.z))-
               ((B.y*M.y*N.z*O.x)/(-M.z*N.x+M.x*N.z))+
               ((A.z*M.z*N.z*O.x)/(-M.z*N.x+M.x*N.z))-
               ((B.z*M.z*N.z*O.x)/(-M.z*N.x+M.x*N.z))-
               A.z*O.z-((A.x*N.x*O.z)/N.z)-((A.y*N.y*O.z)/N.z)+
               ((B.x*M.x*N.x*O.z)/(-M.z*N.x+M.x*N.z))+
               ((B.y*M.y*N.x*O.z)/(-M.z*N.x+M.x*N.z))-
               ((A.z*M.z*N.x*O.z)/(-M.z*N.x+M.x*N.z))+
               ((B.z*M.z*N.x*O.z)/(-M.z*N.x+M.x*N.z))-
               ((A.x*M.z*pow(N.x,2)*O.z)/(N.z*(-M.z*N.x+M.x*N.z)))-
               ((A.y*M.z*N.x*N.y*O.z)/(N.z*(-M.z*N.x+M.x*N.z)))+
               O.x*P.x+O.y*P.y+O.z*P.z)/
               (((M.z*N.y*O.x)/(-M.z*N.x+M.x*N.z))-
               ((M.y*N.z*O.x)/(-M.z*N.x+M.x*N.z))+O.y-
               ((N.y*O.z)/N.z)+((M.y*N.x*O.z)/(-M.z*N.x+M.x*N.z))-
               ((M.z*N.x*N.y*O.z)/(N.z*(-M.z*N.x+M.x*N.z))));
  
  cenQuad.x = (1/(-M.z*N.x+M.x*N.z))*(-A.x*M.z*N.x-A.y*M.z*N.y+B.x*M.x*N.z+
                                     B.y*M.y*N.z-A.z*M.z*N.z+B.z*M.z*N.z+M.z*
                                     N.y*cenQuad.y-M.y*N.z*cenQuad.y);
    
  cenQuad.z = (A.x*N.x+A.y*N.y+A.z*N.z-N.x*cenQuad.x-N.y*cenQuad.y)/N.z;
    
}

void MigrationProcess::setLineCenter(Quad& aQuad)
{
  for (int nodeCount(0);nodeCount<4;nodeCount++)
    {
      if (nodeCount==3)
        {
          aQuad.cPoint[nodeCount]
            = addDivide(nodeList[aQuad.nodeID[nodeCount]].point,
                        nodeList[aQuad.nodeID[nodeCount-3]].point);
        }
      else
        {
          aQuad.cPoint[nodeCount] 
            = addDivide(nodeList[aQuad.nodeID[nodeCount]].point,
                        nodeList[aQuad.nodeID[nodeCount+1]].point);
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
  aQuad.push_back(Point(((hvec[isol[i][0]-1][j][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][0]-1][j][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][0]-1][j][2]+(-minhvecZ))*
                         scalingFactor+translate)/(2*voxelRadius)));
  aQuad.push_back(Point(((hvec[isol[i][1]-1][j][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j][2]+(-minhvecZ))*
                         scalingFactor+translate)/(2*voxelRadius)));
  aQuad.push_back(Point(((hvec[isol[i][1]-1][j+1][0]+(-minhvecX))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j+1][1]+(-minhvecY))*
                         scalingFactor+translate)/(2*voxelRadius),
                        ((hvec[isol[i][1]-1][j+1][2]+(-minhvecZ))*
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

unsigned MigrationProcess::countMol(std::vector<Point>& aQuad,unsigned count,
                                unsigned m,unsigned checkRulesCount,
                                unsigned parentQuadID,Point& fixsurfaceNormal, 
                                double fixsurfaceDisplace)
{
  Point n(theSpatiocyteStepper->coord2point(m));
  if(isOnBelowInternalPlane(aQuad[1],aQuad[2],n,fixsurfaceNormal)
  && isAboveInternalPlane(aQuad[3],aQuad[2],n,fixsurfaceNormal))
    { 	
      setMolCount(m,count,parentQuadID,fixsurfaceNormal,
                  fixsurfaceDisplace);
      checkRulesCount++;
      return checkRulesCount;
    }
  return checkRulesCount;
}

void MigrationProcess::setMolCount(unsigned m,unsigned count,
                                   unsigned parentQuadID,
                                   Point& fixsurfaceNormal,
                                   double fixsurfaceDisplace)
{
  Voxel* aVoxel(&(*theLattice)[m]);
  quadList[parentQuadID].voxel[count].push_back(m);
}

void MigrationProcess::setQuadVoxels(std::vector<Point>& aQuad,
                                     Point& bottomLeft, Point& topRight,
                                     bool flagTmpSurface, unsigned parentQuadID)
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
                  setVacantSpecies(m,fixsurfaceNormal,fixsurfaceDisplace,
                                   flagTmpSurface,parentQuadID);
                }
            }
        }
    }    
}

void MigrationProcess::setVacantSpecies(unsigned m,Point& fixsurfaceNormal,
                                        double fixsurfaceDisplace,
                                        bool flagTmpSurface,
                                        unsigned parentQuadID)
{
  Voxel* aVoxel(&(*theLattice)[m]);	
  Species* aSpecies(getSpecies(aVoxel));
  for (unsigned i(0); i!=theAdjoiningCoordSize; i++)
    {
      unsigned coord(aVoxel->adjoiningCoords[i]);           
      Point o(theSpatiocyteStepper->coord2point(coord));
      if(isOnAboveSurface(o,fixsurfaceNormal,
         fixsurfaceDisplace)==false)
        {
          if (flagTmpSurface == true)
            {
              if (getID(aVoxel) != theOverlapSpecies->getID() &&
                  getID(aVoxel) != theAddedSpecies->getID())
                {
                  setFunctionsForConc(parentQuadID,m);
                  if (aSpecies->getComp() == theComp)
                    {
                      theOverlapSpecies->addMoleculeDirect(aVoxel); 
                      return;
                    }
                  else 
                    {
                    theAddedSpecies->addMoleculeDirect(aVoxel);
                    return;
                    }                    
                }
            }
          if (flagTmpSurface == false)
            {
              if(getID(aVoxel) != theVacantSpecies->getID())
                {  
                  theVacantSpecies->addCompVoxel(m);
                  return;
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
          //Unmoved aVoxel already reassigned as aSpecies.
          //check the rest of aVoxel that move after 1 tstp.
          if(getID(aVoxel) != aSpecies->getID())
            {
              advectSurfaceMolecule(aSpecies, j);
            }
        }
    }
  for (int i(0);i<nodeList.size();i++)
    {
      nodeList[i].conc.resize(theVacantCompSpecies.size()+1);
      for (int j(0);j<nodeList[i].conc.size();j++)
        {
          nodeList[i].conc[j] = 0;
        }
    }
  for (int i(0);i<quadList.size();i++)
    {
      Quad& aQuad(quadList[i]);
      for (int j(0);j<4;j++)
        {
          for (int k(0);k<aQuad.voxel[j].size();k++)
            {
              Voxel* aVoxel(&(*theLattice)[aQuad.voxel[j][k]]);
              if (getID(aVoxel) == theVacantSpecies->getID())
                {
                  nodeList[aQuad.nodeID[j]].conc
                    [theVacantCompSpecies.size()] +=1;
                  vacantSpeciesCount[i][j].push_back(aQuad.voxel[j][k]);
                  totalVacantSpeciesCount.push_back(aQuad.voxel[j][k]);
                  continue;
                } 
              for (unsigned l(0);l!=theVacantCompSpecies.size();++l)
                {
                  Species* aSpecies(theVacantCompSpecies[l]);
                  if (getID(aVoxel) == aSpecies->getID())
                    {
                      nodeList[aQuad.nodeID[j]].conc[l] +=1;
                      vacantCompSpeciesCount[i][j][l].push_back
                        (aQuad.voxel[j][k]);
                      totalVacantCompSpeciesCount[l].push_back
                        (aQuad.voxel[j][k]);
                      continue;
                    }
                }
             }
        }
    }

  unsigned checkAfterSpeciesCount(0);
  for (int i(0);i<nodeList.size();i++)
    {
      for (int j(0);j<nodeList[i].conc.size();j++)
        {
          checkAfterSpeciesCount += nodeList[i].conc[j];
        }
    }
  unsigned checkPI3KCount(0);
  for (int i(0);i<nodeList.size();i++)
    {
      checkPI3KCount += nodeList[i].conc[0];
    }
  std::cout<<"TOTAL NUMBER COUNT: "<<checkPI3KCount<<std::endl;
  for (int i(0);i<nodeList.size();i++)
    {
      nodeList[i].conc[theVacantCompSpecies.size()] /= nodeList[i].area;
      for (int m(0);m<theVacantCompSpecies.size();++m)
        {
          nodeList[i].conc[m] /= nodeList[i].area;
        }
    }
  std::cout<<"total Vacant Species Count: "<<totalVacantSpeciesCount.size()<<
    std::endl;
  std::cout<<"Original Vacant Species Count: "<<theVacantSpecies->size()<<std::endl;
  unsigned checkTotalCount(theVacantSpecies->size());
  for (int i(0);i<theVacantCompSpecies.size();i++)
    {
      std::cout<<"total VCS["<<i<<"]: "<<getIDString(theVacantCompSpecies[i])
        <<totalVacantCompSpeciesCount[i].size()<<std::endl;
      std::cout<<"Original VCS["<<i<<"]: "<<theVacantCompSpecies[i]->size()
        <<std::endl;
      checkTotalCount -= totalVacantCompSpeciesCount[i].size();
    }
  if (checkTotalCount != totalVacantSpeciesCount.size())
    {
      std::cout<<"total of VS and VCS are not equal to the original VS count: "
       <<checkTotalCount<<std::endl; 
      exit(0);
    }
  if (checkAfterSpeciesCount != theVacantSpecies->size())
    {
      std::cout<<"Error after Species count before concentration calculation."
        <<std::endl;
      std::cout<<checkAfterSpeciesCount<<std::endl;
      exit(0);
    }
  theVacantSpecies->setIsPopulated();
  for (int i(0);i<theVacantCompSpecies.size();i++)
    {
      theVacantCompSpecies[i]->setIsPopulated();
    }
}

void MigrationProcess::updateConc()
{
  for (int i(0);i<ns;i++)
    {
      for (int j(0);j<3;j+=2)
        {
          for (int k(0);k<theVacantCompSpecies.size();k++)
            {
              if(!j)
                {
                  svec[i][j][7+k]=nodeList[i].conc[k];
                }
              else
                {
                  svec[i][j][7+k]=nodeList[i+ns].conc[k];
                }
            }
        }
    }
  for (int i(0);i<nl;i++)
    {
      for (int j(0);j<2;j++)
        {
          int stackID(isol[i][j]-1);
          for (int k(0);k<theVacantCompSpecies.size();k++)
            {
              svec[i][1][7+k]=nodeList[stackID+(2*ns)].conc[k];
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
          Voxel* adjadj(&(*theLattice)[adjadjCoord]);                                     
          if (getSpecies(adjadj)==theVacantSpecies || 
              getSpecies(adjadj)->getVacantSpecies()==theVacantSpecies)
            {
              return adjadjCoord;
            }                    
        }
    }
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
      currSpecies->softReplaceMolecule(lastIndex, adj);
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
                      currSpecies->softReplaceMolecule(lastIndex, adj);
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

unsigned MigrationProcess::getLatticeResizeCoord(unsigned aStarCoord)
{
  theComp->dimension = 2;
  for (int i(0);i<theVacantCompSpecies.size();i++)
    {
      theVacantCompSpecies[i]->setDimension(theComp->dimension);
      theVacantCompSpecies[i]->setIsFixedAdjoins(false);
    }
  theVacantSpecies->setDimension(theComp->dimension);
  theVacantSpecies->setIsFixedAdjoins(false);
}

void MigrationProcess::optimizeSurfaceVoxel()
{
  Species* aSpecies(theVacantSpecies);
  for (int j(0);j!=aSpecies->size();++j)
    {
      unsigned aCoord(aSpecies->getCoord(j));
      Voxel& aVoxel((*theLattice)[aCoord]);
      unsigned* forward(aVoxel.adjoiningCoords);
      std::vector<unsigned> adjoiningCopy;
      for(unsigned k(0); k != theAdjoiningCoordSize; ++k)
        {
          adjoiningCopy.push_back(forward[k]);
        }
      //Separate adjoining surface voxels and adjoining volume voxels.
      //Put the adjoining surface voxels at the beginning of the
      //adjoiningCoords list while the volume voxels are put at the end:
      for(std::vector<unsigned>::iterator l(adjoiningCopy.begin());
          l != adjoiningCopy.end(); ++l)
        {
          unsigned anID(getID((*theLattice)[*l]));
          if(anID == theVacantSpecies->getID())
            {
              aVoxel.adjoiningCoords[l-adjoiningCopy.begin()]=*forward;
              (*forward) = (*l);
              ++forward;
            }
        } 
      aVoxel.diffuseSize = forward-aVoxel.adjoiningCoords;
    }    
}

void MigrationProcess::fireOptimizeSurfaceVoxel(Voxel& aVoxel)
{
  unsigned* forward(aVoxel.adjoiningCoords);
  std::vector<unsigned> adjoiningCopy;
  for(unsigned k(0); k != theAdjoiningCoordSize; ++k)
    {
      adjoiningCopy.push_back(forward[k]);
    }
  //Separate adjoining surface voxels and adjoining volume voxels.
  //Put the adjoining surface voxels at the beginning of the
  //adjoiningCoords list while the volume voxels are put at the end:
  for(std::vector<unsigned>::iterator l(adjoiningCopy.begin());
      l != adjoiningCopy.end(); ++l)
    {
      Species* aSpecies(getSpecies((*theLattice)[*l]));
      if (aSpecies->getComp() == theComp)
        {
          aVoxel.adjoiningCoords[l-adjoiningCopy.begin()]=*forward;
          (*forward) = (*l);
          ++forward;
        }
    } 
  aVoxel.diffuseSize = forward-aVoxel.adjoiningCoords;
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

bool MigrationProcess::isOnBelowInternalPlane(Point& node1, Point& node2,  
Point& aPoint, Point& fixsurfaceNormal)
{
  Point AB;
  Point sideSurfaceNormal;
  double sideSurfaceDisplace(0);
  AB = sub(node2,node1);
  sideSurfaceNormal = cross(AB,fixsurfaceNormal);
  sideSurfaceNormal = norm(sideSurfaceNormal);
  sideSurfaceDisplace = dot(sideSurfaceNormal, node1);
  double disp(point2planeDisp(aPoint, sideSurfaceNormal, 
              sideSurfaceDisplace));
  if(disp >= 0)
    {
      return true;
    }
  return false;
}

bool MigrationProcess::isAboveInternalPlane(Point& node1, Point& node2,  
Point& aPoint, Point& fixsurfaceNormal)
{
  Point AB;
  Point sideSurfaceNormal;
  double sideSurfaceDisplace(0);
  AB = sub(node2,node1);
  sideSurfaceNormal = cross(fixsurfaceNormal,AB);
  sideSurfaceNormal = norm(sideSurfaceNormal);
  sideSurfaceDisplace = dot(sideSurfaceNormal, node1);
  double disp(point2planeDisp(aPoint, sideSurfaceNormal, 
              sideSurfaceDisplace));
  if(disp > 0)
    {
      return true;
    }
  return false;
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

void MigrationProcess::initUpdateComp()
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
  constructComp(false);
  optimizeSurfaceVoxel();
  updateArea();
}

void MigrationProcess::updateArea()
{
  theComp->actualArea =  (72*pow(voxelRadius,2))*
    theVacantSpecies->size()/(6*pow(2,0.5)+4*pow(3,0.5)+3*pow(6, 0.5));
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
  constructComp(true);
  removeOldSurfaceVoxel();
  getNewSurfaceVoxel();
  updateArea();
  populateMolecules();
  updateConc();
}

void MigrationProcess::removeOldSurfaceVoxel()
{
  Species* aSpecies(theVacantSpecies);
  for (int i(0);i<aSpecies->size();++i)
    {
      unsigned coord(aSpecies->getCoord(i));
      Voxel& aVoxel((*theLattice)[coord]);
      if(getID(aVoxel) != theOverlapSpecies->getID())
        {
          aSpecies->removeCompVoxel(i);
          for (int j(0);j<theAdjoiningCoordSize;++j)
            {
              unsigned adjCoord(aVoxel.adjoiningCoords[j]);
              Voxel& adj((*theLattice)[adjCoord]);
              if (getID(adj) == theOverlapSpecies->getID() || 
                  getID(adj) == theAddedSpecies->getID())
                {
                  fireOptimizeSurfaceVoxel(adj);
                }
            }
          --i;
        }
    }
}

void MigrationProcess::getNewSurfaceVoxel()
{
  std::vector<unsigned> overlapList;
  for (int i(0);i<theAddedSpecies->size();++i)
    {
      Voxel& aVoxel(*theAddedSpecies->getMolecule(i));
      for (int j(0);j<theAdjoiningCoordSize;++j)
        {
          unsigned adjCoord(aVoxel.adjoiningCoords[j]);
          Voxel& adj((*theLattice)[adjCoord]);
          if (getID(adj) == theOverlapSpecies->getID())  
            {
              fireOptimizeSurfaceVoxel(adj);
            }
        }
      fireOptimizeSurfaceVoxel(aVoxel);
      theVacantSpecies->addCompVoxel(theAddedSpecies->getCoord(i));          
    }
  for (int i(0);i<theOverlapSpecies->size();++i)
    {
      Voxel* aVoxel(theOverlapSpecies->getMolecule(i));
      aVoxel->idx = theVacantSpecies->getID()*theStride;
    }
  theOverlapSpecies->clearMolecules();
  theAddedSpecies->clearMolecules();
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
  lengthX = theSuperComp->lengthX;
  lengthY = theSuperComp->lengthY;
  lengthZ = theSuperComp->lengthZ;
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

void MigrationProcess::printParameters()
{
  cout << getPropertyInterface().getClassName() << "[" <<
    getFullID().asString() << "]" << std::endl;
  cout << "  " << getIDString(theVacantSpecies) << 
    " number:" << theVacantSpecies->size() << std::endl;
  cout << "  [" << theComp->actualArea << " m^2] Actual surface area " <<
    "{S = (72*r_v^2)*n_s/(6*2^0.5+4*3^0.5+3*6^0.5)}" << std::endl;
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      cout << "    " << getIDString(theVacantCompSpecies[i]) <<
        " number:" << theVacantCompSpecies[i]->size() << std::endl;
    }
}

}

