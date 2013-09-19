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

#include <libecs/MechanicsProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(MechanicsProcess, Process); 

  void MechanicsProcess::initializeThird()
    {
   //std::cout<<"1"<<std::endl;
   ipf = new char[10];
   ipf = "fitz.000";
   id1 = 11;
   id2 = 12;
   id3 = 1*10^-2;
   cmdt = new double[12];
   eul = new char[10];
   eul = "eulerian";
   surfaceDisplace = 0;
   double VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
   double NormVoxelRadius = theSpatiocyteStepper->getNormalizedVoxelRadius();  

   //std::cout<<"2"<<std::endl;
   Comp* theComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
   //std::cout<<"3"<<std::endl;  
   //theComp = aComp;
   //std::cout<<"4"<<std::endl;
   std::cout<<"lengthX: "<<theComp->lengthX<<std::endl;
   std::cout<<"lengthY: "<<theComp->lengthY<<std::endl;
   std::cout<<"lengthZ: "<<theComp->lengthZ<<std::endl;
   std::cout<<"minRow: "<<theComp->minRow<<std::endl;
   std::cout<<"minCol: "<<theComp->minCol<<std::endl;
   std::cout<<"minLayer: "<<theComp->minLayer<<std::endl;
   std::cout<<"maxRow: "<<theComp->maxRow<<std::endl;
   std::cout<<"maxCol: "<<theComp->maxCol<<std::endl;
   std::cout<<"maxLayer: "<<theComp->maxLayer<<std::endl;
   std::cout<<"VoxelRad: "<<VoxelRadius<<std::endl;
   std::cout<<"NormVoxelRad: "<<NormVoxelRadius<<std::endl;

   std::vector<Point> node;
   std::vector<unsigned>row;
   std::vector<unsigned>col;
   std::vector<unsigned>lay;
   std::vector<unsigned>corn;

   node.resize(4);
   node[0].x = 33;
   node[0].y = 35;
   node[0].z = 34;
   node[1].x = 45;
   node[1].y = 64;
   node[1].z = 46;
   node[2].x = 57;
   node[2].y = 36;
   node[2].z = 58;
   node[3].x = 69;
   node[3].y = 68;
   node[3].z = 70;

   row.resize(4);
   col.resize(4);
   lay.resize(4);
   corn.resize(4);

   for(unsigned i(0);i<3;i++)
	{
   	theSpatiocyteStepper->point2global(node[i], row[i], col[i], lay[i]);
   	corn[i]=theSpatiocyteStepper->global2coord(row[i], col[i], lay[i]);
	}


   double minimumX=node[0].x;
   double minimumY=node[0].y;
   double minimumZ=node[0].z;
   double maximumX=node[0].x;
   double maximumY=node[0].y;
   double maximumZ=node[0].z;

   for(int a=1;a<3;a++)
  {
   minimumX=std::min(node[a].x,minimumX);
   minimumY=std::min(node[a].y,minimumY);
   minimumZ=std::min(node[a].z,minimumZ);
   maximumX=std::max(node[a].x,maximumX);
   maximumY=std::max(node[a].y,maximumY);
   maximumZ=std::max(node[a].z,maximumZ);   
  }

  /*std::cout<<"minX: "<<minimumX<<std::endl;
  std::cout<<"minY: "<<minimumY<<std::endl;
  std::cout<<"minZ: "<<minimumZ<<std::endl;
  std::cout<<"minX: "<<maximumX<<std::endl;
  std::cout<<"minY: "<<maximumY<<std::endl;
  std::cout<<"minZ: "<<maximumZ<<std::endl;*/
   	
   bottomLeft.x = minimumX-(NormVoxelRadius+theSpatiocyteStepper->getColLength());
   bottomLeft.y = minimumY-(NormVoxelRadius+theSpatiocyteStepper->getLayerLength());
   bottomLeft.z = minimumZ-(NormVoxelRadius+theSpatiocyteStepper->getRowLength());
   topRight.x = maximumX+(NormVoxelRadius+theSpatiocyteStepper->getColLength());
   topRight.y = maximumY+(NormVoxelRadius+theSpatiocyteStepper->getLayerLength());
   topRight.z = maximumZ+(NormVoxelRadius+theSpatiocyteStepper->getRowLength());
   if(topRight.x>theComp->lengthX)topRight.x=theComp->lengthX;
   if(topRight.y>theComp->lengthY)topRight.y=theComp->lengthY;
   if(topRight.z>theComp->lengthZ)topRight.z=theComp->lengthZ;
   if(bottomLeft.x<0)bottomLeft.x=0;
   if(bottomLeft.y<0)bottomLeft.y=0;
   if(bottomLeft.z<0)bottomLeft.z=0;
   std::cout<<"bottomLeft.x: "<<bottomLeft.x<<std::endl;
   std::cout<<"bottomLeft.y: "<<bottomLeft.y<<std::endl;
   std::cout<<"bottomLeft.z: "<<bottomLeft.z<<std::endl;
   std::cout<<"topRight.x: "<<topRight.x<<std::endl;
   std::cout<<"topRight.y: "<<topRight.y<<std::endl;
   std::cout<<"topRight.z: "<<topRight.z<<std::endl;

  unsigned blRow(0);
  unsigned blLayer(0);
  unsigned blCol(0);
  theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
  unsigned trRow(0);
  unsigned trLayer(0);
  unsigned trCol(0);
  theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
   std::cout<<blRow<<std::endl;
   std::cout<<blLayer<<std::endl;
   std::cout<<blCol<<std::endl;
   std::cout<<trRow<<std::endl;
   std::cout<<trLayer<<std::endl;
   std::cout<<trCol<<std::endl;

for(unsigned d(blRow); d <= trRow; ++d)
    {
      for(unsigned e(blLayer); e <= trLayer; ++e)
        {
          for(unsigned f(blCol); f <= trCol; ++f)
            {
              unsigned m(theSpatiocyteStepper->global2coord(d, e, f));

              //std::cout<<"Lattice Index: "<<m<<std::endl;
              Point n(theSpatiocyteStepper->coord2point(m));
              AB = sub(node[1],node[0]);
              AC = sub(node[2],node[0]);
              surfaceNormal = cross(AB,AC);
              surfaceNormal = norm(surfaceNormal);
              surfaceDisplace = dot(surfaceNormal, node[0]);
              //isOnAboveSurface(n);   
              if(isOnAboveSurface(n))
		{ 
              	Voxel& aVoxel((*theLattice)[m]);

		//std::cout<<o.x<<std::endl;
		//std::cout<<o.y<<std::endl;
		//std::cout<<o.z<<std::endl;
		for (unsigned i(0); i!=theAdjoiningCoordSize; i++)
			{
			unsigned coord(aVoxel.adjoiningCoords[i]);           	
			//std::cout<<m<<"\t"<<coord<<std::endl;
			Point o(theSpatiocyteStepper->coord2point(coord));
                        //isOnAboveSurface(o);  
			if(isOnAboveSurface(o)==false)
			{
			surfaceCoords.push_back(m); 
                        break;
			}
			}
	
		}
           }
        }
    }
for(unsigned i(0);i<corn.size();i++)
	{
	theVacantCompSpecies[0]->addMolecule(&(*theLattice)[corn[i]]);
	}

for (unsigned count(0);count<surfaceCoords.size();count++)
	{
	std::cout<<count<<"\t"<<surfaceCoords[count]<<std::endl;
     	theVacantSpecies->addCompVoxel(surfaceCoords[count]);
	}
theVacantSpecies->setIsPopulated();
theVacantCompSpecies[0]->setIsPopulated();


   std::cout<<(*theLattice).size()<<std::endl;
   //std::cout<<ipf<<std::endl;
   idot1=strchr(ipf,'.');
   //std::cout<<idot1<<std::endl; 
   idot=(idot1-ipf)+1;
   //std::cout<<idot<<std::endl;
   if(idot<=0 || idot>20)
   std::cout<<"invalid inputfile name"<<std::endl;

   extractnumeric(ipf,idot);
   openfile(ipf,id1,idot);
   logInt = 0;
   delt=dscp[idfrm-1][2];
   logInt=dscp[idfrm-1][3];
   //cout<<delt<<"\n"<<logInt<<endl;
   timedump(delt,logInt);
   //assignQuad();
   
   idphi=0; 
   idvis=0; 
   idvfx=0; 
   idgam=0; 
   idpsi=0; 
   idsfr=0; 
   idbfr=0; 
   iddrg=0; 
   idtrc=0; 
   idhyc=0; 
   initsvec(); 
   openunit12();
   idebug=1;
        //for(int a=0;a<5;a++){
   initarea(area);
   //std::cout<<"time1: "<<time_<<std::endl;
   //std::cout<<"tstp1: "<<tstp<<std::endl;
   clchm(area);
   //std::cout<<"time2: "<<time_<<std::endl;
   //std::cout<<"tstp2: "<<tstp<<std::endl;
   cmstpsiz(cmdt,id3);
   //std::cout<<"time3: "<<time_<<std::endl;
   //std::cout<<"tstp3: "<<tstp<<std::endl;
   tstp=std::min(cmdt[1],cmdt[3]);
   //std::cout<<"time4: "<<time_<<std::endl;
   //std::cout<<"tstp4: "<<tstp<<std::endl;
   if(tstp>=(tnext-time_))
   {
	tstp = tnext-time_;
	isve = 1;
   }
   else isve=0;
   //std::cout<<"time5: "<<time_<<std::endl;
   //std::cout<<"tstp5: "<<tstp<<std::endl;
   dfdriver(eul);
   //std::cout<<"time6: "<<time_<<std::endl;
   //std::cout<<"tstp6: "<<tstp<<std::endl;
   time_=time_+tstp;
   //std::cout<<"time7: "<<time_<<std::endl;
   //std::cout<<"tstp7: "<<tstp<<std::endl;

   	//}
  
    }

  void MechanicsProcess::initializeFifth()
    {
	theInterval = tstp;
   	theTime= time_; 
   	thePriorityQueue->move(theQueueID);
    }

  void MechanicsProcess::fire()
    {

   chksurmol();
   //std::cout<<"time8: "<<time_<<std::endl;
   //std::cout<<"tstp8: "<<tstp<<std::endl;
   wrfile(id2,isve,ipf,idot,delt,logInt);
   //std::cout<<"time9: "<<time_<<std::endl;
   //std::cout<<"tstp9: "<<tstp<<"\n\n"<<std::endl;
   initarea(area);
   clchm(area);
   cmstpsiz(cmdt,id3);
   tstp=std::min(cmdt[1],cmdt[3]);
   if(tstp>=(tnext-time_))
   {
	tstp = tnext-time_;
	isve = 1;
   }
   else isve=0;
   dfdriver(eul);
   time_=time_+tstp;
   theInterval = tstp;
   theTime= time_; 
   thePriorityQueue->moveTop();
/*
      logSpecies();
      theLogFile.flush();
      if(LogInterval > 0)
        {
          theTime += LogInterval;
          thePriorityQueue->moveTop();
        }
      else
        {
          //get the next step interval of the SpatiocyteStepper:
          double aTime(theTime);
          theTime = libecs::INF;
          thePriorityQueue->moveTop();
          if(thePriorityQueue->getTop() != 
             dynamic_cast<SpatiocyteProcess*>(this))
            {
              theInterval = thePriorityQueue->getTop()->getTime() -
                theSpatiocyteStepper->getCurrentTime();
              setPriority(thePriorityQueue->getTop()->getPriority()-1);
            }
          theTime = aTime + theInterval;
          thePriorityQueue->move(theQueueID);
        }
*/
    }

bool MechanicsProcess::isOnAboveSurface(Point& aPoint)
{
  double disp(point2planeDisp(aPoint, surfaceNormal, surfaceDisplace));
  if(disp >= 0)
    {
      return true;
    }
  return false;
}

/*void MechanicsProcess::assignQuad()
{
	quadIndex.resize(nq);
        std::cout<<"quadsize: "<<quadIndex.size()<<std::endl;
	//for(unsigned i(0);i<quadIndex.size();i++)	
	//quadIndex[i]	
}*/

}
