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
   ipf = new char[20];
   for(int i(0);i<FileName.size();i++)
   {
     ipf[i] = FileName[i];
   }
   ipf[FileName.size()] = '\0';
   id1 = 11;
   id2 = 12;
   id3 = 1e-2;
   cmdt = new double[12];
   eul = new char[10];
   eul = "eulerian";
   fixsurfaceDisplace = 0;
   voxelRadius = theSpatiocyteStepper->getVoxelRadius();
   normVoxelRadius = theSpatiocyteStepper->getNormalizedVoxelRadius();  
   Comp* theComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
   lengthX = theComp->lengthX;
   lengthY = theComp->lengthY;
   lengthZ = theComp->lengthZ;

   
   idot1=strchr(ipf,'.'); 
   idot=(idot1-ipf)+1;
   if(idot<=0 || idot>20)
   std::cout<<"invalid inputfile name"<<std::endl;

   extractnumeric(ipf,idot);
   openfile(ipf,id1,idot);
   logInt = 0;
   delt=dscp[idfrm-1][2];
   logInt=dscp[idfrm-1][3];
   timedump(delt,logInt);
   assignQuad();
   assignNeigh();
   fitMechanotoSpatio();
   for (int i(0);i<quadIndex.size();i++)
    {
     getBLTR(i);
     getSurfaceCoords();
    }
   populateSurface();

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
   idebug=1;
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

   wrfile(id2,isve,ipf,idot,delt,logInt);
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

    }

void MechanicsProcess::assignQuad()
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


void MechanicsProcess::assignNeigh()
{
	neigh.resize(nq);
	for (unsigned a=0;a<nq;a++)
	{
	  neigh[a].reserve(15);
	}

	for(int i(0);i<quadIndex.size();i++)	
	{		
	  for (int j(0);j<4;j++)
	    {
		for (int k(0);k<quadIndex.size();k++)	
		{
		  for (int l(0);l<4;l++)
		    {
			int nodeID(quadIndex[k][l]);
			if(nodeID==quadIndex[i][j])
			{
			  std::vector<int>::iterator it(std::find(neigh[i].begin(),
			  neigh[i].end(), k));	
			  if(k!=i && it==neigh[i].end())
			  neigh[i].push_back(k);  
			}
		    }
		}
	    }
	}

}


void MechanicsProcess::fitMechanotoSpatio()
{

	minhvecX=hvec[0][0][0];
	minhvecY=hvec[0][0][1];
	minhvecZ=hvec[0][0][2];
	maxhvecX=hvec[0][0][0];
	maxhvecY=hvec[0][0][1];
	maxhvecZ=hvec[0][0][2];

	for (int i(0);i<ns;i++)
	  {
	    for (int j(0);j<3;j++)
	    	{
			minhvecX=std::min(hvec[i][j][0],minhvecX);
			minhvecY=std::min(hvec[i][j][1],minhvecY);
			minhvecZ=std::min(hvec[i][j][2],minhvecZ);
			maxhvecX=std::max(hvec[i][j][0],maxhvecX);
			maxhvecY=std::max(hvec[i][j][1],maxhvecY);
			maxhvecZ=std::max(hvec[i][j][2],maxhvecZ);
		    
		}
	  }

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
	double minRatio(std::min(ratioX,(std::min(ratioY,ratioZ))));
	double translate(4*voxelRadius);

	for (int i(0);i<ns;i++)
	  {
	    for (int j(0);j<3;j++)
	    	{
		realHvec[i][j][0] = hvec[i][j][0]+(-minhvecX);
		realHvec[i][j][0]=(realHvec[i][j][0]*minRatio+translate)/(2*voxelRadius);

		realHvec[i][j][1] = hvec[i][j][1]+(-minhvecY);
		realHvec[i][j][1]=(realHvec[i][j][1]*minRatio+translate)/(2*voxelRadius);

		realHvec[i][j][2] = hvec[i][j][2]+(-minhvecZ);
		realHvec[i][j][2]=(realHvec[i][j][2]*minRatio+translate)/(2*voxelRadius);
		}
	  }
}


void MechanicsProcess::getBLTR(int i)
{

	newNode.resize(4);
	newNode[0].x=realHvec[quadIndex[i][0]-1][2][0];
	newNode[0].y=realHvec[quadIndex[i][0]-1][2][1];
	newNode[0].z=realHvec[quadIndex[i][0]-1][2][2];
	newNode[1].x=realHvec[quadIndex[i][1]-1][2][0];
	newNode[1].y=realHvec[quadIndex[i][1]-1][2][1];
	newNode[1].z=realHvec[quadIndex[i][1]-1][2][2];
	newNode[2].x=realHvec[quadIndex[i][2]-1][2][0];
	newNode[2].y=realHvec[quadIndex[i][2]-1][2][1];
	newNode[2].z=realHvec[quadIndex[i][2]-1][2][2];
	newNode[3].x=realHvec[quadIndex[i][3]-1][2][0];
	newNode[3].y=realHvec[quadIndex[i][3]-1][2][1];
	newNode[3].z=realHvec[quadIndex[i][3]-1][2][2];

   	row.resize(4);
   	col.resize(4);
  	lay.resize(4);
	corn.resize(4);

   	for(unsigned i(0);i<4;i++)
	{
   	theSpatiocyteStepper->point2global(newNode[i], row[i], col[i], lay[i]);
   	corn[i]=theSpatiocyteStepper->global2coord(row[i], col[i], lay[i]);
	}

	for(unsigned i(0);i<corn.size();i++)
		{
		theVacantCompSpecies[0]->addMolecule(&(*theLattice)[corn[i]]);
		}


   	double minimumX=newNode[0].x;
   	double minimumY=newNode[0].y;
   	double minimumZ=newNode[0].z;
  	double maximumX=newNode[0].x;
   	double maximumY=newNode[0].y;
   	double maximumZ=newNode[0].z;

  	for(int a=1;a<4;a++)
 	{
  	minimumX=std::min(newNode[a].x,minimumX);
  	minimumY=std::min(newNode[a].y,minimumY);
   	minimumZ=std::min(newNode[a].z,minimumZ);
   	maximumX=std::max(newNode[a].x,maximumX);
  	maximumY=std::max(newNode[a].y,maximumY);
 	maximumZ=std::max(newNode[a].z,maximumZ);   	
 	}


   	bottomLeft.x = minimumX-(normVoxelRadius+theSpatiocyteStepper->getColLength());
   	bottomLeft.y = minimumY-(normVoxelRadius+theSpatiocyteStepper->getLayerLength());
   	bottomLeft.z = minimumZ-(normVoxelRadius+theSpatiocyteStepper->getRowLength());
   	topRight.x = maximumX+(normVoxelRadius+theSpatiocyteStepper->getColLength());
   	topRight.y = maximumY+(normVoxelRadius+theSpatiocyteStepper->getLayerLength());
   	topRight.z = maximumZ+(normVoxelRadius+theSpatiocyteStepper->getRowLength());
   	if(topRight.x>lengthX)topRight.x=lengthX;
   	if(topRight.y>lengthY)topRight.y=lengthY;
   	if(topRight.z>lengthZ)topRight.z=lengthZ;
  	if(bottomLeft.x<0)bottomLeft.x=0;
  	if(bottomLeft.y<0)bottomLeft.y=0;
  	if(bottomLeft.z<0)bottomLeft.z=0;



   	


}


void MechanicsProcess::getSurfaceCoords()
{

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
		calculateSurfaceNormal(newNode[0],newNode[1],newNode[2]);
        	if(isOnAboveSurface(n) && isOnBelowSideSurface(newNode[0],newNode[1],n)
		&& isOnBelowSideSurface(newNode[1],newNode[2],n) 
		&& isOnBelowSideSurface(newNode[2],newNode[3],n)
		&& isOnBelowSideSurface(newNode[3],newNode[0],n))
		{ 	
              		Voxel& aVoxel((*theLattice)[m]);	
			for (unsigned i(0); i!=theAdjoiningCoordSize; i++)
	  		{
			unsigned coord(aVoxel.adjoiningCoords[i]);           	
			Point o(theSpatiocyteStepper->coord2point(coord));

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
}

void MechanicsProcess::calculateSurfaceNormal(Point& node1, Point& node2, Point& node3)
{
	AB = sub(node2,node1);
	AC = sub(node3,node1);
	fixsurfaceNormal = cross(AB,AC);
	fixsurfaceNormal = norm(fixsurfaceNormal);
	fixsurfaceDisplace = dot(fixsurfaceNormal, node1);		
}	

bool MechanicsProcess::isOnBelowSideSurface(Point& node1, Point& node2, Point& aPoint)
{
	Point sideSurfaceNormal;
	double sideSurfaceDisplace(0);
	AB = sub(node2,node1);
	sideSurfaceNormal = cross(AB, fixsurfaceNormal);
	sideSurfaceNormal = norm(sideSurfaceNormal);
	sideSurfaceDisplace = dot(sideSurfaceNormal, node1);
  	double disp(point2planeDisp(aPoint, sideSurfaceNormal, sideSurfaceDisplace));
  	if(disp-normVoxelRadius <= 0)
   	{
     	return true;
   	}
 	return false;
}

bool MechanicsProcess::isOnAboveSurface(Point& aPoint)
{
  double disp(point2planeDisp(aPoint, fixsurfaceNormal, fixsurfaceDisplace));
  if(disp >= 0)
    {
      return true;
    }
  return false;
}

void MechanicsProcess::populateSurface()
{
	for (unsigned count(0);count<surfaceCoords.size();count++)
		{
     		theVacantSpecies->addCompVoxel(surfaceCoords[count]);
		}
	theVacantSpecies->setIsPopulated();
	theVacantCompSpecies[0]->setIsPopulated();
}



}



