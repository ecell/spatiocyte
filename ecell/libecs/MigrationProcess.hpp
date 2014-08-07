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


#ifndef __MigrationProcess_hpp
#define __MigrationProcess_hpp

#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <libecs/SpatiocyteProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{
extern "C" void extractnumeric(char *ipf, int &idot);
extern "C" void openfile(char *ipf,int &idiv,int &idot);
extern "C" void timedump(double &delt);
extern "C" double dscp[][12];
extern "C" int idfrm,idphi,idvis,idvfx,idgam,idpsi,idsfr,idbfr,iddrg,idtrc,
       idhyc,nq,ns,nl;
extern "C" double vbig;
extern "C" double vfixv[][3];
extern "C" double tnext,time_,tstop,tstp;
extern "C" double svec[][3][12];
extern "C" double cxprm[][12];
extern "C" double cxval[][12];
extern "C" int iqol[];
extern "C" int isoq[][4];
extern "C" int isol[][2];
extern "C" double hvec[][3][12];
extern "C" double evec[][4];
extern "C" double arean[][3];
extern "C" void clphi();
extern "C" void clvis();
extern "C" void clpsi();
extern "C" void govolint(double cnode[][3],double &volint);
extern "C" void gosurfintn(double cnode[][3],double &surfintv,double &surfintd,
                           double &surfinte);
extern "C" void clgam(double &aratio);
extern "C" void vvec1();
extern "C" void clsfr();
extern "C" void modriver(int&,double&,int&);
extern "C" void openunit12();
extern "C" void avtstep(double&);
extern "C" void clchm();
extern "C" void bc();
extern "C" void cmstpsiz(double*, double&);
extern "C" void dfdriver(const char *lag);
extern "C" void avgridmo(const char *lag);
extern "C" void avfield(int&,const char *nw);
extern "C" void wrfile(int&,char *ipf,int&,double&,double&,double&,double&,
                       double&); 

  LIBECS_DM_CLASS(MigrationProcess, SpatiocyteProcess)
  { 
  public:
    LIBECS_DM_OBJECT(MigrationProcess, Process)
      {
        INHERIT_PROPERTIES(Process);
        PROPERTYSLOT_SET_GET(String, FileName);
        PROPERTYSLOT_SET_GET(Real, minhvecX);
        PROPERTYSLOT_SET_GET(Real, minhvecY);
        PROPERTYSLOT_SET_GET(Real, minhvecZ);
        PROPERTYSLOT_SET_GET(Real, maxhvecX);
        PROPERTYSLOT_SET_GET(Real, maxhvecY);
        PROPERTYSLOT_SET_GET(Real, maxhvecZ);
      }
    MigrationProcess() {}
    virtual ~MigrationProcess() {}
    SIMPLE_SET_GET_METHOD(String, FileName); 
    SIMPLE_SET_GET_METHOD(Real, minhvecX); 
    SIMPLE_SET_GET_METHOD(Real, minhvecY); 
    SIMPLE_SET_GET_METHOD(Real, minhvecZ); 
    SIMPLE_SET_GET_METHOD(Real, maxhvecX); 
    SIMPLE_SET_GET_METHOD(Real, maxhvecY); 
    SIMPLE_SET_GET_METHOD(Real, maxhvecZ); 
    virtual void prepreinitialize();
    virtual void preinitialize();
    virtual void initialize();
    virtual void initializeFirst();
    virtual void initializeThird();
    virtual void initializeFifth();
    virtual void fire();
    virtual bool isOnAboveSurface(Point&,Point&,double&);
    virtual bool isOnBelowSideSurface(Point&,Point&,Point&,Point&);
    virtual bool isOnBelowInternalPlane(Point&,Point&,Point&,Point&);
    virtual bool isAboveInternalPlane(Point&,Point&,Point&,Point&);
    void setScalingFactor();
    void constructComp(bool);
    void getBox(std::vector<Point>&,Point&,Point&);
    void setQuadVoxels(std::vector<Point>&,Point&,Point&,bool,unsigned);
    void calculateSurfaceNormal(Point&,Point&,Point&,Point&,double&);
    void initValue();
    void initVmaxCnwmin();
    void initAvdtTstp();
    void writeNewFile();
    std::vector<Point> getQuad(int,int);
    std::vector<Point> getEdgeQuad(int,int);
    std::vector<Point> getGaussQuad(Point&,Point&,Point&,Point&);
    void initForces();
    void updateComp();
    void getCompartmentLength();
    void setCenterPoint();
    void populateMolecules();
    void advectSurfaceMolecule(Species*,unsigned);
    void replaceMolecules(std::vector<unsigned>,std::vector<unsigned>,Voxel*,
                          Voxel*,Species*,unsigned);
    unsigned getSurfaceAdjCoord(unsigned);
    void optimizeSurfaceVoxel();
    unsigned getLatticeResizeCoord(unsigned);
    void removeOldSurfaceVoxel();
    void getNewSurfaceVoxel();
    void fireOptimizeSurfaceVoxel(Voxel&);
    void setVacantSpecies(unsigned,Point&,double,bool,unsigned);
    void initUpdateComp();
    void setArea();
    void createList();
    void setQuadList();
    void setEdgeList();
    void setLineCenter(Quad&);//,Point[]);
    void setCenQuad(Point&,Point&,Point&,Point&,Point&,Point&,Point&);
    void setConc(unsigned,unsigned);
    //bool countMol(std::vector<Point>&,unsigned,unsigned);
    unsigned countMol(std::vector<Point>&,unsigned,unsigned,unsigned,unsigned,
                      Point&,double);
    unsigned getConc(Quad&,unsigned);
    void setMolCount(unsigned,unsigned,unsigned,Point&,double);
    void setFunctionsForConc(unsigned,unsigned);
    void updateConc();
    //Point setCenQuad();
    virtual void printParameters();
  private:
    void updateArea();

  private:
    String FileName;
    char *idot1;
    char *ipf;
    char *lag;
    char *nw;
    int idot;
    int idebug;
    int isve;
    int id1;
    int id2;
    int id4;
    int icyc;
    std::vector<std::vector<std::vector<unsigned> > > vacantSpeciesCount;
    std::vector<unsigned> totalVacantSpeciesCount;
    //int bigQuadCount;
    double id3;
    double volint;
    double surfintv;
    double surfintd;
    double surfinte;
    double tsurf;
    double r0;
    double a0;
    double aratio;
    double avdt;
    double cnode[2048][3];
    double delt;
    double *cmdt;
    double minhvecX;
    double minhvecY;
    double minhvecZ;
    double maxhvecX;
    double maxhvecY;
    double maxhvecZ;
    double initminX;
    double initminY;
    double initminZ;
    double initmaxX;
    double initmaxY;
    double initmaxZ;
    double voxelRadius; 
    double normVoxelRadius; 
    double lengthX;
    double lengthY;
    double lengthZ;
    double epsl;
    double scalingFactor;
    double translate;
    double nodeArea[2048][3];
    //std::vector<std::vector<int> > nodeArea;
    std::vector<Quad> quadList;
    std::vector<Node> nodeList;
    //std::vector<Point> cenQuad;
    Point cenQuad;
    std::vector<std::vector<std::vector<std::vector<unsigned> > > > 
      vacantCompSpeciesCount;
    std::vector<std::vector<unsigned> > totalVacantCompSpeciesCount;
    //std::vector<std::vector<std::vector<unsigned> > > tmpSurface;
    std::vector<unsigned> surfaceCoords;
    std::vector<Species*> theVacantCompSpecies;
    std::vector<Variable*> theVacantCompVariables;
    std::vector<std::vector<int> > quadIndex;
    std::vector<std::vector<int> > edgeIndex;
    std::vector<std::vector<int> > neigh;
    std::vector<unsigned> check;
    //std::vector<std::vector<double> > quadArea;
    Comp* theComp;  
    Comp* theSuperComp;  
    Variable* theVacantVariable;
    Variable* theAddedVariable;
    Variable* theOverlapVariable;
    Species* theVacantSpecies;
    Species* theAddedSpecies;
    Species* theOverlapSpecies;
    System* theSystem;
  };
}


#endif /* __MigrationProcess_hpp */
