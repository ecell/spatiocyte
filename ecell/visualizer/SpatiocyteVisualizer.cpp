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


#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <gtkmm.h>
#include <gtkglmm.h>
#include <pango/pangoft2.h>

#ifdef G_OS_WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#endif

#include <png.h>
#include <zlib.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <gtkmm/main.h>
#include <gtkmm/table.h>
#include <gtkmm/window.h>
#include <gtkmm/ruler.h>
#include <gtkmm/drawingarea.h>
#include <netinet/in.h>
#include <SpatiocyteVisualizer.hpp>

#define PI 3.1415926535897932384626433832795028841971693993751
#define MAX_COLORS 20
#define PNG_NUM_MAX 9999999
const unsigned int GLScene::TIMEOUT_INTERVAL = 20;
const unsigned int SCREEN_WIDTH = 460;
const unsigned int SCREEN_HEIGHT = 360;

double hue2rgb( double a, double b, double h )
{
  if( h < 0 )
    {
      h += 1;
    }
  if( h > 1 )
    {
      h -= 1;
    }
  if( 6*h < 1 )
    {
      return (b-a)*6*h+a;
    }
  if( 2*h < 1 )
    {
      return b;
    }
  if( 3*h < 2 )
    {
      return a+(b-a)*((2./3 )-h)*6;
    }
  return a;
}

// hsl[0:1]
// rgb[0:255]
void hsl2rgb(double h, double l, float* r, float* g, float* b)
{
  double s(1);
  double a;
  double d;
  if( s == 0 )
    {
      *r = l;
      *g = l;
      *b = l;
    }
  else
    {
      if( l<0.5 ) 
        {
          d = l*(1+s);
        }
      else           
        {
          d = l+s-s*l;
        }
      a = 2*l-d; 
      *r = hue2rgb( a, d, h + ( 1./3 ) );
      *g = hue2rgb( a, d, h );
      *b = hue2rgb( a, d, h - ( 1./3 ) );
    } 
}

GLScene::GLScene(const Glib::RefPtr<const Gdk::GL::Config>& config,
                 const char* aBaseName)
: Gtk::GL::DrawingArea(config),
  theRotateAngle(5.0),
  isInvertBound(false),
  isShownSurface(false),
  is_playing_(false),
  is_forward_(true),
  show3DMolecule(false),
  showSurface(false),
  showTime(true),
  startRecord(false),
  theResetTime(0),
  xAngle(0),
  yAngle(0),
  zAngle(0),
  frame_cnt_(0),
  font_size(24),
  theMeanPointSize(0),
  thePngNumber(1),
  theScreenHeight(SCREEN_HEIGHT),
  theScreenWidth(SCREEN_WIDTH),
  is_mouse_rotate_(false),
  is_mouse_zoom_(false),
  is_mouse_pan_(false),
  mouse_drag_pos_x_(0),
  mouse_drag_pos_y_(0),
  mouse_drag_pos_z_(0),
  mouse_x_(0),
  mouse_y_(0),
  init_zoom_(1.8)
{
  add_events(Gdk::VISIBILITY_NOTIFY_MASK | Gdk::BUTTON_PRESS_MASK |
             Gdk::BUTTON_RELEASE_MASK | Gdk::POINTER_MOTION_MASK | 
             Gdk::POINTER_MOTION_HINT_MASK | Gdk::KEY_PRESS_MASK);
    /*
             GDK_EXPOSURE_MASK|
    GDK_POINTER_MOTION_MASK|GDK_BUTTON_MOTION_MASK|
    GDK_BUTTON_PRESS_MASK|GDK_BUTTON_RELEASE_MASK|
    GDK_KEY_PRESS_MASK|GDK_KEY_RELEASE_MASK|
    GDK_ENTER_NOTIFY_MASK|GDK_LEAVE_NOTIFY_MASK|
    GDK_FOCUS_CHANGE_MASK|GDK_STRUCTURE_MASK);
    */


  std::ostringstream aFileName;
  aFileName << aBaseName << std::ends;
  theFile.open( aFileName.str().c_str(), std::ios::binary );
  theFile.read((char*) (&theLatticeType), sizeof(theLatticeType));
  theFile.read((char*) (&theMeanCount), sizeof(theMeanCount));
  theFile.read((char*) (&theStartCoord), sizeof(theStartCoord));
  theFile.read((char*) (&theRowSize), sizeof(theRowSize));
  theFile.read((char*) (&theLayerSize), sizeof(theLayerSize));
  theFile.read((char*) (&theColSize), sizeof(theColSize));
  theFile.read((char*) (&min_point_), sizeof(min_point_));
  theFile.read((char*) (&max_point_), sizeof(max_point_));
  theFile.read((char*) (&theLatticeSpSize), sizeof(theLatticeSpSize));
  theFile.read((char*) (&thePolymerSize), sizeof(thePolymerSize));
  theFile.read((char*) (&theReservedSize), sizeof(theReservedSize));
  theFile.read((char*) (&theOffLatticeSpSize), sizeof(theOffLatticeSpSize));
  theFile.read((char*) (&theLogMarker), sizeof(theLogMarker));
  theFile.read((char*) (&theVoxelRadius), sizeof(theVoxelRadius));
  unsigned int aSourceSize(thePolymerSize);
  unsigned int aTargetSize(thePolymerSize);
  unsigned int aSharedSize(thePolymerSize);
  /*The following is the order of the logged species:
   * theLatticeSpSize : coords
   * aSourceSize : coords
   * aTargetSize : coords
   * aSharedSize : coords
   * theReservedSize : coords
   * thePolymerSize : points
   * theOffLatticeSpSize : points
   */
  theTotalLatticeSpSize = theLatticeSpSize+aSourceSize+aTargetSize+aSharedSize+
    theReservedSize;
  theTotalOffLatticeSpSize = thePolymerSize+theOffLatticeSpSize; 
  theTotalSpeciesSize = theTotalLatticeSpSize+theTotalOffLatticeSpSize;
  theSpeciesNameList = new char*[theTotalSpeciesSize];
  theRadii = new double[theTotalSpeciesSize];
  //Set up the names of normal lattice species:
  for(unsigned int i(0); i!=theLatticeSpSize; ++i)
    {
      unsigned int aStringSize;
      theFile.read((char*) (&aStringSize), sizeof(aStringSize));
      theSpeciesNameList[i] = new char[aStringSize];
      char* buffer;
      buffer = new char[aStringSize+1];
      theFile.read(buffer, aStringSize);
      buffer[aStringSize] = '\0';
      sscanf(buffer, "Variable:%s", theSpeciesNameList[i]);
      theFile.read((char*) (&theRadii[i]), sizeof(theRadii[i]));
      std::cout << theSpeciesNameList[i] << " radius:" <<
        theRadii[i] << std::endl;
    }
  //Set up the names of polymer species:
  //source : coord
  //target : coord
  //shared : coord
  //.... <- skip space for reserved species : coord
  //poly : points
  for(unsigned int i(theLatticeSpSize); i!=theLatticeSpSize+thePolymerSize; ++i)
    {
      unsigned int aPolySpeciesIndex;
      theFile.read((char*) (&aPolySpeciesIndex), sizeof(aPolySpeciesIndex));
      thePolySpeciesList.push_back(aPolySpeciesIndex);
      theSpeciesNameList[i] = new char[20];
      sprintf(theSpeciesNameList[i], "source");
      theFile.read((char*) (&theRadii[i]), sizeof(theRadii[i]));
      theSpeciesNameList[thePolymerSize+i] = new char[20];
      sprintf(theSpeciesNameList[thePolymerSize+i], "target");
      theRadii[thePolymerSize+i] = theRadii[i];
      theSpeciesNameList[thePolymerSize*2+i] = new char[20];
      sprintf(theSpeciesNameList[thePolymerSize*2+i], "shared");
      theRadii[thePolymerSize*2+i] = theRadii[i];
      theSpeciesNameList[thePolymerSize*3+theReservedSize+i] = new char[20];
      sprintf(theSpeciesNameList[thePolymerSize*3+theReservedSize+i], "poly");
      theRadii[thePolymerSize*2+theReservedSize+i] = theRadii[i];
    }
  //Set up the tmp names of reserved species:
  for(unsigned int i(theLatticeSpSize+thePolymerSize*3);
      i!=theLatticeSpSize+thePolymerSize*3+theReservedSize; ++i)
    {
      theSpeciesNameList[i] = new char[20];
      sprintf(theSpeciesNameList[i], "tmp %d", i);
      theRadii[i] = theRadii[0];
    }
  //Set up the names of off lattice species:
  for(unsigned int i(theLatticeSpSize+thePolymerSize*4+theReservedSize);
      i!=theLatticeSpSize+thePolymerSize*4+theReservedSize+theOffLatticeSpSize;
      ++i)
    {
      unsigned int aStringSize;
      theFile.read((char*) (&aStringSize), sizeof(aStringSize));
      theSpeciesNameList[i] = new char[aStringSize];
      char* buffer;
      buffer = new char[aStringSize+1];
      theFile.read(buffer, aStringSize);
      buffer[aStringSize] = '\0';
      sscanf(buffer, "Variable:%s", theSpeciesNameList[i]);
      theFile.read((char*) (&theRadii[i]), sizeof(theRadii[i]));
      std::cout << theSpeciesNameList[i] << " radius:" <<
        theRadii[i] << std::endl;
    }
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i)
    {
      theRadii[i] /= theVoxelRadius*2;
    }
  theOriCol = theStartCoord/(theRowSize*theLayerSize);
  theSpeciesColor = new Color[theTotalSpeciesSize];
  theSpeciesVisibility = new bool[theTotalSpeciesSize];
  theXUpBound = new unsigned int[theTotalSpeciesSize];
  theXLowBound = new unsigned int[theTotalSpeciesSize];;
  theYUpBound = new unsigned int[theTotalSpeciesSize];
  theYLowBound = new unsigned int[theTotalSpeciesSize];
  theZUpBound = new unsigned int[theTotalSpeciesSize];
  theZLowBound = new unsigned int[theTotalSpeciesSize];
  theCoords = new unsigned int*[theTotalLatticeSpSize];
  theMeanPoints = new Point[1];
  theFrequency = new unsigned int*[theLatticeSpSize];
  theMoleculeSize = new unsigned int[theTotalLatticeSpSize];
  for(unsigned int j(0); j!=theTotalSpeciesSize; ++j)
    {
      theSpeciesVisibility[j] = true; 
      theXUpBound[j] = 0;
      theXLowBound[j] = 0;
      theYUpBound[j] = theLayerSize;
      theYLowBound[j] = 0;
      theZUpBound[j] = theRowSize;
      theZLowBound[j] = 0;
    }
  for(unsigned int j(0); j!=theTotalLatticeSpSize; ++j)
    {
      theMoleculeSize[j] = 0;
      theCoords[j] = new unsigned int[1];
    }
  for(unsigned int j(0); j!=theLatticeSpSize; ++j )
    {
      theFrequency[j] = new unsigned int[1];
    }
  thePoints = new Point*[theTotalOffLatticeSpSize];
  theOffLatticeMoleculeSize = new unsigned int[theTotalOffLatticeSpSize];
  for(unsigned int j(0); j!=theTotalOffLatticeSpSize; ++j )
    {
      theOffLatticeMoleculeSize[j] = 0;
      thePoints[j] = new Point[1];
    }
  double hueInterval(1.0/double(theLatticeSpSize+theOffLatticeSpSize));
  double speciesLuminosity(0.4);
  double sourceLuminosity(0.6);
  double targetLuminosity(0.75);
  double sharedLuminosity(0.8);
  double polyLuminosity(0.3);
  /*
  theSpeciesColor[0].r = 0.9;
  theSpeciesColor[0].g = 0.9;
  theSpeciesColor[0].b = 0.9;
  */
  for(unsigned int i(0); i!=theLatticeSpSize; ++i)
    {
      hsl2rgb(hueInterval*i, speciesLuminosity,
              &theSpeciesColor[i].r,
              &theSpeciesColor[i].g,
              &theSpeciesColor[i].b);
    }
  for(unsigned int i(0); i!=thePolymerSize; ++i)
    {
      hsl2rgb(hueInterval*thePolySpeciesList[i],
              sourceLuminosity,
              &theSpeciesColor[theLatticeSpSize+i].r,
              &theSpeciesColor[theLatticeSpSize+i].g,
              &theSpeciesColor[theLatticeSpSize+i].b);
      hsl2rgb(hueInterval*thePolySpeciesList[i],
              targetLuminosity,
              &theSpeciesColor[theLatticeSpSize+thePolymerSize+i].r,
              &theSpeciesColor[theLatticeSpSize+thePolymerSize+i].g,
              &theSpeciesColor[theLatticeSpSize+thePolymerSize+i].b);
      hsl2rgb(hueInterval*thePolySpeciesList[i],
              sharedLuminosity,
              &theSpeciesColor[theLatticeSpSize+thePolymerSize*2+i].r,
              &theSpeciesColor[theLatticeSpSize+thePolymerSize*2+i].g,
              &theSpeciesColor[theLatticeSpSize+thePolymerSize*2+i].b);
      hsl2rgb(hueInterval*thePolySpeciesList[i],
              polyLuminosity,
              &theSpeciesColor[theTotalLatticeSpSize+i].r,
              &theSpeciesColor[theTotalLatticeSpSize+i].g,
              &theSpeciesColor[theTotalLatticeSpSize+i].b);
    }
  for(unsigned int i(0); i!=theOffLatticeSpSize; ++i)
    {
      hsl2rgb(hueInterval*(theLatticeSpSize+i), speciesLuminosity,
              &theSpeciesColor[theTotalLatticeSpSize+thePolymerSize+i].r,
              &theSpeciesColor[theTotalLatticeSpSize+thePolymerSize+i].g,
              &theSpeciesColor[theTotalLatticeSpSize+thePolymerSize+i].b);
    }
  hueInterval = 1.0/double(theReservedSize);
  speciesLuminosity = 0.6;
  for(unsigned int i(theTotalLatticeSpSize-theReservedSize);
      i!=theTotalLatticeSpSize; ++i )
    {
      hsl2rgb(hueInterval*(i-(theTotalLatticeSpSize-theReservedSize)),
              speciesLuminosity,
              &theSpeciesColor[i].r,
              &theSpeciesColor[i].g,
              &theSpeciesColor[i].b); 
      theSpeciesVisibility[i] = false;
    }
  /*
  //Set all as black:
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i)
    {
      theSpeciesColor[i].r = 0;
      theSpeciesColor[i].g = 0;
      theSpeciesColor[i].b = 0;
    }
    */
  std::cout << "row:" << theRowSize << " col:" << theColSize  <<
    " layer:" << theLayerSize << " marker:" <<
    theLogMarker << std::endl << std::flush;
  theOriRow = 0;
  theOriLayer = 0;
  theRadius = 0.5;
  switch(theLatticeType)
    {
    case HCP_LATTICE: 
      theHCPl = theRadius/sqrt(3); 
      theHCPy = theRadius*sqrt(3);
      theHCPx = theRadius*sqrt(8.0/3.0); // for division require .0
      if(theMeanCount)
        {
          thePlot3DFunction = &GLScene::plotMean3DHCPMolecules;
          thePlotFunction = &GLScene::plotMeanHCPPoints;
          theLoadCoordsFunction = &GLScene::loadMeanCoords;
        }
      else
        {
          thePlot3DFunction = &GLScene::plot3DHCPMolecules;
          thePlotFunction = &GLScene::plotHCPPoints;
          theLoadCoordsFunction = &GLScene::loadCoords;
        }
      break;
    case CUBIC_LATTICE:
      if(theMeanCount)
        {
          thePlot3DFunction = &GLScene::plotMean3DCubicMolecules;
          theLoadCoordsFunction = &GLScene::loadMeanCoords;
        }
      else
        {
          thePlotFunction = &GLScene::plotCubicPoints;
          thePlot3DFunction = &GLScene::plot3DCubicMolecules;
          theLoadCoordsFunction = &GLScene::loadCoords;
        }
      break;
    }
  init_frames();
  dimensions_.x = max_point_.x - min_point_.x;
  dimensions_.y = max_point_.y - min_point_.y;
  dimensions_.z = max_point_.z - min_point_.z;
  mid_point_.x = dimensions_.x/2 + min_point_.x;
  mid_point_.y = dimensions_.y/2 + min_point_.y;
  mid_point_.z = dimensions_.z/2 + min_point_.z;
  ViewSize = 1.05*sqrt(dimensions_.x*dimensions_.x+
                       dimensions_.y*dimensions_.y+
                       dimensions_.z*dimensions_.z);
  if(ViewSize==0)
    { 
      ViewSize=1.0;
    }
  FieldOfView=45;
  Xtrans=Ytrans=0;
  Near=-ViewSize/2.0;
  //Near=-ViewSize;
  Aspect=1.0;
  set_size_request(theScreenWidth, theScreenHeight);
  std::cout << "done" << std::endl;
}

unsigned GLScene::get_frame_size() {
  return frames_.size();
}

void GLScene::init_frames() {
  std::streampos aStreamPos;
  frames_.push_back(theFile.tellg());
  loadCoords(aStreamPos);
  int curr_frame(frame_cnt_);
  theFile.seekg(frames_[++frame_cnt_]);
  while((this->*theLoadCoordsFunction)(aStreamPos)) { 
    ++frame_cnt_;
  }
  theFile.clear();
  frame_cnt_ = curr_frame;
  theFile.seekg(frames_[1]);
  (this->*theLoadCoordsFunction)(aStreamPos);
}

bool GLScene::on_button_press_event(GdkEventButton* event) {
  return on_button_release_event(event);
}

bool GLScene::get_is_event_masked(GdkEventButton* event, int mask) {
  return (event->state & mask) == mask;
} 

bool GLScene::get_is_button(GdkEventButton* event, int button, int mask) {
  if(event->button == button) {
    return (event->type == GDK_BUTTON_PRESS);
  }
  return get_is_event_masked(event, mask); 
}

void GLScene::set_position(double x, double y, double& px, double& py,
                           double& pz) {
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  px = float(x-viewport[0])/float(viewport[2]);
  py = float(y-viewport[1])/float(viewport[3]);      
  px = left_ + px*(right_-left_);
  py = top_  + py*(bottom_-top_);
  pz = Near;
}

/*
bool GLScene::on_motion_notify_event(GdkEventMotion* event) {
  if(event && event->window) {
    const Glib::RefPtr<Gdk::Window> refWindow =
        Glib::wrap(event->window, true); // true == take_copy
    if(refWindow) {
      int x(0);
      int y(0);
      Gdk::ModifierType state(Gdk::ModifierType(0));
      const Glib::RefPtr<const Gdk::Device> device =
        Glib::wrap(event->device, true); // true == take_copy

      if(event->is_hint) {
        //refWindow->get_device_position(device, x, y, state);
        refWindow->get_pointer(x, y, state);
      }
      else {
        x = event->x;
        y = event->y;
        //state = event->state;
      }
      double dx(x-mouse_x_);
      double dy(y-mouse_y_);
      if(dx == 0 && dy == 0) {
        return true;
      }
      mouse_x_ = x;
      mouse_y_ = y;
      bool changed(false);
      if(is_mouse_zoom_) {
        double s(exp(float(dy)*0.01));
        glTranslatef(mid_point_.x, mid_point_.y, mid_point_.z);
        glScalef(s,s,s);
        glTranslatef(-mid_point_.x, -mid_point_.y, -mid_point_.z);
        changed = true;
      }
      else if(is_mouse_rotate_) {
        double ax(dx);
        double ay(dy);
        double az(0);
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        double angle(sqrt(ax*ax+ay*ay+az*az)/float(viewport[2]+1)*180.0);
        double m[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, m);
        double bx(m[0]*ax + m[3]*ay + m[6]*az);
        double by(m[1]*ax + m[4]*ay + m[7]*az);
        double bz(m[2]*ax + m[5]*ay + m[8]*az);
        glTranslatef(mid_point_.x, mid_point_.y, mid_point_.z);
        glRotatef(angle,bx,by,bz);
        glTranslatef(-mid_point_.x, -mid_point_.y, -mid_point_.z);
        changed = true;
      }
      else if(is_mouse_pan_) {
        double px, py, pz;
        set_position(x, y, px, py, pz);
        double m[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, m);
        glLoadIdentity();
        glTranslatef(px-mouse_drag_pos_x_,py-mouse_drag_pos_y_,
                     pz-mouse_drag_pos_z_);
        glMultMatrixd(m);
        mouse_drag_pos_x_ = px;
        mouse_drag_pos_y_ = py;
        mouse_drag_pos_z_ = pz;
        changed = true;
      }
      if(changed) {
        queue_draw();
      }
    }
  }
  return true;
}
*/


static void
invertMatrix(const GLdouble *m, GLdouble *out )
{

/* NB. OpenGL Matrices are COLUMN major. */
#define MAT(m,r,c) (m)[(c)*4+(r)]

/* Here's some shorthand converting standard (row,column) to index. */
#define m11 MAT(m,0,0)
#define m12 MAT(m,0,1)
#define m13 MAT(m,0,2)
#define m14 MAT(m,0,3)
#define m21 MAT(m,1,0)
#define m22 MAT(m,1,1)
#define m23 MAT(m,1,2)
#define m24 MAT(m,1,3)
#define m31 MAT(m,2,0)
#define m32 MAT(m,2,1)
#define m33 MAT(m,2,2)
#define m34 MAT(m,2,3)
#define m41 MAT(m,3,0)
#define m42 MAT(m,3,1)
#define m43 MAT(m,3,2)
#define m44 MAT(m,3,3)

   GLdouble det;
   GLdouble d12, d13, d23, d24, d34, d41;
   GLdouble tmp[16]; /* Allow out == in. */

   /* Inverse = adjoint / det. (See linear algebra texts.)*/

   /* pre-compute 2x2 dets for last two rows when computing */
   /* cofactors of first two rows. */
   d12 = (m31*m42-m41*m32);
   d13 = (m31*m43-m41*m33);
   d23 = (m32*m43-m42*m33);
   d24 = (m32*m44-m42*m34);
   d34 = (m33*m44-m43*m34);
   d41 = (m34*m41-m44*m31);

   tmp[0] =  (m22 * d34 - m23 * d24 + m24 * d23);
   tmp[1] = -(m21 * d34 + m23 * d41 + m24 * d13);
   tmp[2] =  (m21 * d24 + m22 * d41 + m24 * d12);
   tmp[3] = -(m21 * d23 - m22 * d13 + m23 * d12);

   /* Compute determinant as early as possible using these cofactors. */
   det = m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2] + m14 * tmp[3];

   /* Run singularity test. */
   if (det == 0.0) {
     std::cout << "det 0" << std::endl;
      /* printf("invert_matrix: Warning: Singular matrix.\n"); */
/*    memcpy(out,_identity,16*sizeof(double)); */
   }
   else {
      GLdouble invDet = 1.0 / det;
      /* Compute rest of inverse. */
      tmp[0] *= invDet;
      tmp[1] *= invDet;
      tmp[2] *= invDet;
      tmp[3] *= invDet;

      tmp[4] = -(m12 * d34 - m13 * d24 + m14 * d23) * invDet;
      tmp[5] =  (m11 * d34 + m13 * d41 + m14 * d13) * invDet;
      tmp[6] = -(m11 * d24 + m12 * d41 + m14 * d12) * invDet;
      tmp[7] =  (m11 * d23 - m12 * d13 + m13 * d12) * invDet;

      /* Pre-compute 2x2 dets for first two rows when computing */
      /* cofactors of last two rows. */
      d12 = m11*m22-m21*m12;
      d13 = m11*m23-m21*m13;
      d23 = m12*m23-m22*m13;
      d24 = m12*m24-m22*m14;
      d34 = m13*m24-m23*m14;
      d41 = m14*m21-m24*m11;

      tmp[8] =  (m42 * d34 - m43 * d24 + m44 * d23) * invDet;
      tmp[9] = -(m41 * d34 + m43 * d41 + m44 * d13) * invDet;
      tmp[10] =  (m41 * d24 + m42 * d41 + m44 * d12) * invDet;
      tmp[11] = -(m41 * d23 - m42 * d13 + m43 * d12) * invDet;
      tmp[12] = -(m32 * d34 - m33 * d24 + m34 * d23) * invDet;
      tmp[13] =  (m31 * d34 + m33 * d41 + m34 * d13) * invDet;
      tmp[14] = -(m31 * d24 + m32 * d41 + m34 * d12) * invDet;
      tmp[15] =  (m31 * d23 - m32 * d13 + m33 * d12) * invDet;

      memcpy(out, tmp, 16*sizeof(GLdouble));
   }

#undef m11
#undef m12
#undef m13
#undef m14
#undef m21
#undef m22
#undef m23
#undef m24
#undef m31
#undef m32
#undef m33
#undef m34
#undef m41
#undef m42
#undef m43
#undef m44
#undef MAT
}

bool GLScene::on_motion_notify_event(GdkEventMotion* event) {
  double x(0);
  double y(0);
  /*
  Gdk::ModifierType state(Gdk::ModifierType(0));
  if(event->is_hint) {
    event->window->get_pointer(x, y, state);
    //event->window->get_device_position(event->device, x, y, state);
  }
  else {
    x = event->x;
    y = event->y;
    state = event->state;
  }
  */
  x = event->x;
  y = event->y;
  double dx(x-mouse_x_);
  double dy(y-mouse_y_);
  if(dx == 0 && dy == 0) {
    return true;
  }
  mouse_x_ = x;
  mouse_y_ = y;
  bool changed(false);
  if(is_mouse_zoom_) {
    FieldOfView *= pow(1.05, dy/4.0);
    if(FieldOfView > 180) {
      FieldOfView = 180;
    }
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(FieldOfView,Aspect,Near,ViewSize+Near);
    glMatrixMode(GL_MODELVIEW);
    changed = true;
  }
  else if(is_mouse_rotate_) {
    double ax(dy);
    double ay(dx);
    double az(0);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double len(sqrt(ax*ax+ay*ay+az*az));
    double angle(len/float(viewport[2]+1)*180.0);
    GLdouble m[16];
    GLdouble inv[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, m);
    invertMatrix(m, inv);
    double bx(inv[0]*ax + inv[4]*ay + inv[8]*az);
    double by(inv[1]*ax + inv[5]*ay + inv[9]*az);
    double bz(inv[2]*ax + inv[6]*ay + inv[10]*az);

    xAngle += bx*angle/len;
    normalizeAngle(xAngle);
    yAngle += by*angle/len;
    normalizeAngle(yAngle);
    zAngle += bz*angle/len;
    normalizeAngle(zAngle);

    glTranslatef(mid_point_.x, mid_point_.y, mid_point_.z);
    glRotatef(angle,bx,by,bz);
    glTranslatef(-mid_point_.x, -mid_point_.y, -mid_point_.z);

    m_control_->setXangle(xAngle);
    m_control_->setYangle(yAngle);
    m_control_->setZangle(zAngle);
    isShownSurface = false;
    changed = true;
  }
  else if(is_mouse_pan_) {
    double px, py, pz;
    set_position(x, y, px, py, pz);
    double m[16];
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX, m);
    glLoadIdentity();
    glTranslatef(px-mouse_drag_pos_x_,py-mouse_drag_pos_y_,
                 pz-mouse_drag_pos_z_);
    glMultMatrixd(m);
    mouse_drag_pos_x_ = px;
    mouse_drag_pos_y_ = py;
    mouse_drag_pos_z_ = pz;
    changed = true;
  }
  if(changed) {
    queue_draw();
  }
  return true;
}


bool GLScene::on_button_release_event(GdkEventButton* event) {
  bool left(get_is_button(event, 1, GDK_BUTTON1_MASK));
  bool middle(get_is_button(event, 2, GDK_BUTTON2_MASK));
  bool right(get_is_button(event, 3, GDK_BUTTON3_MASK));
  is_mouse_rotate_ = left and not (middle or right);
  is_mouse_zoom_ = middle or (left and right);
  is_mouse_pan_ = right and get_is_event_masked(event, GDK_CONTROL_MASK);
  double x(event->x);
  mouse_x_ = x;
  double y(event->y);
  mouse_y_ = y;
  set_position(x, y, mouse_drag_pos_x_, mouse_drag_pos_y_, mouse_drag_pos_z_);
  queue_draw();
  return true;
}

bool GLScene::on_scroll_event(GdkEventScroll* scroll_event) {
  double s(0); 
  if(scroll_event->direction == GDK_SCROLL_UP) {
    zoomIn();
  }
  else {
    zoomOut();
  }
}

bool GLScene::on_key_press_event(GdkEventKey* key_event) {

}

GLScene::~GLScene()
{
}

void GLScene::setXUpBound(unsigned int aBound )
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      if(theSpeciesVisibility[i])
        {
          theXUpBound[i] = aBound;
        }
    }
  queue_draw();
}

void GLScene::setXLowBound(unsigned int aBound )
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      if(theSpeciesVisibility[i])
        {
          theXLowBound[i] = aBound;
        }
    }
  queue_draw();
}

void GLScene::setYUpBound(unsigned int aBound )
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      if(theSpeciesVisibility[i])
        {
          theYUpBound[i] = aBound;
        }
    }
  queue_draw();
}

void GLScene::setYLowBound(unsigned int aBound )
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      if(theSpeciesVisibility[i])
        {
          theYLowBound[i] = aBound;
        }
    }
  queue_draw();
}

void GLScene::setZUpBound(unsigned int aBound )
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      if(theSpeciesVisibility[i])
        {
          theZUpBound[i] = aBound;
        }
    }
  queue_draw();
}

void GLScene::setZLowBound(unsigned int aBound )
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      if(theSpeciesVisibility[i])
        {
          theZLowBound[i] = aBound;
        }
    }
  queue_draw();
}

void GLScene::set3DMolecule(bool is3D)
{
  show3DMolecule = is3D;
  queue_draw();
}

void GLScene::setInvertBound(bool isInvert)
{
  isInvertBound = isInvert;
  queue_draw();
}

void GLScene::setShowTime(bool isShowTime)
{
  showTime = isShowTime;
  queue_draw();
}

void GLScene::setShowSurface(bool isShowSurface)
{
  showSurface = isShowSurface;
  queue_draw();
}

void GLScene::setRecord(bool isRecord)
{
  std::cout << "starting to record frames..." << std::endl;
  startRecord = isRecord;
}

void GLScene::resetTime()
{
  theResetTime = theCurrentTime;
  queue_draw();
}

void GLScene::setSpeciesVisibility(unsigned int id, bool isVisible)
{
  theSpeciesVisibility[id] = isVisible;
  queue_draw();
}

bool GLScene::getSpeciesVisibility(unsigned int id)
{
  return theSpeciesVisibility[id];
}

void GLScene::setControlBox(ControlBox* aControl)
{
  m_control_ = aControl;
}

void GLScene::set_is_forward(bool is_forward)
{
  is_forward_ = is_forward;
}

Color GLScene::getSpeciesColor(unsigned int id)
{
  return theSpeciesColor[id];
}

void GLScene::setSpeciesColor(unsigned int id, Color aColor)
{
  theSpeciesColor[id] = aColor;
  queue_draw();
}

void GLScene::setBackgroundColor(Color aColor)
{
  glClearColor (aColor.r, aColor.g, aColor.b, 0);
  queue_draw();
}

char* GLScene::getSpeciesName(unsigned int id)
{
  return theSpeciesNameList[id];
}

void GLScene::on_realize()
{
  Gtk::GL::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::GL::Window> glwindow = get_gl_window();
  if (!glwindow->gl_begin(get_gl_context()))
    {
      return;
    }
  //background color3D:
  glClearColor (0, 0, 0, 0);
  //glClearColor (1, 1, 1, 0);
  glClearDepth (1);
  if(!theMeanCount)
    {
      glEnable(GL_DEPTH_TEST); //To darken molecules farther away
      glDepthFunc(GL_LESS); //To show the molecule only if it is nearer (less)
    }
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  if(!theMeanCount)
    {
      // This hint is for antialiasing
      glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST); 
    }
  glEnable(GL_TEXTURE_2D);
  glColorMaterial( GL_FRONT_AND_BACK, GL_DIFFUSE );
  glEnable(GL_COLOR_MATERIAL);
  glMatrixMode(GL_MODELVIEW);
  glTranslatef(-mid_point_.x,-mid_point_.y,-mid_point_.z); 
  if(theMeanCount)
    {
      //for GFP visualization:
      glClearColor (0, 0, 0, 0);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
  else
    {
      GLfloat LightAmbient[]= {0.,0.,0.,1.}; 
      GLfloat LightDiffuse[]= {1.,1.,1.,1.};
      GLfloat LightSpecular[]= {0.7,0.7,0.7,1.};
      GLfloat LightPosition[]= {1.,1.,1.,0.};
      GLfloat MaterialAmbient[]= {.7,.7,.7,1.}; 
      GLfloat MaterialDiffuse[]= {.8,.8,.8,1.};
      GLfloat MaterialSpecular[]= {1.,1.,1.,1.};
      GLfloat Shininess(100.0);
      //for 3D molecules:
      glLightfv(GL_LIGHT0,GL_AMBIENT, LightAmbient);
      glLightfv(GL_LIGHT0,GL_DIFFUSE, LightDiffuse);
      glLightfv(GL_LIGHT0,GL_SPECULAR, LightSpecular);
      glLightfv(GL_LIGHT0,GL_POSITION, LightPosition);
      glMaterialfv(GL_FRONT,GL_AMBIENT, MaterialAmbient);
      glMaterialfv(GL_FRONT,GL_DIFFUSE, MaterialDiffuse);
      glMaterialfv(GL_FRONT,GL_SPECULAR, MaterialSpecular);
      glMaterialf(GL_FRONT,GL_SHININESS, Shininess);
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glDepthFunc(GL_LESS);
      glEnable(GL_DEPTH_TEST);
      glEnable(GL_NORMALIZE);
      glEnable(GL_COLOR_MATERIAL);

      /*
      glEnable(GL_LIGHTING);
      GLfloat LightAmbient[]= { 0.8, 0.8, 0.8, 1 }; 
      GLfloat LightDiffuse[]= { 1, 1, 1, 1 };
      //GLfloat LightPosition[]= { theLayerSize/2, theRowSize/2, theColSize, 1 };
      glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
      //glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
      glEnable(GL_LIGHT0);
      /*
      GLfloat   whiteSpecular[] = { 1.0, 1.0, 1.0, 1.0 };
      glMaterialfv( GL_FRONT, GL_SPECULAR, whiteSpecular );
      glMaterialf( GL_FRONT, GL_SHININESS, 50.0 );
      */
    }
  GLUquadricObj* qobj = gluNewQuadric();
  gluQuadricDrawStyle(qobj, GLU_FILL);
  theGLIndex = glGenLists(theTotalSpeciesSize+1);
  for(unsigned int i(theGLIndex); i != theTotalSpeciesSize+theGLIndex; ++i)
    {
      glNewList(i, GL_COMPILE);
      if(!theMeanCount)
        {
          gluSphere(qobj, theRadii[i-theGLIndex], 30, 30);
        }
      else
        {
          gluSphere(qobj, theRadii[i-theGLIndex], 10, 10);
        }
      glEndList();
    }
  PangoFontMap *fontmap;
  fontmap = pango_ft2_font_map_new();
  pango_ft2_font_map_set_resolution (PANGO_FT2_FONT_MAP(fontmap), 72, 72);
  ft2_context = Glib::wrap(pango_font_map_create_context(fontmap));
  /*
  glNewList(theTotalSpeciesSize+theGLIndex, GL_COMPILE);
  drawBox(0,theRealColSize,0,theRealLayerSize,0,theRealRowSize);
  glEndList();
  */
  glwindow->gl_end();
  step();
}


void GLScene::renderLayout(Glib::RefPtr<Pango::Layout> layout)
{
	unsigned char* begin_bitmap_buffer;
  Pango::Rectangle pango_extents = layout->get_pixel_logical_extents();
  pixel_extent_width = pango_extents.get_width();
  pixel_extent_height = pango_extents.get_height(); 

  FT_Bitmap bitmap;
  bitmap.rows = pixel_extent_height;
  bitmap.width = pixel_extent_width;
  bitmap.pitch = -bitmap.width;
  begin_bitmap_buffer = new unsigned char[bitmap.rows*bitmap.width];
  memset(begin_bitmap_buffer, 0, bitmap.rows * bitmap.width );
  bitmap.buffer = begin_bitmap_buffer + ( bitmap.rows - 1 ) * bitmap.width;
  bitmap.num_grays = 256;
  bitmap.pixel_mode = ft_pixel_mode_grays; 
  pango_ft2_render_layout_subpixel( &bitmap, layout->gobj(), 
                                    -pango_extents.get_x(), 0); 

  GLfloat bg[4];
  glGetFloatv(GL_COLOR_CLEAR_VALUE, bg); 
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glPixelTransferf(GL_RED_BIAS,  (1.0-bg[0]));
  glPixelTransferf(GL_GREEN_BIAS,  (1.0-bg[1]));
  glPixelTransferf(GL_BLUE_BIAS, (1.0-bg[2]));
  glPixelTransferf(GL_ALPHA_SCALE, 1 );
  glDrawPixels(bitmap.width, bitmap.rows,
               GL_ALPHA, GL_UNSIGNED_BYTE, begin_bitmap_buffer);
  //reset
  if(!theMeanCount)
    {
      glDisable(GL_BLEND);
    }
  glPixelTransferf(GL_RED_BIAS, 0.0f);
  glPixelTransferf(GL_GREEN_BIAS, 0.0f);
  glPixelTransferf(GL_BLUE_BIAS, 0.0f);
  delete[] begin_bitmap_buffer;
}


bool GLScene::on_expose_event(GdkEventExpose* event)
{
  Glib::RefPtr<Gdk::GL::Window> glwindow = get_gl_window();
  if (!glwindow->gl_begin(get_gl_context()))
    {
      return false;
    }

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if(theMeanCount)
    {
      glDisable(GL_LIGHTING);
      if(show3DMolecule)
        {
          (this->*thePlot3DFunction)();
        }
      else
        {
          (this->*thePlotFunction)();
        }
    }
  else if(show3DMolecule)
    {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      (this->*thePlot3DFunction)();
    }
  else
    {
      glDisable(GL_LIGHTING);
      (this->*thePlotFunction)();
    }
  if(showSurface && !isShownSurface) {
    resetView();
    rotateMidAxisAbs(90, 0, 1, 0);
    isShownSurface = true;
    }
  else if(!showSurface && isShownSurface) {
    resetView();
    rotateMidAxisAbs(0, 0, 1, 0);
    isShownSurface = false;
  }

  if(showTime)
    {
      drawTime();
    }

  glwindow->swap_buffers();
  glwindow->gl_end();
  return true;
}

void GLScene::drawTime()
{
 // Get the right font
  Glib::RefPtr<Pango::Context> widget_context = get_pango_context();
  Pango::FontDescription font_desc = widget_context->get_font_description();
  font_desc.set_size( font_size * PANGO_SCALE);
  ft2_context->set_font_description( font_desc); 
  
  // Compute the layout
  Glib::RefPtr<Pango::Layout> layout = Pango::Layout::create(ft2_context);
  layout->set_width(PANGO_SCALE*get_allocation().get_width());
  layout->set_alignment(Pango::ALIGN_LEFT);
  double aTime(theCurrentTime-theResetTime);
  std::stringstream ss;
  if(fabs(aTime) < 1e-3)
    {
      ss << aTime*1e+6 << " us";
    }
  else if(fabs(aTime) < 1)
    {
      ss << aTime*1e+3 << " ms";
    }
  else if(fabs(aTime) < 60)
    {
      ss << aTime << " s";
    }
  else if(fabs(aTime) < 3600)
    {
      ss << int(aTime/60) << "m " << int(aTime)%60  << "s";
    }
  else if(fabs(aTime) < 86400)
    {
      ss << int(aTime/3600) << "h " << int(aTime)%3600/60 << "m";
        //"m " << int(aTime)%3600%60  << "s";
    }
  else
    {
      ss << int(aTime/86400) << "d " << int(aTime)%86400/3600 << "h " 
        << int(aTime)%86400%3600/60 << "m";
        //"m " << int(aTime)%86400%3600%60  << "s";
    }
  layout->set_text(ss.str().c_str());

  Pango::Rectangle pango_extents = layout->get_pixel_logical_extents();
  pixel_extent_width = pango_extents.get_width();
  pixel_extent_height = pango_extents.get_height(); 
  GLfloat text_h;

  /* Text position */
  text_h = pixel_extent_height;

  GLfloat w(get_width());
  GLfloat h(get_height());
  unsigned int screenWidth(static_cast<GLsizei>(w));
  unsigned int screenHeight(static_cast<GLsizei>(h));

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity(); 
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity(); 

  glOrtho( 0, screenWidth, 0, screenHeight, 0, 1 );
  glDisable(GL_LIGHTING);

  glColor3f(1, 1, 1);
  glRasterPos2i(10,screenHeight-text_h-3);
  renderLayout(layout);

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glEnable(GL_LIGHTING);
}

bool GLScene::writePng()
{
  char filename[256];
  char str[256];
  sprintf(str,"image%%0%ii.png",(int)log10(PNG_NUM_MAX)+1);
  sprintf(filename,str,thePngNumber);
  ++thePngNumber; 
  GLfloat w(get_width());
  GLfloat h(get_height());
  unsigned int screenWidth(static_cast<GLsizei>(w));
  unsigned int screenHeight(static_cast<GLsizei>(h));

  FILE *outFile;
  outFile = fopen(filename, "wb");
  if(outFile == NULL)
    {
      return false;
    } 
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                NULL, NULL, NULL);
  if(!png_ptr)
    {
      fclose(outFile);
      return false;
    } 
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      fclose(outFile);
      return false;
    } 
  if(setjmp(png_jmpbuf(png_ptr)))
    {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(outFile);
      return false;
   }
  png_init_io(png_ptr, outFile);
  /* set the zlib compression level */
  /*png_set_compression_level(png_ptr, Z_NO_COMPRESSION);*/
  png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);
  /* set other zlib parameters */
  png_set_compression_mem_level(png_ptr, 8);
  png_set_compression_strategy(png_ptr, Z_DEFAULT_STRATEGY);
  png_set_compression_window_bits(png_ptr, 15);
  png_set_compression_method(png_ptr, 8);
  png_set_compression_buffer_size(png_ptr, 8192);
  png_set_IHDR(png_ptr, info_ptr, screenWidth, screenHeight, 8,
               PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr); 
  GLubyte *image =
    (GLubyte*)malloc(screenWidth*screenHeight*sizeof(GLubyte)*3);
  if(!image)
    {
      return false;
    }
  glPixelStorei(GL_PACK_ALIGNMENT,1);
  glReadPixels(0, 0, screenWidth, screenHeight, GL_RGB,
               GL_UNSIGNED_BYTE, image); 
  GLubyte** rowPtrs = new GLubyte*[screenHeight];
  for(GLuint i = 0; i < screenHeight; i++)
    {
      rowPtrs[screenHeight - i - 1] = &(image[i * screenWidth * 3]);
    } 
  png_write_image(png_ptr, rowPtrs); 
  png_write_end(png_ptr, info_ptr); 
  png_destroy_write_struct(&png_ptr, &info_ptr); 
  free(rowPtrs);
  fclose(outFile);
  //delete []rowPtrs;
  return true;
}


void GLScene::drawBox(GLfloat xlo, GLfloat xhi, GLfloat ylo, GLfloat yhi,
                      GLfloat zlo, GLfloat zhi)
{
  glBegin(GL_LINE_STRIP);
  glVertex3f(xlo,ylo,zlo);
  glVertex3f(xlo,ylo,zhi);
  glVertex3f(xlo,yhi,zhi);
  glVertex3f(xlo,yhi,zlo);
  glVertex3f(xlo,ylo,zlo);
  glVertex3f(xhi,ylo,zlo);
  glVertex3f(xhi,yhi,zlo);
  glVertex3f(xhi,yhi,zhi);
  glVertex3f(xhi,ylo,zhi);
  glVertex3f(xhi,ylo,zlo);
  glEnd();

  glBegin(GL_LINES);
  glVertex3f(xlo,ylo,zhi);
  glVertex3f(xhi,ylo,zhi);
  glVertex3f(xlo,yhi,zhi);
  glVertex3f(xhi,yhi,zhi);
  glVertex3f(xlo,yhi,zlo);
  glVertex3f(xhi,yhi,zlo);
  glEnd();
}

/*
void GLScene::reset() {
  GLfloat w = get_width();
  GLfloat h = get_height();
  Aspect = w/h;
  glViewport(0, 0, static_cast<GLsizei>(w), static_cast<GLsizei>(h));
  FieldOfView=45;
  Xtrans=Ytrans=0;
  if(w>=h) Near=ViewSize/2.0/tan(FieldOfView*PI/180.0/2.0);
  else Near=ViewSize/2.0/tan(FieldOfView*Aspect*PI/180.0/2.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(FieldOfView,Aspect,Near,ViewSize+Near);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-mid_point_.x,-mid_point_.y,-mid_point_.z);
  glTranslatef(0,0,-ViewSize/2.0-Near);
}
*/

void GLScene::resetView()
{
  GLfloat w = get_width();
  GLfloat h = get_height();
  Aspect = w/h;
  glViewport(0, 0, static_cast<GLsizei>(w), static_cast<GLsizei>(h));
  FieldOfView=45;
  Xtrans=Ytrans=0;
  if(w>=h) Near=ViewSize/2.0/tan(FieldOfView*PI/180.0/2.0);
  else Near=ViewSize/2.0/tan(FieldOfView*Aspect*PI/180.0/2.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //FieldOfView /= init_zoom_;
  gluPerspective(FieldOfView,Aspect,Near,ViewSize+Near);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-mid_point_.x, -mid_point_.y, -mid_point_.z-ViewSize/2.0-Near);
  invalidate();
  xAngle = 0;
  yAngle = 0;
  zAngle = 0;
  m_control_->setXangle(xAngle);
  m_control_->setYangle(yAngle);
  m_control_->setZangle(zAngle);
  isShownSurface = false;
  double tangent(tan(FieldOfView*PI/180.0/2));
  double height(2*ViewSize*tangent);
  double width(height*Aspect);
  left_ = -width/2;
  right_ = width/2;
  bottom_ = -height/2;
  top_ = height/2;
  /*
  top_ = h/2;
  bottom_ = -top_;
  left_ = -w/2;
  right_ = -left_;
  */
  std::cout << "viewsize:" << ViewSize << " near:" << Near << " tangent:" << tangent << " w:" << width << " h:" << height << std::endl;
  std::cout << "l:" << left_ << " r:" << right_ << " b:" << bottom_ << " t:" << top_ << std::endl;
}

/*
void GLScene::configure() {
  GLfloat width(get_width());
  GLfloat height(get_height());
  glViewport(0, 0, static_cast<GLsizei>(width), static_cast<GLsizei>(height));
  double zoom_out(1.2);
  double ratio(width/height);
  if(width >= height) {
    top_ = mid_point_.y*zoom_out;
    right_ = mid_point_.x*zoom_out*width/height;
  }
  else {
    top_ = mid_point_.y*zoom_out*height/width;
    right_ = mid_point_.x*zoom_out;
  }
  bottom_ = -top_;
  left_ = -right_;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(left_, right_, bottom_, top_, -40, 40);
  glMatrixMode(GL_MODELVIEW);
}
*/

void GLScene::configure() {
  GLfloat w = get_width();
  GLfloat h = get_height();
  GLfloat nearold = Near;
  GLfloat m[16];
  Aspect = w/h;
  glViewport(0, 0, static_cast<GLsizei>(w), static_cast<GLsizei>(h));
  if(w>=h) Near=ViewSize/2.0/tan(FieldOfView*PI/180.0/2.0);
  else Near=ViewSize/2.0/tan(FieldOfView*Aspect*PI/180.0/2.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //FieldOfView /= init_zoom_;
  gluPerspective(FieldOfView,Aspect,Near,ViewSize+Near); 
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX,m);
  glLoadIdentity();
  glTranslatef(0,0,nearold-Near);
  glMultMatrixf(m);
  double tangent(tan(FieldOfView*PI/180.0/2));
  double height(2*ViewSize*tangent);
  double width(height*Aspect);
  left_ = -width/2;
  right_ = width/2;
  bottom_ = -height/2;
  top_ = height/2;
  /*
  top_ = h/2;
  bottom_ = -top_;
  left_ = -w/2;
  right_ = -left_;
  */
  std::cout << "viewsize:" << ViewSize << " near:" << Near << " tangent:" << tangent << " w:" << width << " h:" << height << std::endl;
  std::cout << "l:" << left_ << " r:" << right_ << " b:" << bottom_ << " t:" << top_ << std::endl;
}


bool GLScene::on_configure_event(GdkEventConfigure* event)
{
  Glib::RefPtr<Gdk::GL::Window> glwindow = get_gl_window();
  if (!glwindow->gl_begin(get_gl_context()))
    {
      return false;
    }
  configure();

  glwindow->gl_end();
  return false;
}

bool GLScene::loadCoords(std::streampos& begin_pos)
{
  begin_pos = theFile.tellg();
  if(theFile.read((char*) (&theCurrentTime), sizeof(theCurrentTime)).fail())
    {
      return false;
    }
  unsigned int anIndex;
  // Get the species index
  if(theFile.read((char*) (&anIndex), sizeof(anIndex)).fail())
    {
      return false;
    }
  //Read lattice species:
  while(anIndex != theLogMarker)
    {
      unsigned int aMoleculeSize(0);
      if(theFile.read((char*) (&aMoleculeSize), sizeof(unsigned int)).fail())
        {
          return false;
        }
      //First read the coords into a temporary holder, aCoords
      //Only if we have successfully read it, we store the coords
      //in theCoords:
      unsigned int* aCoords(NULL);
      if(aMoleculeSize)
        {
          aCoords = new unsigned int[aMoleculeSize];
          if(theFile.read((char*) (aCoords), 
                          sizeof(unsigned int)*aMoleculeSize).fail())
            {
              delete []aCoords;
              return false;
            }
        }
      delete []theCoords[anIndex];
      theMoleculeSize[anIndex] = aMoleculeSize;
      theCoords[anIndex] = new unsigned int[aMoleculeSize];
      for(unsigned int j(0); j != aMoleculeSize; ++j)
        {
          theCoords[anIndex][j] = aCoords[j];
        }
      if(aMoleculeSize)
        {
          delete []aCoords;
        }
      if(theFile.read((char*) (&anIndex), sizeof(anIndex)).fail())
        {
          return false;
        }
    }
  if(theFile.read((char*) (&anIndex), sizeof(anIndex)).fail())
    {
      return false;
    }
  //Read off lattice species:
  while(anIndex != theLogMarker)
    { 
      unsigned int aMoleculeSize(0);
      if(theFile.read((char*) (&aMoleculeSize), sizeof(unsigned int)).fail())
        {
          return false;
        }
      //First read the points into a temporary holder, aPoints
      //Only if we have successfully read it, we store the points
      //in thePoints:
      Point* aPoints(NULL);
      if(aMoleculeSize)
        { 
          aPoints = new Point[aMoleculeSize];
          if(theFile.read((char*) (aPoints),
                          sizeof(Point)*aMoleculeSize).fail())
            {
              delete []aPoints;
              return false;
            }
        }
      delete []thePoints[anIndex];
      theOffLatticeMoleculeSize[anIndex] = aMoleculeSize;
      thePoints[anIndex] = new Point[aMoleculeSize];
      for(unsigned int j(0); j != aMoleculeSize; ++j)
        {
          thePoints[anIndex][j] = aPoints[j];
        }
      if(aMoleculeSize)
        {
          delete []aPoints;
        }
      if(theFile.read((char*) (&anIndex), sizeof(anIndex)).fail())
        {
          return false;
        }
    }
  if(frame_cnt_ == int(frames_.size())-1)
    {
      frames_.push_back(theFile.tellg());
    }
  return true;
}

bool GLScene::loadMeanCoords(std::streampos& begin_pos)
{
  begin_pos = theFile.tellg();
  if(theFile.read((char*) (&theCurrentTime), sizeof(theCurrentTime)).fail())
    {
      return false;
    }
  unsigned int aMoleculeSize(0);
  if(theFile.read((char*) (&aMoleculeSize), sizeof(aMoleculeSize)).fail())
    {
      return false;
    }
  Point* aPoints(NULL);
  if(aMoleculeSize)
    {
      aPoints = new Point[aMoleculeSize];
      if(theFile.read((char*) (aPoints), sizeof(Point)*aMoleculeSize).fail())
        {
          delete []aPoints;
          return false;
        }
    }
  delete []theMeanPoints;
  theMeanPointSize = aMoleculeSize;
  theMeanPoints = new Point[aMoleculeSize];
  for(unsigned int j(0); j != aMoleculeSize; ++j)
    {
      theMeanPoints[j] = aPoints[j];
    }
  if(aMoleculeSize)
    {
      delete []aPoints;
    }
  for(unsigned int j(0); j != theLatticeSpSize; ++j)
    {
      delete []theFrequency[j];
      theFrequency[j] = new unsigned int[theMeanPointSize];
      if(theFile.read((char*) (theFrequency[j]), 
                      sizeof(unsigned int)*theMeanPointSize).fail())
        {
          return false;
        }
    }
  if(frame_cnt_ == int(frames_.size())-1)
    {
      frames_.push_back(theFile.tellg());
    }
  return true;
}


void GLScene::plotMeanHCPPoints()
{
  plotHCPPoints();
  glBegin(GL_POINTS);
  double x,y,z;
  for(unsigned int k(0); k != theMeanPointSize; ++k)
    {
      y = theMeanPoints[k].y;
      z = theMeanPoints[k].z;
      x = theMeanPoints[k].x;
      for(unsigned int j(0); j!=theLatticeSpSize; ++j)
        {
          if(theSpeciesVisibility[j])
            {
              Color clr(theSpeciesColor[j]);
              double intensity((double)(theFrequency[j][k])/
                               (double)(theMeanCount/4));
              //glColor3f(clr.r*intensity, clr.g*intensity, clr.b*intensity); 
              glColor4f(clr.r, clr.g, clr.b, intensity);
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glVertex3f(x, y, z);
                }
            }
        }
    }
  glEnd();
}

void GLScene::plotMean3DHCPMolecules()
{
  plot3DHCPMolecules();
  double x,y,z;
  for(unsigned int k(0); k != theMeanPointSize; ++k)
    {
      y = theMeanPoints[k].y;
      z = theMeanPoints[k].z;
      x = theMeanPoints[k].x;
      for(unsigned int j(0); j!=theLatticeSpSize; ++j)
        {
          if(theSpeciesVisibility[j])
            {
              Color clr(theSpeciesColor[j]);
              double intensity((double)(theFrequency[j][k])/
                               (double)(theMeanCount/4));
              //glColor3f(clr.r*intensity, clr.g*intensity, clr.b*intensity); 
              glColor4f(clr.r, clr.g, clr.b, intensity);
              //glColor4f(clr.r*intensity, clr.g*intensity, clr.b*intensity,
               //         0.5f); 
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+j);
                  glPopMatrix();
                }
            }
        }
    }
}

void GLScene::plotMean3DCubicMolecules()
{
  unsigned int col, layer, row;
  double x,y,z;
  for(unsigned int k(0); k != theMeanPointSize; ++k)
    {
      y = theMeanPoints[k].y;
      z = theMeanPoints[k].z;
      x = theMeanPoints[k].x;
      for(unsigned int j(0); j!=theLatticeSpSize; ++j)
        {
          if(theSpeciesVisibility[j])
            {
              Color clr(theSpeciesColor[j]);
              double intensity((double)(theFrequency[j][k])/
                               (double)(theMeanCount/4));
              //glColor3f(clr.r*intensity, clr.g*intensity, clr.b*intensity); 
              glColor4f(clr.r, clr.g, clr.b, intensity);
              //glColor4f(clr.r*intensity, clr.g*intensity, clr.b*intensity,
               //         0.5f); 
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+j);
                  glPopMatrix();
                }
            }
        }
    }
  for(unsigned int j(theLatticeSpSize); j!=theTotalLatticeSpSize; ++j )
    {
      if(theSpeciesVisibility[j])
        {
          Color clr(theSpeciesColor[j]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theMoleculeSize[j]; ++k )
            {
              col = theCoords[j][k]/(theRowSize*theLayerSize)-theOriCol; 
              layer =
                (theCoords[j][k]%(theRowSize*theLayerSize))/theRowSize;
              row =
                (theCoords[j][k]%(theRowSize*theLayerSize))%theRowSize;
              y = layer*2*theRadius;
              z = row*2*theRadius;
              x = col*2*theRadius; 
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+j);
                  glPopMatrix();
                }
            }
        }
    }
}

void GLScene::plot3DHCPMolecules()
{
  unsigned int col, layer, row;
  double x,y,z;
  for(int j(theTotalLatticeSpSize-1); j!=-1; --j )
    {
      if( theSpeciesVisibility[j] )
        {
          Color clr(theSpeciesColor[j]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theMoleculeSize[j]; ++k )
            {
              col = theCoords[j][k]/(theRowSize*theLayerSize)-theOriCol; 
              layer =
                (theCoords[j][k]%(theRowSize*theLayerSize))/theRowSize;
              row =
                (theCoords[j][k]%(theRowSize*theLayerSize))%theRowSize;
              y = (col%2)*theHCPl + theHCPy*layer;
              z = row*2*theRadius + ((layer+col)%2)*theRadius;
              x = col*theHCPx;
              bool isBound(col <= theXUpBound[j] &&
                           col >= theXLowBound[j] &&
                           layer <= theYUpBound[j] &&
                           layer >= theYLowBound[j] &&
                           row <= theZUpBound[j] &&
                           row >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+j);
                  glPopMatrix();
                }
            }
        }
    }
  for(int j(theTotalOffLatticeSpSize-1); j!=-1; --j )
    {
      int m(j+theTotalLatticeSpSize);
      if(theSpeciesVisibility[m])
        {
          Color clr(theSpeciesColor[m]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theOffLatticeMoleculeSize[j]; ++k )
            {
              x = (thePoints[j][k].x);
              y = (thePoints[j][k].y);
              z = (thePoints[j][k].z);
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+m);
                  glPopMatrix();
                }
            }
        }
    }
  /*
  glPushMatrix();
  glTranslatef(0,0,0);
  glCallList(theTotalSpeciesSize+theGLIndex);
  glPopMatrix();
  */
}

void GLScene::plot3DCubicMolecules()
{
  unsigned int col, layer, row;
  double x,y,z;
  for( unsigned int j(0); j!=theTotalLatticeSpSize; ++j )
    {
      if( theSpeciesVisibility[j] )
        {
          Color clr(theSpeciesColor[j]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theMoleculeSize[j]; ++k )
            {
              col = theCoords[j][k]/(theRowSize*theLayerSize)-theOriCol; 
              layer =
                (theCoords[j][k]%(theRowSize*theLayerSize))/theRowSize;
              row =
                (theCoords[j][k]%(theRowSize*theLayerSize))%theRowSize;
              y = layer*2*theRadius;
              z = row*2*theRadius;
              x = col*2*theRadius; 
              if(!( x <= theXUpBound[j] && x >= theXLowBound[j] &&
                  y <= theYUpBound[j] && y >= theYLowBound[j] &&
                  z <= theZUpBound[j] && z >= theZLowBound[j]))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+j);
                  glPopMatrix();
                }
            }
        }
    }
  for(int j(theTotalOffLatticeSpSize-1); j!=-1; --j )
    {
      int m(j+theTotalLatticeSpSize);
      if( theSpeciesVisibility[m] )
        {
          Color clr(theSpeciesColor[m]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theOffLatticeMoleculeSize[j]; ++k )
            {
              x = (thePoints[j][k].x);
              y = (thePoints[j][k].y);
              z = (thePoints[j][k].z);
              if(!( x <= theXUpBound[m] && x >= theXLowBound[m] &&
                  y <= theYUpBound[m] && y >= theYLowBound[m] &&
                  z <= theZUpBound[m] && z >= theZLowBound[m]))
                {
                  glPushMatrix();
                  glTranslatef(x,y,z);
                  glCallList(theGLIndex+m);
                  glPopMatrix();
                }
            }
        }
    }
}

void GLScene::plotHCPPoints()
{
  glBegin(GL_POINTS);
  unsigned int col, layer, row;
  double x,y,z;
  for( unsigned int j(0); j!=theTotalLatticeSpSize; ++j )
    {
      if( theSpeciesVisibility[j] )
        {
          Color clr(theSpeciesColor[j]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theMoleculeSize[j]; ++k )
            {
              col = theCoords[j][k]/(theRowSize*theLayerSize)-theOriCol; 
              layer =
                (theCoords[j][k]%(theRowSize*theLayerSize))/theRowSize;
              row =
                (theCoords[j][k]%(theRowSize*theLayerSize))%theRowSize;
              y = (col%2)*theHCPl + theHCPy*layer;
              z = row*2*theRadius + ((layer+col)%2)*theRadius;
              x = col*theHCPx;
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glVertex3f(x, y, z);
                }
            }
        }
    }
  for(int j(theTotalOffLatticeSpSize-1); j!=-1; --j )
    {
      int m(j+theTotalLatticeSpSize);
      if( theSpeciesVisibility[m] )
        {
          Color clr(theSpeciesColor[m]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theOffLatticeMoleculeSize[j]; ++k )
            {
              x = (thePoints[j][k].x);
              y = (thePoints[j][k].y);
              z = (thePoints[j][k].z);
              bool isBound(x <= theXUpBound[j] && x >= theXLowBound[j] &&
                           y <= theYUpBound[j] && y >= theYLowBound[j] &&
                           z <= theZUpBound[j] && z >= theZLowBound[j]);
              if((isInvertBound && isBound) || (!isInvertBound && !isBound))
                {
                  glVertex3f(x, y, z);
                }
            }
        }
    }
  glEnd();
}

void GLScene::plotCubicPoints()
{
  glBegin(GL_POINTS);
  unsigned int col, layer, row;
  double x,y,z;
  for( unsigned int j(0); j!=theTotalLatticeSpSize; ++j )
    {
      if( theSpeciesVisibility[j] )
        {
          Color clr(theSpeciesColor[j]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theMoleculeSize[j]; ++k )
            {
              col = theCoords[j][k]/(theRowSize*theLayerSize)-theOriCol; 
              layer =
                (theCoords[j][k]%(theRowSize*theLayerSize))/theRowSize;
              row =
                (theCoords[j][k]%(theRowSize*theLayerSize))%theRowSize;
              y = layer*2*theRadius;
              z = row*2*theRadius;
              x = col*2*theRadius; 
              if(!( x <= theXUpBound[j] && x >= theXLowBound[j] &&
                  y <= theYUpBound[j] && y >= theYLowBound[j] &&
                  z <= theZUpBound[j] && z >= theZLowBound[j]))
                {
                  glVertex3f(x, y, z);
                }
            }
        }
    }
  for(int j(theTotalOffLatticeSpSize-1); j!=-1; --j )
    {
      int m(j+theTotalLatticeSpSize);
      if( theSpeciesVisibility[m] )
        {
          Color clr(theSpeciesColor[m]);
          glColor3f(clr.r, clr.g, clr.b); 
          for( unsigned int k(0); k!=theOffLatticeMoleculeSize[j]; ++k )
            {
              x = (thePoints[j][k].x);
              y = (thePoints[j][k].y);
              z = (thePoints[j][k].z);
              if(!( x <= theXUpBound[m] && x >= theXLowBound[m] &&
                  y <= theYUpBound[m] && y >= theYLowBound[m] &&
                  z <= theZUpBound[m] && z >= theZLowBound[m]))
                {
                  glVertex3f(x, y, z);
                }
            }
        }
    }
  glEnd();
}

void GLScene::setLayerColor( unsigned int aLayer )
{

  GLfloat a( theLayerSize/3.0 );
  GLfloat aRes( 1.0/a );
  GLfloat r( aLayer*aRes );
  GLfloat g( 0.0 );
  GLfloat b( 0.0 );
  if( r > 1 )
    {
      g = r - 1;
      r = 1.0;
      if( g > 1 )
        {
          b = g - 1;
          g = 1.0;
         }
    }
  glColor3f(r,g,b);
}

void GLScene::setTranslucentColor( unsigned int i, GLfloat j )
{  
  switch(i)
    {
    case 0:
      glColor4f(0.2,0.2,0.2,j);
      break;
    case 1:
      glColor4f(0.9,0.2,0.0,j);
      break;
    case 2:
      glColor4f(0,0.9,0.2,j);
      break;
    case 3:
      glColor4f(0.5,0.5,0.5,j);
      break;
    case 4:
      glColor4f(0.2,0,0.9,j);
      break;
    case 5:
      glColor4f(0.9,0,0.8,j);
      break;
    case 6:
      glColor4f(0,0.7,0.8,j);
      break;
    case 7:
      glColor4f(0.5,0,0,j);
      break;
    case 8:
      glColor4f(0,0.5,0,j);
      break;
    case 9:
      glColor4f(0,0,0.5,j);
      break;
    case 10:
      glColor4f(0.5,0.5,0,j);
      break;
    default :
      glColor4f((float)rand()/RAND_MAX,
                (float)rand()/RAND_MAX,
                (float)rand()/RAND_MAX,j);
    }
}

void GLScene::rotate(int aMult, int x, int y, int z)
{
  glMatrixMode(GL_MODELVIEW);
  glRotatef(theRotateAngle*aMult,x,y,z);
  invalidate();
}

void GLScene::translate(int x, int y, int z)
{
  glMatrixMode(GL_MODELVIEW);
  glTranslatef(x,y,z);
  invalidate();
}


void GLScene::zoomIn()
{ 
  FieldOfView/=1.05;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(FieldOfView,Aspect,Near,ViewSize+Near);
  glMatrixMode(GL_MODELVIEW);
  invalidate();
}

void GLScene::zoomOut()
{
  FieldOfView*=1.05;
  if(FieldOfView>180)
    {
      FieldOfView=180;
    }
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(FieldOfView,Aspect,Near,ViewSize+Near);
  glMatrixMode(GL_MODELVIEW);
  invalidate();
}

int GLScene::get_frame_cnt() {
  return frame_cnt_;
}

void GLScene::set_frame_cnt(int frame_cnt) {
  if(is_playing_) {
    is_playing_ = false;
    timeout_remove();
  }
  if(is_forward_) {
    frame_cnt_ = std::max(frame_cnt-1, 0);
  }
  else {
    frame_cnt_ = std::min(frame_cnt+1, int(frames_.size())+1);
  }
  on_timeout();
}

bool GLScene::on_timeout()
{
  inc_dec_frame_cnt();
  theFile.seekg(frames_[frame_cnt_]);
  std::streampos aCurrStreamPos;
  if(!(this->*theLoadCoordsFunction)(aCurrStreamPos)) { 
    theFile.clear();
    theFile.seekg(aCurrStreamPos);
    frame_cnt_ = std::max(frame_cnt_-1, 1);
  }
  m_control_->set_frame_cnt(frame_cnt_);
  m_control_->setTime(theCurrentTime);
  invalidate();
  if(startRecord)
    {
      timeout_remove();
      writePng();
      std::cout << "wrote png:" << thePngNumber << std::endl; 
      timeout_add();
    }
  return true;
}

void GLScene::inc_dec_frame_cnt() {
  if(is_forward_) {
    frame_cnt_ = std::min(frame_cnt_+1, int(frames_.size()));
  }
  else {
    frame_cnt_ = std::max(frame_cnt_-1, 1);
  }
}

void GLScene::step()
{
  if(is_playing_) {
    is_playing_ = false;
    timeout_remove();
  }
  on_timeout();
}

void GLScene::timeout_add()
{
  if (!m_ConnectionTimeout.connected())
    m_ConnectionTimeout = Glib::signal_timeout().connect(
      sigc::mem_fun(*this, &GLScene::on_timeout), TIMEOUT_INTERVAL);
}

void GLScene::timeout_remove()
{
  if (m_ConnectionTimeout.connected())
    m_ConnectionTimeout.disconnect();
}

bool GLScene::on_map_event(GdkEventAny* event)
{
  if(is_playing_) {
    timeout_add();
  }
  return true;
}

bool GLScene::on_unmap_event(GdkEventAny* event)
{
  timeout_remove();
  return true;
}

bool GLScene::on_visibility_notify_event(GdkEventVisibility* event)
{
  if(is_playing_) {
    if (event->state == GDK_VISIBILITY_FULLY_OBSCURED) {
      timeout_remove();
    }
    else {
      timeout_add();
    }
  }
  return true;
}

void GLScene::resetBound()
{
  for(unsigned int i(0); i!=theTotalSpeciesSize; ++i )
    {
      theXUpBound[i] = 0;
      theXLowBound[i] = 0;
      theYUpBound[i] = theLayerSize;
      theYLowBound[i] = 0;
      theZUpBound[i] = theRowSize;
      theZLowBound[i] = 0;
    }
  queue_draw();
}

void GLScene::rotateMidAxis(int aMult, int x, int y, int z)
{
  GLfloat m[16];
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX,m);
  glLoadIdentity();
  glTranslatef(Xtrans,Ytrans,-(Near+ViewSize/2.0));
  glRotatef(theRotateAngle*aMult,x,y,z);
  glTranslatef(-Xtrans,-Ytrans,+(Near+ViewSize/2.0));
  glMultMatrixf(m);
  invalidate();
  if(x)
    {
      xAngle += aMult*theRotateAngle;
      normalizeAngle(xAngle);
      m_control_->setXangle(xAngle);
    }
  if(y)
    {
      yAngle += aMult*theRotateAngle;
      normalizeAngle(yAngle);
      m_control_->setYangle(yAngle);
    }
  if(z)
    {
      zAngle += aMult*theRotateAngle;
      normalizeAngle(zAngle);
      m_control_->setZangle(zAngle);
    }
}

void GLScene::rotateMidAxisAbs(double angle, int x, int y, int z)
{
  if(x && angle == xAngle) {
    return;
  }
  if(y && angle == yAngle) {
    return;
  }
  if(z && angle == zAngle) {
    return;
  }
  double a(0);
  if(x)
    {
      a = angle-xAngle;
      xAngle = angle;
    }
  else if(y)
    {
      a = angle-yAngle;
      yAngle = angle;
    }
  else if(z)
    {
      a = angle-zAngle;
      zAngle = angle;
    }
  double ax(x);
  double ay(y);
  double az(z);
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  GLdouble m[16];
  GLdouble inv[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, m);
  invertMatrix(m, inv);
  double bx(inv[0]*ax + inv[4]*ay + inv[8]*az);
  double by(inv[1]*ax + inv[5]*ay + inv[9]*az);
  double bz(inv[2]*ax + inv[6]*ay + inv[10]*az);

  glTranslatef(mid_point_.x, mid_point_.y, mid_point_.z);
  glRotatef(a,bx,by,bz);
  glTranslatef(-mid_point_.x, -mid_point_.y, -mid_point_.z);

  m_control_->setXangle(xAngle);
  m_control_->setYangle(yAngle);
  m_control_->setZangle(zAngle);
  isShownSurface = false;
  queue_draw();


  /*
  GLfloat m[16];
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX,m);
  glLoadIdentity();
  glTranslatef(Xtrans,Ytrans,-(Near+ViewSize/2.0));
  if(x)
    {
      glRotatef(angle-xAngle,x,y,z);
      xAngle = angle;
      m_control_->setXangle(xAngle);
    }
  if(y)
    {
      glRotatef(angle-yAngle,x,y,z);
      yAngle = angle;
      m_control_->setYangle(yAngle);
    }
  if(z)
    {
      glRotatef(angle-zAngle,x,y,z);
      zAngle = angle;
      m_control_->setZangle(zAngle);
    }
  glTranslatef(-Xtrans,-Ytrans,+(Near+ViewSize/2.0));
  glMultMatrixf(m);
  invalidate();
  */
}

void GLScene::normalizeAngle(double &angle)
{
  while(angle > 180)
    {
      angle -= 360;
    }
  while(angle < -180)
    {
      angle += 360;
    }
}

bool GLScene::get_is_playing() {
  return is_playing_;
}

void GLScene::pause()
{
  is_playing_ = !is_playing_;
  if(is_playing_) {
    timeout_add();
  }
  else {
    timeout_remove();
    invalidate();
  }
}

void GLScene::play()
{
  is_playing_ = true;
  timeout_add();
}

ControlBox::ControlBox(GLScene& anArea, Gtk::Table& aTable) :
  m_area_(anArea),
  m_table(10, 10),
  m_area_table_(aTable),
  theDepthAdj( 0, -200, 130, 5, 0, 0 ),
  theXAdj( 0, -180, 180, 1, 20, 0 ),
  theXLowBoundAdj( 0, 0, 100, 1, 0, 0 ),
  theXUpBoundAdj( 100, 0, 0, 1, 0, 0 ),
  theYAdj( 0, -180, 180, 1, 20, 0 ),
  theYLowBoundAdj( 0, 0, 100, 1, 0, 0 ),
  theYUpBoundAdj( 100, 0, 100, 1, 0, 0 ),
  theZAdj( 0, -180, 180, 1, 20, 0 ),
  theZLowBoundAdj( 0, 0, 100, 1, 0, 0 ),
  theZUpBoundAdj( 100, 0, 100, 1, 0, 0 ),
  theResetTimeButton( "Reset Time" ),
  theResetBoundButton( "Reset" ),
  theResetDepthButton( "Reset" ),
  theResetRotButton( "Reset View" ),
  theCheck3DMolecule( "Show 3D Molecules" ),
  //theCheckFix( "Fix rotation" ),
  theCheckInvertBound( "Invert Bounding" ),
  theCheckShowSurface( "Show Surface" ),
  theCheckShowTime( "Show Time" ),
  theFrameBoundAdj("Bounding"),
  theFrameLatticeAdj("Zoom"),
  theFrameRotAdj( "Rotation" ),
  theFrameScreen("Resolution"),
  theDepthScale( theDepthAdj ),
  theXLowBoundScale( theXLowBoundAdj ),
  theXScale( theXAdj ),
  theXUpBoundScale( theXUpBoundAdj ),
  theYLowBoundScale( theYLowBoundAdj ),
  theYScale( theYAdj ),
  theYUpBoundScale( theYUpBoundAdj ),
  theZLowBoundScale( theZLowBoundAdj ),
  theZScale( theZAdj ),
  theZUpBoundScale( theZUpBoundAdj ),
  theDepthLabel( "Depth" ),
  theHeightLabel( "Height:" ),
  theWidthLabel( "Width:" ),
  theXLabel( "x" ),
  theXLowBoundLabel( "-x" ),
  theXUpBoundLabel( "+x" ),
  theYLabel( "y" ),
  theYLowBoundLabel( "-y" ),
  theYUpBoundLabel( "+y" ),
  theZLabel( "z" ),
  theZLowBoundLabel( "-z" ),
  theZUpBoundLabel( "+z" ),
  theDepthSpin( theDepthAdj, 0, 0  ),
  theXLowBoundSpin( theXLowBoundAdj, 0, 0  ),
  theXSpin( theXAdj, 0, 0  ),
  theXUpBoundSpin( theXUpBoundAdj, 0, 0  ),
  theYLowBoundSpin( theYLowBoundAdj, 0, 0  ),
  theYSpin( theYAdj, 0, 0  ),
  theYUpBoundSpin( theYUpBoundAdj, 0, 0  ),
  theZLowBoundSpin( theZLowBoundAdj, 0, 0  ),
  theZSpin( theZAdj, 0, 0  ),
  theZUpBoundSpin( theZUpBoundAdj, 0, 0  ),
  theRecordButton( "Record Frames" ),
  progress_adj_(0, 1, std::max(1, int(m_area_.get_frame_size()-2)), 1, 0, 0),
  progress_spin_( progress_adj_, 0, 0  ),
  progress_bar_(progress_adj_)
{
  m_area_table_.attach(progress_box_, 0, 1, 1, 2, Gtk::SHRINK | Gtk::FILL,
                       Gtk::SHRINK | Gtk::FILL, 0, 0);
  progress_box_.pack_start(play_button_, false, false, 2);
  play_button_.set_stock_id(Gtk::Stock::MEDIA_PLAY);
  play_button_.signal_clicked().connect( sigc::mem_fun(*this,
                            &ControlBox::play_or_pause) );

  progress_bar_.set_draw_value(false);
  progress_box_.pack_start(progress_bar_);
  progress_spin_.set_width_chars(4);
  progress_spin_.set_has_frame(false);
  frame_cnt_label_.set_text("Frame:");
  progress_box_.pack_start(frame_cnt_label_, false, false, 2);
  progress_box_.pack_start(progress_spin_, false, false, 2);
  progress_adj_.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::progress_changed ) );


  set_border_width(2);
  set_size_request(475, SCREEN_HEIGHT);
  set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
  add(m_rightBox);
  m_table.set_row_spacings(1);
  m_table.set_col_spacings(2);
  m_rightBox.pack_start(m_table, Gtk::PACK_SHRINK);
  //Create a table on the right that holds the control.
  m_rightBox.pack_start(theBoxCtrl, Gtk::PACK_SHRINK); 
  // Create a frame the will have the rotation adjusters
  // and the Fix and Reset buttons.


  /*
  // Create a frame the will have the lattice depth adjuster
  // and background color selector
  theBoxCtrl.pack_start( theFrameLatticeAdj, false, false, 1 );
  theBoxInLattice.set_border_width( 3 );
  theFrameLatticeAdj.add( theBoxInLattice );
  theDepthBox.set_homogeneous( false );
  theBoxInLattice.pack_start( theDepthBox, false, false, 1 );
  */
  theBoxInLattice.pack_start( the3DMoleculeBox, false, false, 1 );
  //theResetDepthButton.connect( 'clicked', resetDepth );
  the3DMoleculeBox.pack_start( theResetDepthButton ); 


  theCheck3DMolecule.signal_toggled().connect( sigc::mem_fun(*this,
                            &ControlBox::on_3DMolecule_toggled) );
  //theCheck3DMolecule.set_active();
  theCheck3DMolecule.set_active(false);
  theBoxCtrl.pack_start( theCheck3DMolecule, false, false, 2 );


  theCheckShowTime.signal_toggled().connect( sigc::mem_fun(*this,
                            &ControlBox::on_showTime_toggled) );
  theCheckShowTime.set_active();
  theBoxCtrl.pack_start( theCheckShowTime, false, false, 2 );

  /*
  theCheckShowSurface.signal_toggled().connect( sigc::mem_fun(*this,
                            &ControlBox::on_showSurface_toggled) );
  //theCheckShowSurface.set_active();
  theCheckShowSurface.set_active(false);
  theBoxCtrl.pack_start( theCheckShowSurface, false, false, 2 );
  */

  theResetRotButton.signal_clicked().connect( sigc::mem_fun(*this,
                            &ControlBox::onResetRotation) );
  theBoxCtrl.pack_start( theResetRotButton, false, false, 2 );

  theResetTimeButton.signal_clicked().connect( sigc::mem_fun(*this,
                            &ControlBox::on_resetTime_clicked) );
  theBoxCtrl.pack_start( theResetTimeButton, false, false, 2 );

  theRecordButton.signal_toggled().connect( sigc::mem_fun(*this,
                            &ControlBox::on_record_toggled) );
  theBoxCtrl.pack_start( theRecordButton, false, false, 2 );



  // screen resolution display
  theBoxCtrl.pack_start( theFrameScreen, false, false, 1 );
  theBoxInScreen.set_border_width( 3 );
  theFrameScreen.add(theBoxInScreen);
  theBoxInScreen.pack_start( theWidthBox, false, false, 1 ); 
  theBoxInScreen.pack_start( theHeightBox, false, false, 1 ); 

  theHeightLabel.set_width_chars( 6 );
  theHeightBox.pack_start( theHeightLabel, false, false, 2 );
  theHeightBox.pack_start( m_height );

  theWidthLabel.set_width_chars( 6 );
  theWidthBox.pack_start( theWidthLabel, false, false, 2 );
  theWidthBox.pack_start( m_width );


  theDepthLabel.set_width_chars( 1 );
  theDepthBox.pack_start( theDepthLabel, false, false, 2 );
  //theDepthAdj.connect( 'value_changed', depthChanged );
  theDepthScale.set_draw_value( false );
  theDepthBox.pack_start( theDepthScale );
  theDepthSpin.set_width_chars( 3 );
  theDepthSpin.set_has_frame( false );
  theDepthBox.pack_start( theDepthSpin, false, false, 2 );

  // Create a frame the will have the lattice boundary adjusters,
  // Fix and Reset buttons.
  theBoxCtrl.pack_start( theFrameBoundAdj, false, false, 1 );
  theBoxInBound.set_border_width( 3 );
  theFrameBoundAdj.add( theBoxInBound );
  theBoxInBound.pack_start( theXUpBoundBox, false, false, 1 ); 
  theBoxInBound.pack_start( theXLowBoundBox, false, false, 1 ); 
  theBoxInBound.pack_start( theYUpBoundBox, false, false, 1 ); 
  theBoxInBound.pack_start( theYLowBoundBox, false, false, 1 ); 
  theBoxInBound.pack_start( theZUpBoundBox, false, false, 1 ); 
  theBoxInBound.pack_start( theZLowBoundBox, false, false, 1 ); 
  theBoxInBound.pack_start( theBoxBoundFixReset, false, false, 1 ); 
  theCheckInvertBound.signal_toggled().connect( sigc::mem_fun(*this,
                            &ControlBox::on_InvertBound_toggled) );
  theBoxBoundFixReset.pack_start( theCheckInvertBound );
  theResetBoundButton.signal_clicked().connect( sigc::mem_fun(*this,
                            &ControlBox::onResetBound) );
  theBoxBoundFixReset.pack_start( theResetBoundButton );

  unsigned int aLayerSize(m_area_.getLayerSize());
  unsigned int aColSize(m_area_.getColSize());
  unsigned int aRowSize(m_area_.getRowSize());

  // x up bound
  theXUpBoundLabel.set_width_chars( 2 );
  theXUpBoundBox.pack_start( theXUpBoundLabel, false, false, 2 );
  theXUpBoundAdj.set_value( aColSize );
  theXUpBoundAdj.set_upper( aColSize );
  theXUpBoundAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::xUpBoundChanged ) );
  theXUpBoundScale.set_draw_value( false );
  theXUpBoundBox.pack_start( theXUpBoundScale );
  theXUpBoundSpin.set_width_chars( 3 );
  theXUpBoundSpin.set_has_frame( false );
  theXUpBoundBox.pack_start( theXUpBoundSpin, false, false, 2 );

  // x low bound
  theXLowBoundLabel.set_width_chars( 2 );
  theXLowBoundBox.pack_start( theXLowBoundLabel, false, false, 2 );
  theXLowBoundAdj.set_upper( aColSize );
  theXLowBoundAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::xLowBoundChanged ) );
  theXLowBoundScale.set_draw_value( false );
  theXLowBoundBox.pack_start( theXLowBoundScale );
  theXLowBoundSpin.set_width_chars( 3 );
  theXLowBoundSpin.set_has_frame( false );
  theXLowBoundBox.pack_start( theXLowBoundSpin, false, false, 2 );

  // y up bound
  theYUpBoundLabel.set_width_chars( 2 );
  theYUpBoundBox.pack_start( theYUpBoundLabel, false, false, 2 );
  theYUpBoundAdj.set_value( aLayerSize );
  theYUpBoundAdj.set_upper( aLayerSize );
  theYUpBoundAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::yUpBoundChanged ) );
  theYUpBoundScale.set_draw_value( false );
  theYUpBoundBox.pack_start( theYUpBoundScale );
  theYUpBoundSpin.set_width_chars( 3 );
  theYUpBoundSpin.set_has_frame( false );
  theYUpBoundBox.pack_start( theYUpBoundSpin, false, false, 2 );

  // y low bound
  theYLowBoundLabel.set_width_chars( 2 );
  theYLowBoundBox.pack_start( theYLowBoundLabel, false, false, 2 );
  theYLowBoundAdj.set_upper( aLayerSize );
  theYLowBoundAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::yLowBoundChanged ) );
  theYLowBoundScale.set_draw_value( false );
  theYLowBoundBox.pack_start( theYLowBoundScale );
  theYLowBoundSpin.set_width_chars( 3 );
  theYLowBoundSpin.set_has_frame( false );
  theYLowBoundBox.pack_start( theYLowBoundSpin, false, false, 2 );

  // z up bound
  theZUpBoundLabel.set_width_chars( 2 );
  theZUpBoundBox.pack_start( theZUpBoundLabel, false, false, 2 );
  theZUpBoundAdj.set_value( aRowSize );
  theZUpBoundAdj.set_upper( aRowSize );
  theZUpBoundAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::zUpBoundChanged ) );
  theZUpBoundScale.set_draw_value( false );
  theZUpBoundBox.pack_start( theZUpBoundScale );
  theZUpBoundSpin.set_width_chars( 3 );
  theZUpBoundSpin.set_has_frame( false );
  theZUpBoundBox.pack_start( theZUpBoundSpin, false, false, 2 );

  // z low bound
  theZLowBoundLabel.set_width_chars( 2 );
  theZLowBoundBox.pack_start( theZLowBoundLabel, false, false, 2 );
  theZLowBoundAdj.set_upper( aRowSize );
  theZLowBoundAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::zLowBoundChanged ) );
  theZLowBoundScale.set_draw_value( false );
  theZLowBoundBox.pack_start( theZLowBoundScale );
  theZLowBoundSpin.set_width_chars( 3 );
  theZLowBoundSpin.set_has_frame( false );
  theZLowBoundBox.pack_start( theZLowBoundSpin, false, false, 2 );

  // rotation adjusters
  theBoxCtrl.pack_start( theFrameRotAdj, false, false, 1 );
  theBoxInFrame.set_border_width( 3 );
  theFrameRotAdj.add( theBoxInFrame );
  theXBox.set_homogeneous( false );
  theBoxInFrame.pack_start( theXBox, false, false, 1 );
  theBoxInFrame.pack_start( theYBox, false, false, 1 );
  theBoxInFrame.pack_start( theZBox, false, false, 1 );
  theBoxInFrame.pack_start( theBoxRotFixReset, false, false, 1 );
  //theCheckFix.connect( 'toggled', fixRotToggled );
  //theBoxRotFixReset.pack_start( theCheckFix );

  // X
  theXLabel.set_width_chars( 1 );
  theXBox.pack_start( theXLabel, false, false, 2 ); 
  theXScale.set_draw_value( false );
  theXBox.pack_start( theXScale );
  theXSpin.set_width_chars( 3 );
  theXSpin.set_wrap( true );
  theXSpin.set_has_frame( false );
  theXBox.pack_start( theXSpin, false, false, 2 );
  theXAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::xRotateChanged ) );

  // Y
  theYLabel.set_width_chars( 1 );
  theYBox.pack_start( theYLabel, false, false, 2 ); 
  theYScale.set_draw_value( false );
  theYBox.pack_start( theYScale );
  theYSpin.set_width_chars( 3 );
  theYSpin.set_wrap( true );
  theYSpin.set_has_frame( false );
  theYBox.pack_start( theYSpin, false, false, 2 );
  theYAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::yRotateChanged ) );

  // Z
  theZLabel.set_width_chars( 1 );
  theZBox.pack_start( theZLabel, false, false, 2 ); 
  theZScale.set_draw_value( false );
  theZBox.pack_start( theZScale );
  theZSpin.set_width_chars( 3 );
  theZSpin.set_wrap( true );
  theZSpin.set_has_frame( false );
  theZBox.pack_start( theZSpin, false, false, 2 );
  theZAdj.signal_value_changed().connect( sigc::mem_fun(*this, 
                           &ControlBox::zRotateChanged ) );

  // Time and step counts
  m_sizeGroup = Gtk::SizeGroup::create(Gtk::SIZE_GROUP_HORIZONTAL);
  m_table.attach(m_timeBox, 0, 1, 0, 1, Gtk::FILL,
                 Gtk::SHRINK | Gtk::FILL, 0, 0 );

  m_bgColor.set_text("Background Color");
  Gtk::EventBox* eventBox=Gtk::manage(new Gtk::EventBox);
  eventBox->set_events(Gdk::BUTTON_PRESS_MASK);
  eventBox->signal_button_press_event().connect(
    sigc::mem_fun(*this, &ControlBox::on_background_clicked));
  eventBox->add(m_bgColor);
  eventBox->set_tooltip_text("Click to change background color");
  m_table.attach(*eventBox, 0, 1, 1, 2, Gtk::FILL,
                 Gtk::SHRINK | Gtk::FILL, 0, 0 );
  theBgColor.set_rgb(0, 0, 0);

  m_timeLabel.set_text("Time:");
  m_sizeGroup->add_widget(m_timeLabel);
  m_timeBox.pack_start(m_timeLabel, Gtk::PACK_SHRINK);
  m_timeBox.pack_start(m_time, Gtk::PACK_SHRINK);
  theSpeciesSize = m_area_.getSpeciesSize();
  theButtonList = new Gtk::CheckButton*[theSpeciesSize]; 
  theLabelList = new Gtk::Label*[theSpeciesSize]; 
  for(unsigned int i(0); i!=theSpeciesSize; ++i)
    {
      theButtonList[i] = Gtk::manage(new Gtk::CheckButton());
      theLabelList[i] = Gtk::manage(new Gtk::Label);
      Gtk::HBox* aHBox = Gtk::manage(new Gtk::HBox);
      Gtk::EventBox* m_EventBox=Gtk::manage(new Gtk::EventBox);
      m_EventBox->set_events(Gdk::BUTTON_RELEASE_MASK);
      theLabelList[i]->set_text(m_area_.getSpeciesName(i));
      Color clr(m_area_.getSpeciesColor(i));
      Gdk::Color aColor;
      aColor.set_rgb(int(clr.r*65535), int(clr.g*65535), int(clr.b*65535));
      theLabelList[i]->modify_fg(Gtk::STATE_NORMAL, aColor);
      theLabelList[i]->modify_fg(Gtk::STATE_ACTIVE, aColor);
      theLabelList[i]->modify_fg(Gtk::STATE_PRELIGHT, aColor);
      theLabelList[i]->modify_fg(Gtk::STATE_SELECTED, aColor);
      theButtonList[i]->signal_toggled().connect(sigc::bind( 
              sigc::mem_fun(*this, &ControlBox::on_checkbutton_toggled), i ) );
      m_EventBox->signal_button_release_event().connect(sigc::bind( 
              sigc::mem_fun(*this, &ControlBox::on_checkbutton_clicked), i ) );
      m_EventBox->add(*theLabelList[i]);
      m_EventBox->set_tooltip_text("Right click to change species color");
      theButtonList[i]->set_active(m_area_.getSpeciesVisibility(i));
      aHBox->pack_start(*theButtonList[i], false, false, 2);
      aHBox->pack_start(*m_EventBox, false, false, 2);
      m_table.attach(*aHBox, 0, 1, i+2, i+3, Gtk::FILL,
                     Gtk::SHRINK | Gtk::FILL, 0, 0 );
    }
  std::cout << "theSpeciesSize:" << theSpeciesSize << std::endl;
}

void ControlBox::play_or_pause() {
  if(m_area_.get_is_playing()) {
    pause();
  }
  else {
    play();
  }
}

void ControlBox::play() {
  m_area_.play();
  play_button_.set_stock_id(Gtk::Stock::MEDIA_PAUSE);
}

void ControlBox::pause() {
  m_area_.pause();
  play_button_.set_stock_id(Gtk::Stock::MEDIA_PLAY);
}


void
ControlBox::on_checkbutton_toggled(unsigned int id)
{
  m_area_.setSpeciesVisibility(id, theButtonList[id]->get_active());
}

bool ControlBox::on_background_clicked(GdkEventButton* event)
{
  if(event)
    {
      Gtk::ColorSelectionDialog dlg("Select a color"); 
      Gtk::ColorSelection* colorSel(dlg.get_colorsel());
      colorSel->set_current_color(theBgColor);
      colorSel->signal_color_changed().connect(sigc::bind( 
        sigc::mem_fun(*this, &ControlBox::update_background_color), colorSel));
      if(dlg.run() == Gtk::RESPONSE_CANCEL)
        {
          Color clr;
          clr.r = theBgColor.get_red()/65535.0;
          clr.g = theBgColor.get_green()/65535.0;
          clr.b = theBgColor.get_blue()/65535.0;
          m_area_.setBackgroundColor(clr);
          m_bgColor.modify_fg(Gtk::STATE_NORMAL, theBgColor);
          m_bgColor.modify_fg(Gtk::STATE_ACTIVE, theBgColor);
          m_bgColor.modify_fg(Gtk::STATE_PRELIGHT, theBgColor);
          m_bgColor.modify_fg(Gtk::STATE_SELECTED, theBgColor);
        }
      else
        {
          theBgColor = colorSel->get_current_color();
          update_background_color(colorSel);
        }
      return true;
    }
  return false;
}

bool ControlBox::on_checkbutton_clicked(GdkEventButton* event, unsigned int id)
{
  if(event->type == GDK_BUTTON_RELEASE && event->button == 3)
    {
      Color clr(m_area_.getSpeciesColor(id));
      Gdk::Color aColor;
      aColor.set_rgb(int(clr.r*65535), int(clr.g*65535), int(clr.b*65535));
      Gtk::ColorSelectionDialog dlg("Select a color"); 
      Gtk::ColorSelection* colorSel(dlg.get_colorsel());
      colorSel->set_current_color(aColor);
      colorSel->signal_color_changed().connect(sigc::bind( 
        sigc::mem_fun(*this, &ControlBox::update_species_color), id, colorSel));
      if(dlg.run() == Gtk::RESPONSE_CANCEL)
        {
          m_area_.setSpeciesColor(id, clr);
          theLabelList[id]->modify_fg(Gtk::STATE_NORMAL, aColor);
          theLabelList[id]->modify_fg(Gtk::STATE_ACTIVE, aColor);
          theLabelList[id]->modify_fg(Gtk::STATE_PRELIGHT, aColor);
          theLabelList[id]->modify_fg(Gtk::STATE_SELECTED, aColor);
        }
      else
        {
          update_species_color(id, colorSel);
        }
      return true;
    }
  else
    {
      if(theButtonList[id]->get_active())
        {
          theButtonList[id]->set_active(false);
        }
      else
        {
          theButtonList[id]->set_active(true);
        }
    }
  return false;
}

void ControlBox::update_species_color(unsigned int id,
                                      Gtk::ColorSelection* colorSel)
{
  Gdk::Color aColor(colorSel->get_current_color());
  Color clr;
  clr.r = aColor.get_red()/65535.0;
  clr.g = aColor.get_green()/65535.0;
  clr.b = aColor.get_blue()/65535.0;
  m_area_.setSpeciesColor(id, clr);
  theLabelList[id]->modify_fg(Gtk::STATE_NORMAL, aColor);
  theLabelList[id]->modify_fg(Gtk::STATE_ACTIVE, aColor);
  theLabelList[id]->modify_fg(Gtk::STATE_PRELIGHT, aColor);
  theLabelList[id]->modify_fg(Gtk::STATE_SELECTED, aColor);
}

void ControlBox::update_background_color(Gtk::ColorSelection* colorSel)
{
  Gdk::Color aColor(colorSel->get_current_color());
  Color clr;
  clr.r = aColor.get_red()/65535.0;
  clr.g = aColor.get_green()/65535.0;
  clr.b = aColor.get_blue()/65535.0;
  m_area_.setBackgroundColor(clr);
  m_bgColor.modify_fg(Gtk::STATE_NORMAL, aColor);
  m_bgColor.modify_fg(Gtk::STATE_ACTIVE, aColor);
  m_bgColor.modify_fg(Gtk::STATE_PRELIGHT, aColor);
  m_bgColor.modify_fg(Gtk::STATE_SELECTED, aColor);
}

void ControlBox::on_3DMolecule_toggled()
{
  m_area_.set3DMolecule(theCheck3DMolecule.get_active());
}

void ControlBox::on_InvertBound_toggled()
{
  m_area_.setInvertBound(theCheckInvertBound.get_active());
}

void ControlBox::on_showTime_toggled()
{
  m_area_.setShowTime(theCheckShowTime.get_active());
}

void ControlBox::on_showSurface_toggled()
{
  m_area_.setShowSurface(theCheckShowSurface.get_active());
}

void ControlBox::on_resetTime_clicked()
{
  m_area_.resetTime();
}

void ControlBox::onResetRotation()
{
  m_area_.resetView();
}

void ControlBox::onResetBound()
{
  unsigned int aLayerSize(m_area_.getLayerSize());
  unsigned int aColSize(m_area_.getColSize());
  unsigned int aRowSize(m_area_.getRowSize());
  theXUpBoundAdj.set_value( 0 );
  theYUpBoundAdj.set_value( aLayerSize );
  theZUpBoundAdj.set_value( aRowSize );
  theCheckInvertBound.set_active(false);
  m_area_.resetBound();
}

void ControlBox::on_record_toggled()
{
  m_area_.setRecord(theRecordButton.get_active());
}

void ControlBox::setXangle(double angle)
{
  theXAdj.set_value(angle);
}

void ControlBox::setYangle(double angle)
{
  theYAdj.set_value(angle);
}

void ControlBox::setZangle(double angle)
{
  theZAdj.set_value(angle);
}

void ControlBox::xRotateChanged()
{
  m_area_.rotateMidAxisAbs(theXAdj.get_value(), 1, 0, 0);
}

void ControlBox::yRotateChanged()
{
  m_area_.rotateMidAxisAbs(theYAdj.get_value(), 0, 1, 0);
}

void ControlBox::zRotateChanged()
{
  m_area_.rotateMidAxisAbs(theZAdj.get_value(), 0, 0, 1);
}

void ControlBox::resizeScreen(unsigned aWidth, unsigned aHeight)
{
  std::stringstream w, h;
  w << aWidth;
  h << aHeight;
  m_width.set_text(w.str().c_str());
  m_height.set_text(h.str().c_str());
}

void ControlBox::step()
{
  play_button_.set_stock_id(Gtk::Stock::MEDIA_PLAY);
  m_area_.step();
}

void ControlBox::progress_changed()
{
  if(m_area_.get_frame_cnt() != progress_adj_.get_value()) {
    play_button_.set_stock_id(Gtk::Stock::MEDIA_PLAY);
    m_area_.set_frame_cnt(progress_adj_.get_value());
  }
}

void
ControlBox::xUpBoundChanged()
{
  m_area_.setXUpBound((unsigned int)theXUpBoundAdj.get_value());
}

void
ControlBox::xLowBoundChanged()
{
  m_area_.setXLowBound((unsigned int)theXLowBoundAdj.get_value());
}

void
ControlBox::yUpBoundChanged()
{
  m_area_.setYUpBound((unsigned int)theYUpBoundAdj.get_value());
}

void
ControlBox::yLowBoundChanged()
{
  m_area_.setYLowBound((unsigned int)theYLowBoundAdj.get_value());
}

void
ControlBox::zUpBoundChanged()
{
  m_area_.setZUpBound((unsigned int)theZUpBoundAdj.get_value());
}

void
ControlBox::zLowBoundChanged()
{
  m_area_.setZLowBound((unsigned int)theZLowBoundAdj.get_value());
}

void
ControlBox::set_frame_cnt(int frame_cnt)
{
  if(frame_cnt > progress_adj_.get_upper()) {
    progress_adj_.set_upper(frame_cnt);
  }
  if(frame_cnt != progress_adj_.get_value()) { 
    progress_adj_.set_value(frame_cnt);
  }
}

void
ControlBox::setTime(double time)
{
  std::stringstream s;
  s << time;
  m_time.set_text(s.str().c_str());
}

ControlBox::~ControlBox()
{
}


Rulers::Rulers(const Glib::RefPtr<const Gdk::GL::Config>& config,
               const char* aFileName) :
  m_area_(config, aFileName),
  m_hbox(),
  m_table(3, 1, false),
  m_control_(m_area_, m_table),
  isRecord(false)
{
  m_area_.setControlBox(&m_control_);
  set_title("Spatiocyte Visualizer");
  //set_reallocate_redraws(true);
  set_border_width(10);

  add(m_hbox);
  m_hbox.pack1(m_table, Gtk::PACK_EXPAND_WIDGET, 5);
  m_hbox.pack2(m_control_, Gtk::PACK_SHRINK, 5);
  m_table.attach(m_area_, 0,1,0,1,
		 Gtk::EXPAND | Gtk::FILL , Gtk::EXPAND | Gtk::FILL, 0, 0);
  show_all_children();
}

void GLScene::on_size_allocate(Gtk::Allocation& allocation) {
  m_control_->resizeScreen(allocation.get_width(), allocation.get_height());
  Gtk::GL::DrawingArea::on_size_allocate(allocation);
}

bool Rulers::on_key_press_event(GdkEventKey* event)
{
  switch (event->keyval)
    {
    case GDK_x:
      m_area_.rotate(1,1,0,0);
      break;
    case GDK_X:
      m_area_.rotate(-1,1,0,0);
      break;
    case GDK_y:
      m_area_.rotate(1,0,1,0);
      break;
    case GDK_Y:
      m_area_.rotate(-1,0,1,0);
      break;
    case GDK_Home:
      m_area_.resetView();
      break;
    case GDK_Pause:
      m_control_.pause();
      break;
    case GDK_p:
      m_control_.pause();
      break;
    case GDK_P:
      m_control_.pause();
      break;
    case GDK_Return:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.set_is_forward(false);
          m_control_.step();
        }
      else
        {
          m_area_.set_is_forward(true);
          m_control_.step();
        }
      break;
    case GDK_space:
      m_control_.play_or_pause();
      break;
    case GDK_F:
      m_area_.translate(0,0,-1);
      break;
    case GDK_f:
      m_area_.zoomIn();
      break;
    case GDK_B:
      m_area_.translate(0,0,1);
      break;
    case GDK_b:
      m_area_.zoomOut();
      break;
    case GDK_Page_Up:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.translate(0,0,-1);
        }
      else
        {
          m_area_.zoomIn();
        }
      break;
    case GDK_Page_Down:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.translate(0,0,1);
        }
      else
        {
          m_area_.zoomOut();
        }
      break;
    case GDK_0:
      if(event->state&Gdk::CONTROL_MASK)
        {
          m_area_.resetView();
        }
      break;
    case GDK_equal:
      if(event->state&Gdk::CONTROL_MASK)
        {
          m_area_.zoomIn();
        }
      break;
    case GDK_plus:
      if(event->state&Gdk::CONTROL_MASK)
        {
          m_area_.zoomIn();
        }
      break;
    case GDK_minus:
      if(event->state&Gdk::CONTROL_MASK)
        {
          m_area_.zoomOut();
        }
      break;
    case GDK_Down:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.translate(0,-1,0);
        }
      else if (event->state&Gdk::CONTROL_MASK)
        {
          m_area_.rotateMidAxis(1,1,0,0);
        }
      else
        {
          m_area_.set_is_forward(false);
          m_control_.step();
        }
      break;
    case GDK_Up:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.translate(0,1,0);
        }
      else if (event->state&Gdk::CONTROL_MASK)
        {
          m_area_.rotateMidAxis(-1,1,0,0);
        }
      else
        {
          m_area_.set_is_forward(true);
          m_control_.step();
        }
      break;
    case GDK_Right:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.translate(1,0,0);
        }
      else if (event->state&Gdk::CONTROL_MASK)
        {
          m_area_.rotateMidAxis(1,0,1,0);
        }
      else
        {
          m_area_.set_is_forward(true);
          m_control_.play();
        }
      break;
    case GDK_Left:
      if(event->state&Gdk::SHIFT_MASK)
        {
          m_area_.translate(-1,0,0);
        }
      else if (event->state&Gdk::CONTROL_MASK)
        {
          m_area_.rotateMidAxis(-1,0,1,0);
        }
      else
        {
          m_area_.set_is_forward(false);
          m_control_.play();
        }
      break;
    case GDK_z:
      m_area_.rotateMidAxis(-1,0,0,1);
      break;
    case GDK_Z:
      m_area_.rotateMidAxis(1,0,0,1);
      break;
    case GDK_l:
      m_area_.rotate(1,0,0,1);
      break;
    case GDK_r:
      m_area_.rotate(-1,0,0,1);
      break;
    case GDK_s:
      std::cout << "saving frame" << std::endl;
      m_area_.writePng();
      break;
    case GDK_S:
      if(!isRecord)
        {
          isRecord = true;
          std::cout << "Started saving frames" << std::endl; 
        }
      else
        {
          isRecord = false;
          std::cout << "Stopped saving frames" << std::endl; 
        }
      m_area_.setRecord(isRecord);
      break;
    default:
      return true;
    }
  return true;
}

void printUsage( const char* aProgramName )
{
  std::cerr << "usage:" << std::endl;
  std::cerr <<  aProgramName <<
    " <fileBaseName> (for Little Endian binary files)" << std::endl;
  std::cerr <<  aProgramName <<
    " -b <fileBaseName> (for Big Endian binary files)"
    << std::endl << std::flush;
}

void printNotEndian( const char* anEndian, const char* aBaseName)
{
  std::cerr << "file is " << aBaseName << 
    "\nbut I can't open " << aBaseName << " file.\n";
  std::cerr << "Could it be that " << aBaseName << " is not";
  std::cerr << " a " << anEndian << " Endian binary file?\n";
}

int main(int argc, char** argv)
{
  std::string aBaseName;
  if(argc == 1)
    {
      aBaseName = "VisualLog.dat";
    }
  else if(argc == 2)
    {
      aBaseName = argv[1];
    }
  else if(argc == 3)
    {
      if( argv[1][1] == 'b' )
        {
        }
      else
        {
          std::cout << "Unknown option " << argv[1] << " ignored\n";
          printUsage( argv[0] );
          std::exit(1);
        }
    }
  else
    {
      printUsage( argv[0] );
      std::exit(1);
    }
  std::ostringstream aFileName;
  aFileName << aBaseName << std::ends;
  std::ifstream aParentFile( aFileName.str().c_str(), std::ios::binary );
  if ( !aParentFile.is_open() )
    {
      std::cerr << "Could not open file: " << aBaseName <<  
        std::endl;
      printUsage( argv[0] );
      std::exit(1);
    }
  else
    {
      std::ostringstream aFileName;
      aFileName << aBaseName << std::ends;
      std::ifstream aFile( aFileName.str().c_str(), std::ios::binary );
      if( !aFile.is_open() )
        {
          printNotEndian( "Little", aBaseName.c_str());
          printUsage( argv[0] );
          std::exit(1);
        }
    }
  Gtk::Main anApp(argc, argv);
  Gtk::GL::init(argc, argv);
  Glib::RefPtr<Gdk::GL::Config> aGLConfig;
  aGLConfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB |
                                      Gdk::GL::MODE_DEPTH |
                                      Gdk::GL::MODE_DOUBLE);
  if (!aGLConfig)
    {
      std::cerr << "*** Cannot find the double-buffered visual.\n"
                << "*** Trying single-buffered visual.\n";
      aGLConfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB   |
                                         Gdk::GL::MODE_DEPTH);
      if (!aGLConfig)
        {
          std::cerr << "*** Cannot find any OpenGL-capable visual.\n";
          std::exit(1);
        }
    }

  Rulers aRuler(aGLConfig, aBaseName.c_str());
  Gtk::Main::run(aRuler);
  return 0;
}




