
/* wiimote head tracking for linux
   Originally by Steven Thomas Snyder, stsnyder@ucla.edu
   (C) copyright 2008, Steven Snyder, All Rights Reserved
   (C) copyright 2011, Marcus Bannerman, All Rights Reserved

  LICENSING INFORMATION:
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef COIL_wiimote
#include "wiiheadtracking.hpp"
#include <cmath>
#include <stdexcept>

//Hidden namespace for constants
namespace {
  //The wiimote is factory calibrated such that the angle per pixel
  //for both dimensions is equal. Thus, only the calibration data is
  //needed for the X dimension. Although the FOV is stated to be less
  //than 45 degree's on other websites, this value yeilds accurate
  //(+-5mm) measurements on my wiimote
  const double WiiFOVX = (45.0f / 180.0) * M_PI; 
  //The angle corresponding to a pixel (the camera only "sees" angles, not distances)
  const double anglePerPixel =  WiiFOVX / double(CWIID_IR_X_MAX);  
  //Distance between the two IR sources being tracked (in cm)
  const double IRPointSeparation = 16.5;
  const double simlength = 50; //The length of a simulation unit length in cm (the zoom factor)
  //The screen dimensions in cm
  const double ScreenXlength = 41.1;
  const double ScreenYlength = 30.9;
}

TrackWiimote::TrackWiimote(bool wiimoteAboveScreen):
  m_wiimote(NULL),
  calibrate_request(false),
  v_angle(0),
  _wiimoteAboveScreen(wiimoteAboveScreen)
{
  eye_pos[0] = eye_pos[1] = 0;
  eye_pos[2] = 40;
  for (size_t i(0); i < 4; ++i)
    for (size_t j(0); j < 2; ++j)
      ir_positions[i][j] = 0;

  for (size_t i(0); i < 4; ++i)
    ir_sizes[i]= 4;
}

bool TrackWiimote::updateState()
{
  if (!m_wiimote) return false;

  if (cwiid_get_state(m_wiimote, &m_state))
    return false; //Failed to obtain data from the remote

  return updateIRPositions() == 2;//Return success if we have two IR sources
}

bool TrackWiimote::connect(bdaddr_t* bt_address)
{
  if (m_wiimote) return true;

  m_wiimote = cwiid_open(bt_address, CWIID_FLAG_CONTINUOUS|CWIID_FLAG_NONBLOCK);
  
  if (!m_wiimote) return false;//Couldn't connect

  if (cwiid_command(m_wiimote, CWIID_CMD_RPT_MODE, 
		    CWIID_RPT_IR | CWIID_RPT_ACC | CWIID_RPT_BTN
		    ) != 0)
    throw std::runtime_error("Failed to enable Wii functions.");

  if (cwiid_command(m_wiimote, CWIID_CMD_LED,
		    CWIID_LED1_ON) != 0)
    throw std::runtime_error("Failed to enable Wii functions.");

  return true;
}

//! Downloads the ir positions and returns how many were recorded
size_t TrackWiimote::updateIRPositions()
{
  size_t points = 0;

  for (size_t i(0); i < 4; ++i)
    if (m_state.ir_src[i].valid)
      {
	ir_positions[i][0] = m_state.ir_src[i].pos[CWIID_X];
	ir_positions[i][1] = m_state.ir_src[i].pos[CWIID_Y];
	
	if (m_state.ir_src[i].size != -1)
	  ir_sizes[i] = m_state.ir_src[i].size + 1;
	
	++points;
      }
  
  if (points == 2) 
    {
      if (calibrate_request)
	{       
	  calibrate();
	  calibrate_request = false;
	}
      
      updateHeadPos();
    }
  
  return points;
}

/* update the camera position using ir points camerapt_1 and camerapt_2 */
void TrackWiimote::updateHeadPos()
{
  //The positions of the IR points, in angles from the centre of the
  //wiimote's view
  double x1 = (ir_positions[0][0] - CWIID_IR_X_MAX / 2) * anglePerPixel;
  double y1 = (ir_positions[0][1] - CWIID_IR_Y_MAX / 2) * anglePerPixel;
  double x2 = (ir_positions[1][0] - CWIID_IR_X_MAX / 2) * anglePerPixel;
  double y2 = (ir_positions[1][1] - CWIID_IR_Y_MAX / 2) * anglePerPixel;

  //The seperation angles of the two IR points
  double dx = x1 - x2;
  double dy = y1 - y2;

  //The absolute angle between the two points
  double pointsAngle = std::sqrt(dx * dx + dy * dy);

  //The distance to the points from the camera (Z coordinate). All
  //distances in eye_pos are in whatever units IRPointSeparation is
  //in.
  eye_pos[2] = 0.5 * IRPointSeparation / std::tan(pointsAngle / 2);

  //The "position" angle of the IR points
  double Xangle = (x1 + x2) / 2;
  //The Y angle is corrected for the angle of the remote (it should be
  //parallel to the normal of the screen, but that rarely works thanks
  //to the remotes design). v_angle is set in the calibration routines
  double Yangle = (y1 + y2) / 2 + v_angle;

  //Assuming that we worked out the distance to the points, and the
  //distance is a radius we can work out the distance from the centre
  eye_pos[0] = -eye_pos[2] * std::sin(Xangle);

  //We try to correct for the remote being off the center of the
  //screen. We simply add the vector to the camera from the screen
  //(only the y component is done here).
  eye_pos[1] = eye_pos[2] * std::sin(Yangle) + ((_wiimoteAboveScreen) ? 0.5f : -0.5f) * ScreenYlength;
}

void TrackWiimote::glPerspective(const magnet::GL::ViewPort& vp, const Vector offset)
{
  //Build a matrix to rotate from camera to world
  Matrix Transformation 
    = Rodrigues(Vector(0, -vp.getPan() * M_PI/180, 0))
    * Rodrigues(Vector(-vp.getTilt() * M_PI / 180.0, 0, 0));
  
  Vector movement = Transformation * (Vector(eye_pos[0], eye_pos[1], eye_pos[2]) + offset ) ;

  glMatrixMode(GL_MODELVIEW);
  glTranslatef(-movement[0] / simlength,
	       -movement[1] / simlength,
	       -movement[2] / simlength);
  
  //We have moved the camera to the location of the head in sim
  //space. Now we must create a viewing frustrum which, in real
  //space, cuts through the image on the screen. The trick is to
  //take the real world relative coordinates of the screen and
  //head transform them to simulation units.
  //
  //This allows us to calculate the left, right, bottom and top of
  //the frustrum as if the near plane of the frustrum was at the
  //screens location.
  //
  //Finally, all length scales are multiplied by
  //(vp._zNearDist/eye_pos[2]).
  //
  //This is to allow the frustrum's near plane to be placed
  //somewhere other than the screen (this factor places it at
  //_zNearDist)!
  //
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum((-0.5f * ScreenXlength - eye_pos[0]) * vp.getZNear() / eye_pos[2],// left
	    (+0.5f * ScreenXlength - eye_pos[0]) * vp.getZNear() / eye_pos[2],// right
	    (-0.5f * ScreenYlength - eye_pos[1]) * vp.getZNear() / eye_pos[2],// bottom 
	    (+0.5f * ScreenYlength - eye_pos[1]) * vp.getZNear() / eye_pos[2],// top
	    vp.getZNear(),//Near distance
	    vp.getZFar()//Far distance
	    );
  
  glMatrixMode(GL_MODELVIEW);
}

void TrackWiimote::calibrate()
{
  //The person is in the center of the screen, so we can determine the
  //tilt of the camera
  v_angle = std::acos(0.5 *  ScreenYlength / eye_pos[2]) - M_PI / 2;
  if (!_wiimoteAboveScreen)
    v_angle = -v_angle;
}
#endif
