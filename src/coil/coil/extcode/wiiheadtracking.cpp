
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
#include "wiiheadtracking.hpp"
#include <cmath>
#include <stdexcept>

//Hidden namespace for constants
namespace {
  //The wiimote is factory calibrated such that the angle per pixel
  //for both dimensions is equal. Thus, only the calibration data is
  //needed for the X dimension. Although the FOV is stated to be less
  //than 45 degree's on other websites, this value yeilds accurate
  //(+-5mm) measurements on my wiimote for all distances
  const double WiiFOVX = (45.0f / 180.0) * M_PI; 
  //The angle corresponding to a pixel (the camera only sees angles, not distances)
  const double anglePerPixel =  WiiFOVX / double(CWIID_IR_X_MAX);  
}

TrackWiimote::TrackWiimote():
  m_wiimote(NULL),
  calibrate_request(false)
{
  eye_pos[0] = eye_pos[1] = eye_pos[2] = 0;
  
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

  m_wiimote = cwiid_connect(bt_address, CWIID_FLAG_CONTINUOUS|CWIID_FLAG_NONBLOCK);
  
  if (!m_wiimote) return false;//Couldn't connect

  if (cwiid_command(m_wiimote, CWIID_CMD_RPT_MODE, 
		    CWIID_RPT_IR | CWIID_RPT_ACC | CWIID_RPT_BTN
		    ) != 0)
    throw std::runtime_error("Failed to enable Wii functions.");

  if (cwiid_command(m_wiimote, CWIID_CMD_LED_MODE,
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

  //The distance to the points from the camera (Z coordinate). As with
  //all values in eye_pos, the units are such that the distance
  //between the IR sources is equal to 1. Just multiply by the 
  eye_pos[2] = 0.5 / std::tan(pointsAngle / 2);

  //The "position" angle of the IR points
  double Xangle = (x1 + x2) / 2;
  double Yangle = (y1 + y2) / 2;

  //Assuming that we worked out the distance to the points, and the
  //distance is a radius we can work out the distance from the centre
  eye_pos[0] = eye_pos[2] * std::sin(Xangle);
  eye_pos[1] = eye_pos[2] * std::sin(Yangle);
}

void TrackWiimote::calibrate()
{}
