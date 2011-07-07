
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
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

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
  const double IRPointSeparation = 15.3;

  void cwiid_err_hidden(struct wiimote*, const char*, va_list) {}
}

TrackWiimote::TrackWiimote():
  m_wiimote(NULL),
  eye_pos(0, 0, 50),
  v_angle(0),
  _batteryLevel(0)
{
  for (size_t i(0); i < CWIID_IR_SRC_COUNT; ++i)
    _ir_data.src[i].valid = false;

#ifndef MAGNET_DEBUG
  cwiid_set_err(&cwiid_err_hidden);
#endif
}

TrackWiimote::~TrackWiimote()
{
  if (m_wiimote) 
    cwiid_close(m_wiimote);
}

bool TrackWiimote::connect(bdaddr_t* bt_address)
{
  if (m_wiimote) return true;

  m_wiimote = cwiid_open(bt_address, CWIID_FLAG_MESG_IFC);
  
  if (!m_wiimote) return false;//Couldn't connect

  if (cwiid_set_mesg_callback(m_wiimote, &TrackWiimote::cwiid_callback)) 
    {
      if (cwiid_close(m_wiimote))
	M_throw() << "Failed to close the Wiimote, after failing to set the callback.";

      M_throw() << "Failed to set the Wiimote callback.";
    }

  if (cwiid_set_rpt_mode(m_wiimote, CWIID_RPT_IR | CWIID_RPT_BTN | CWIID_RPT_STATUS) != 0)
    throw std::runtime_error("Failed to enable Wii functions.");

  if (cwiid_command(m_wiimote, CWIID_CMD_LED,
		    CWIID_LED1_ON) != 0)
    throw std::runtime_error("Failed to enable Wii functions.");
  
  cwiid_request_status(m_wiimote);
  
  return true;
}

//! Downloads the ir positions and returns how many were recorded
size_t TrackWiimote::updateIRPositions()
{
  _valid_ir_points = 0;

  for (size_t i(0); i < 4; ++i)
    if (getIRState(i).valid)
      ++_valid_ir_points;

  if (_valid_ir_points == 2)
    updateHeadPos();
  
  return _valid_ir_points;
}

/* update the camera position using ir points camerapt_1 and camerapt_2 */
void TrackWiimote::updateHeadPos()
{
  if (_valid_ir_points != 2)
    M_throw() << "Need two IR points to perform the head position calculations.";

  //We have to find the first two valid IR positions and store them here
  uint16_t ir_positions[2][2];
  {
    size_t done(0);
    for (size_t i(0); (i < 4) && (done < 2); ++i)
      if (getIRState(i).valid)
	{
	  ir_positions[done][0] = getIRState(i).pos[CWIID_X];
	  ir_positions[done][1] = getIRState(i).pos[CWIID_Y];
	  ++done;
	}

    if (done != 2)
      M_throw() << "_valid_ir_points does not match valid number of points!";
  }

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
  eye_pos[1] = eye_pos[2] * std::sin(Yangle);
}

void TrackWiimote::calibrate()
{
  if (_valid_ir_points != 2) return;

  //The person is in the center of the screen, so we can determine the
  //tilt of the camera
  //v_angle = std::acos(0.5 *  ScreenYlength / eye_pos[2]) + M_PI / 2;
}

void 
TrackWiimote::cwiid_callback(cwiid_wiimote_t *wiimote, int mesg_count,
			     union cwiid_mesg mesg_array[], struct timespec *timestamp)
{
  for (int i(0); i < mesg_count; ++i) 
    switch (mesg_array[i].type) {
    case CWIID_MESG_STATUS:
      TrackWiimote::getInstance()._batteryLevel =  float(mesg_array[i].status_mesg.battery) / CWIID_BATTERY_MAX;
      break;
    case CWIID_MESG_IR:
      TrackWiimote::getInstance()._ir_data = mesg_array[i].ir_mesg;
      TrackWiimote::getInstance().updateIRPositions();
      break;
    case CWIID_MESG_ERROR:
      if (cwiid_close(wiimote))
	M_throw() << "Error closing the wiimote after Error from the wiimote";

      TrackWiimote::getInstance().m_wiimote = NULL;
      break;
    default:
      break;
    }
}
