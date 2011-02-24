
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
    
namespace {
  // desired IR dot coordinate range.
// the coordinates will range from -RANGE/2 to +RANGE/2
// with 0,0 always being center of the wiimote's field of view
  const float TANWIIFOVH = 0.4142; // tangent of half the horizontal field of view, i.e. tan(WIIFOVH/2)
  const size_t WIIFOVH = 45;       // wiimote's horizontal field of view in degrees
  const size_t WIIFOVV = 37;       // wiimote's vertical field of view
  const float TANWIIFOVV = 0.3346;   // tangent of half the vertical field of view
  const float initial_dist = 65;
  const float IRBarWidth = 16.5;
  const float VRANGE = TANWIIFOVV * 2;

  //Data taken from http://wiibrew.org/wiki/Wiimote#IR_Camera. The
  //wiimote is factory calibrated such that the angle per pixel for
  //both dimensions is correct! Thus, only the calibration data is
  //needed for the X dimension
  const double WiiFOVX = (45.0f / 180.0) * M_PI; 
  //The angle corresponding to a pixel
  const double anglePerPixel =  WiiFOVX / double(CWIID_IR_X_MAX); 
 
  //Distance between the two IR sources being tracked (in cm)
  const double IRPointSeparation = 16.5;

  //We normalize by the screen dimensions
  const double ScreenHeight = 27;
  const double ScreenWidth = 48;

}
// TrackWiimote
// ==================================================================================
TrackWiimote::TrackWiimote():
  m_connect_attempts(0),
  m_max_attempts(10),
  m_bt_address(*BDADDR_ANY),
  m_wiimote(NULL),
  calibrate_request(false),
  ini_wii_to_real(2*TANWIIFOVH*initial_dist),
  cur_wii_to_real(2*TANWIIFOVH*initial_dist),
  ini_point_dist(IRBarWidth/(2 * TANWIIFOVH * initial_dist))
{
  eye_pos[0] = 40; eye_pos[1] = 0; eye_pos[2] = initial_dist;
  ini_pos[0] = 40; ini_pos[1] = 0; ini_pos[2] = initial_dist;
  for (size_t i(0); i < 4; ++i)
    for (size_t j(0); j < 2; ++j)
      ir_positions[i][j] = 0;

  for (size_t i(0); i < 4; ++i)
    ir_sizes[i]= 4;
}

int TrackWiimote::updateState()
{
  int retval = cwiid_get_state(m_wiimote, &m_state);
  updateIRPositions();
  return retval;
}

int TrackWiimote::connect()
{
  printf("Wiimote: Press 1 and 2 simultaneously to connect the Wiimote.\n");

  m_wiimote = cwiid_connect(&m_bt_address, CWIID_FLAG_CONTINUOUS|CWIID_FLAG_NONBLOCK);
  if (!m_wiimote)
    {
      printf("Unable to connect to Wiimote.\n");
      return 0;
    }
  printf("Established connection to Wiimote!\nInitializing...\n");
  if (cwiid_command(m_wiimote, CWIID_CMD_LED, CWIID_LED1_ON|CWIID_LED4_ON)!=0)
    {
      printf("Failed attempting to enable LED1 and LED4.\n");
      return 0;
    }
  if (cwiid_command(m_wiimote, CWIID_CMD_RPT_MODE,CWIID_RPT_IR|CWIID_RPT_ACC|CWIID_RPT_BTN)!=0)
    {
      printf("Failed to enable Wii functions.\n");
      return 0;
    }
  if (cwiid_get_state(m_wiimote, &m_state)!=0)
    {
      printf("Failed to get Wii remote state.\n");
      return 0;
    }
  printf("Initialization complete!\n");
  return 1;
}

//! Downloads the ir positions and returns how many were recorded
int TrackWiimote::updateIRPositions()
{
  if (m_wiimote) /* if the wiimote is still connected */
    {
      int ret = 0;
      int i;
      for (i=0; i<4; i++)
	{
	  if (m_state.ir_src[i].valid)
	    {
	      ir_positions[i][0] = m_state.ir_src[i].pos[CWIID_X];
	      ir_positions[i][1] = m_state.ir_src[i].pos[CWIID_Y];
	      if (m_state.ir_src[i].size != -1)
		ir_sizes[i] = m_state.ir_src[i].size+1;

	      ret++;
	    }
	}            
      if (ret > 1) 
	{
	  if (calibrate_request)
	    {       
	      calibrate();
	      calibrate_request = false;
	    }

	  updateHeadPos();
	}
      return ret;
    }
  return 0;
}

/* update the camera position using ir points camerapt_1 and camerapt_2 */
int TrackWiimote::updateHeadPos()
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

  //The distance to the points from the camera (Z coordinate)
  eye_pos[2] = (IRPointSeparation / 2) / std::tan(pointsAngle / 2);

  //The "position" angle of the IR points
  double Xangle = (x1 + x2) / 2;
  double Yangle = (y1 + y2) / 2;

  //Assuming that we worked out the distance to the points, and the
  //distance is a radius we can work out the distance from the centre
  eye_pos[0] = eye_pos[2] * std::sin(Xangle) / ScreenWidth;
  eye_pos[1] = eye_pos[2] * std::sin(Yangle) / ScreenHeight;

  return 1;
}


void TrackWiimote::calibrate()
{
//  const size_t pt1=0, pt2=1;
//  double cur_pt_dist =  calc_distance(ir_positions[pt1][0],
//				      ir_positions[pt1][1],
//				      ir_positions[pt2][0],
//				      ir_positions[pt2][1]);
//  
//  printf("Old distance between IR sources: %.2fcm. New: %.2fcm.\n", ini_point_dist, cur_pt_dist);
//  ini_point_dist = cur_pt_dist;
//  printf("Assuming your head is %.2f units (usually cm) from the wiimote, calibration is complete.\n",ini_pos[2]);
}
