
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

#include <GL/glut.h> 
#include <GL/glu.h>
    
#include "wiiheadtracking.hpp"
    
namespace {
  // desired IR dot coordinate range.
// the coordinates will range from -RANGE/2 to +RANGE/2
// with 0,0 always being center of the wiimote's field of view
  const size_t HRANGE = 100;
  const float TANWIIFOVH = 0.4142; // tangent of half the horizontal field of view, i.e. tan(WIIFOVH/2)
  const size_t WIIFOVH = 45;       // wiimote's horizontal field of view in degrees
  const size_t WIIFOVV = 37;       // wiimote's vertical field of view
  const float TANWIIFOVV = 0.3346;   // tangent of half the vertical field of view
  const float initial_dist = 65;
  const float IRBarWidth = 16.5;
  const float VRANGE = HRANGE * TANWIIFOVV * 2;

  double calc_distance(double x1, double y1, double x2, double y2)
  {
    double xdiff = x1-x2;
    double ydiff = y1-y2;
    double distsqr = (xdiff*xdiff)+(ydiff*ydiff);
    return std::sqrt(distsqr);
  }
}
// TrackWiimote
// ==================================================================================
TrackWiimote::TrackWiimote():
  m_connect_attempts(0),
  m_max_attempts(10),
  m_bt_address(*BDADDR_ANY),
  m_wiimote(NULL),
  ini_wii_to_real(2*TANWIIFOVH*initial_dist/HRANGE),
  cur_wii_to_real(2*TANWIIFOVH*initial_dist/HRANGE),
  ini_point_dist(HRANGE*IRBarWidth/(2 * TANWIIFOVH * initial_dist))
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
	      ir_positions[i][0] = ((m_state.ir_src[i].pos[CWIID_X] * HRANGE
					   / (double)CWIID_IR_X_MAX) -(HRANGE/2));
	      ir_positions[i][1] = (m_state.ir_src[i].pos[CWIID_Y] * VRANGE
				    / (double)CWIID_IR_Y_MAX) -(VRANGE/2);
	      if (m_state.ir_src[i].size != -1)
		ir_sizes[i] = m_state.ir_src[i].size+1;

	      ret++;
	    }
	}
      if (ret > 1) updateHeadPos();
      
      return ret;
    }
  return 0;
}

/* update the camera position using ir points camerapt_1 and camerapt_2 */
int TrackWiimote::updateHeadPos()
{
  const size_t pt1=0, pt2=1;
  double cur_pt_dist = calc_distance(ir_positions[pt1][0],
				     ir_positions[pt1][1],
				     ir_positions[pt2][0],
				     ir_positions[pt2][1]);
  double x = (ir_positions[pt1][0] + ir_positions[pt2][0]) / 2;
  double y = (ir_positions[pt1][1] + ir_positions[pt2][1]) / 2;

  double z_scale = (ini_point_dist/cur_pt_dist);

  if (z_scale > 100 || cur_pt_dist < 0.01 || cur_pt_dist > 100.0)
    return 0;

  cur_wii_to_real = z_scale * ini_wii_to_real; // is this needed?
  eye_pos[0] = x * cur_wii_to_real;
  eye_pos[1] = y * cur_wii_to_real;
  eye_pos[2] = z_scale * ini_pos[2];
  return 1;
}


void TrackWiimote::calibrate()
{
  const size_t pt1=0, pt2=1;
  double cur_pt_dist =  calc_distance(ir_positions[pt1][0],
				      ir_positions[pt1][1],
				      ir_positions[pt2][0],
				      ir_positions[pt2][1]);

  printf("Old distance between IR sources: %.2fcm. New: %.2fcm.\n", ini_point_dist, cur_pt_dist);
  ini_point_dist = cur_pt_dist;
  printf("Assuming your head is %.2f units (usually cm) from the wiimote, calibration is complete.\n",ini_pos[2]);
}
