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

#pragma once

#include <GL/glut.h> 
#include <GL/glu.h>

#include <string>
#include <cmath>
#include <cwiid.h> /* cwiid wii remote library */

// Head/camera information structure
// ========================================================================
typedef struct {
	double eye_pos[3]; // calculated viewpoint position relative to the screen, in cm.
	double ini_pos[3]; // initial offests of the viewpoint from the center of the screen in real space in cm
	double ini_wii_to_real; //conversion from wii dot coordinates to cm (relative to viewpoint)
				// at the initial head position. ONLY valid when the infrared sources are the initial distance
				// from the screen in real space
	double cur_wii_to_real; // conversion from wii dot coordinates to cm at the distance the camera is
				// calculated to be from the screen
	double ini_point_dist; // distance between IR points in wiis
	bool calibrate_request; // true if calibration requested by user, (default is false)
} head_info_t;

// TrackWiimote
// ==================================================================================
class TrackWiimote {
private:
    int m_connect_attempts; /* number of connection attempts */
    int m_max_attempts; /* max number of connection attempts */
    bdaddr_t m_bt_address; /* bluetooth address */
    cwiid_wiimote_t *m_wiimote; /* wiimote connection through cwiid library */
    struct cwiid_state m_state; /* wiimote state (updated every frame) using cwiid library */

public:
    TrackWiimote();

    int updateState();

    int connect();

    int updateIRPositions(double ir_positions[][2], int ir_sizes[]);

    int updateHeadPos(head_info_t *head, double ir_positions[][2], int ir_sizes[],
				int pt1, int pt2);

    void calibrate(head_info_t *head, double ir_positions[][2], int ir_sizes[],
				    int pt1, int pt2);
};
