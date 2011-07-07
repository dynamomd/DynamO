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

#ifdef COIL_wiimote
#include <magnet/exception.hpp>
#include <magnet/math/vector.hpp>
#include <cwiid.h> /* cwiid wii remote library */

/*! \brief A class to facilitate head tracking using the cwiid
 * library.
 *
 */
class TrackWiimote {
public:
  static TrackWiimote& getInstance()
  {
    static TrackWiimote _singleton;
    return _singleton;
  }

  bool connect(bdaddr_t* bt_address = BDADDR_ANY);

  void calibrate();

  inline const Vector& getHeadPosition() const { return eye_pos; }

  inline bool connected() const { return m_wiimote; }

  inline float getBatteryLevel() const { return _batteryLevel; }

  const cwiid_ir_src& getIRState(size_t IRDotID) const 
  { 
    if (IRDotID >= CWIID_IR_SRC_COUNT)
      M_throw() << "Out of bounds access on wii IR data.";
    return _ir_data.src[IRDotID]; 
  }

  inline size_t getValidIRSources() const { return _valid_ir_points; }

  inline double getCalibrationAngle() const { return v_angle; }

private:
  static void cwiid_callback(cwiid_wiimote_t *wiimote, int mesg_count,
			     union cwiid_mesg mesg_array[], struct timespec *timestamp);

  TrackWiimote();
  ~TrackWiimote();

  size_t updateIRPositions();

  void updateHeadPos();

  cwiid_wiimote_t* m_wiimote; 
  Vector eye_pos;
  double v_angle;
  size_t _valid_ir_points;

  float _batteryLevel;
  struct cwiid_ir_mesg _ir_data;

};
#endif
