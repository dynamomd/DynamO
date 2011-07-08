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

#include <magnet/exception.hpp>
#include <magnet/thread/mutex.hpp>
#include <magnet/math/vector.hpp>
#include <cwiid.h> /* cwiid wii remote library */
#include <vector>

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

  inline Vector getHeadPosition()
  { 
    magnet::thread::ScopedLock lock(_irdatalock);
    return eye_pos;
  }

  inline bool connected() const { return m_wiimote; }

  inline volatile const float& getBatteryLevel() const { return _batteryLevel; }

  struct IRData {
    IRData(int8_t _size, uint16_t _x, uint16_t _y): size(_size), x(_x), y(_y) {}
    bool operator>(const IRData& o) const { return size > o.size; }
    int8_t size; uint16_t x; uint16_t y; 
  };

  inline std::vector<IRData> getSortedIRData()
  { 
    magnet::thread::ScopedLock lock(_irdatalock);
    return _irdata;
  }

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
  volatile double v_angle;
  volatile float _batteryLevel;

  std::vector<IRData> _irdata;

  magnet::thread::Mutex _irdatalock;
};
