#include <dynamo/eventtypes.hpp>
#include <magnet/exception.hpp>

namespace dynamo {

#define printEnum(VAL)                                                         \
  case VAL:                                                                    \
    return os << #VAL;

inline std::ostream &operator<<(std::ostream &os, EEventType etype) {
  switch (etype) {
    ETYPE_ENUM_FACTORY(printEnum)
  case FINAL_ENUM_TO_CATCH_THE_COMMA:
  default:
    break;
  }

  M_throw() << "Failed to find a name for the Event Type! The value must be "
               "uninitialised or memory is being corrupted somehow.";
}
inline std::ostream &operator<<(std::ostream &os, EventSource etype) {
  switch (etype) {
  case INTERACTION:
    os << "Interaction";
    break;
  case LOCAL:
    os << "Local";
    break;
  case GLOBAL:
    os << "Global";
    break;
  case SYSTEM:
    os << "System";
    break;
  case SCHEDULER:
    os << "Scheduler";
    break;
  case NOSOURCE:
    os << "No-Source";
    break;
  default:
    M_throw() << "Failed to find a name for the Event source type! The value "
                 "must be uninitialised or memory is being corrupted somehow.";
    break;
  };
  return os;
}

Event::Event()
    : _dt(std::numeric_limits<float>::infinity()),
      _particle1ID(std::numeric_limits<size_t>::max()),
      _sourceID(std::numeric_limits<size_t>::max()),
      _additionalData1(std::numeric_limits<size_t>::max()),
      _additionalData2(std::numeric_limits<size_t>::max()), _source(NOSOURCE),
      _type(NONE) {}

Event::Event(size_t particle1ID, double dt, EventSource source, EEventType type,
             size_t sourceID, size_t additionalData1,
             size_t additionalData2) throw()
    : _dt(dt), _particle1ID(particle1ID), _sourceID(sourceID),
      _additionalData1(additionalData1), _additionalData2(additionalData2),
      _source(source), _type(type) {}

bool Event::operator<(const Event &o) const throw() { return _dt < o._dt; }

bool Event::operator>(const Event &o) const throw() { return _dt > o._dt; }

bool Event::operator==(const Event &o) const throw() {
  return (_dt == o._dt) && (_particle1ID == o._particle1ID) &&
         (_sourceID == o._sourceID) && (_type == o._type) &&
         (_additionalData1 == o._additionalData1) && (_source == o._source) &&
         ((_source == INTERACTION) || (_additionalData2 == o._additionalData2));
}

std::ostream &operator<<(std::ostream &os, Event event) {
  os << "Event{dt = " << event._dt << ", p1ID = " << event._particle1ID
     << ", sourceID = " << event._sourceID
     << ", data1 = " << event._additionalData1
     << ", data2 = " << event._additionalData2 << ", source = " << event._source
     << ", type = " << event._type << "}";
  return os;
}
} // namespace dynamo