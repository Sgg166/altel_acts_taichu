#pragma once

#include "eudaq/RawEvent.hh"
#include "TelEvent.hpp"

namespace altel{

  std::shared_ptr<TelEvent> createTelEvent(eudaq::EventSPC eudaqEvent);

}
