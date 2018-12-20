#ifndef ANITA_SIM_TRIGGER_STATE_H
#define ANITA_SIM_TRIGGER_STATE_H

#include "anita.h"

namespace anitaSim {

  class Settings;

  /**
   * @class TriggerState
   * @brief A class to sanitize GlobalTrigger function arguments.
   */

  class TriggerState {
  public:
    TriggerState(const Settings* s);

    std::array<int, Anita::NPOL> passes {{0, 0}};
    std::array<int, Anita::NPOL> L3 {{0, 0}}; ///16 bit number which says which phi sectors pass L3 V-POL

    /// These vectors were previously arrays of length Anita::NTRIGGERLAYERS_MAX.
    /// @todo how big do they actually need to be?
    std::array<std::vector<int>, Anita::NPOL> L2; // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs 
    std::array<std::vector<int>, Anita::NPOL> L1; ///For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs

    int disconesPassing; /// A relic from the past that seems to do nothing... but I'll leave it in.
    const int antennaClump; ///From settings...

    static_assert(Anita::NLAYERS_MAX <= 5, "These arrays are getting large. Think about this!");
    std::array<std::array<std::vector<int>, Anita::NLAYERS_MAX>, Anita::NPOL> location;
    std::array<std::vector<int>, Anita::NPOL> locationNadirOnly;
    // int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; /// Keep me?
    // int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; ///Keep me?
  };
}

#endif 
