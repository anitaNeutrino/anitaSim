////////////
// rx.h //
////////////

#ifndef ANITA_SIM_RX_H
#define ANITA_SIM_RX_H

// Standard Library #includes
#include <vector>

// ROOT Library #includes
#include "TObject.h"

// from RVersion.h
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
#include "TClingRuntime.h"
#endif
#else
#include "TCint.h"
#endif

namespace anitaSim{
  class RX : public TObject {
  public:
    unsigned phi_sector;
    unsigned layer;
    double x, y, z;
    std::vector <double>* waveform;
    std::vector <double>* digitized;
    RX (void) : waveform (NULL), digitized (NULL) {}
    ~RX (void) {delete waveform; delete digitized;}
    ClassDef(RX, 1);
  };
}

#endif
