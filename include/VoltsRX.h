#ifndef ANITA_SIM_VOLTS_RX_H
#define ANITA_SIM_VOLTS_RX_H

#include "TObject.h"
#include "anita.h"

namespace anitaSim {


  class ChannelVolts {
  public:
    ChannelVolts();
    ~ChannelVolts();
    void reset();

    std::array<double, Anita::HALFNFOUR> rfcm_lab_all;
    std::array<double, Anita::HALFNFOUR> justNoiseTrig;
    std::array<double, Anita::HALFNFOUR> justSignalTrig;
    std::array<double, Anita::HALFNFOUR> justNoiseDig;
    std::array<double, Anita::HALFNFOUR> justSignalDig;
    ClassDef(ChannelVolts, 1)
  };

  /**
   * @class VoltsRX
   * @brief Voltage seen at the antennas
   */
  class VoltsRX {
  public:
    VoltsRX(int nRX);
    void reset();
    
    double max;		///< max voltage seen on an antenna - just for debugging purposes
    double ave;		///< ave voltage seen on an antenna,  among hit antennas
    double sum;		///< ave voltage seen on an antenna,  among hit antennas
    double max_highband;///< max voltage seen on an antenna - just for debugging purposes
    double max_lowband;	///< max voltage seen on an antenna - just for debugging purposes

    const int fNRX;
    std::vector<ChannelVolts> channelsV;
    std::vector<ChannelVolts> channelsH;
    
    ClassDef(VoltsRX, 3);
  };
}

#endif
