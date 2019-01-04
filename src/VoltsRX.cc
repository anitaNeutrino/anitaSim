#include "VoltsRX.h"

anitaSim::ChannelVolts::ChannelVolts(){

}


anitaSim::ChannelVolts::~ChannelVolts(){

}


void anitaSim::ChannelVolts::reset(){
  rfcm_lab_all.fill(0);
  justNoiseTrig.fill(0);
  justSignalTrig.fill(0);
  justNoiseDig.fill(0);
  justSignalDig.fill(0);
}

anitaSim::VoltsRX::VoltsRX(int nRX) : fNRX(nRX) {

  max = 0;
  ave = 0;
  sum = 0;

  max_highband = 0; // max voltage seen on an antenna - just for debugging purposes
  max_lowband = 0; // max voltage seen on an antenna - just for debugging purposes


  channelsV.reserve(nRX);
  channelsH.reserve(nRX);

  for(int i=0; i < fNRX; i++){
    channelsV.emplace_back(ChannelVolts());
  }
  for(int i=0; i < fNRX; i++){
    channelsH.emplace_back(ChannelVolts());
  }
  reset();
}


void anitaSim::VoltsRX::reset() {

  for(auto& chan : channelsV){chan.reset();}
  for(auto& chan : channelsH){chan.reset();}
  
}
