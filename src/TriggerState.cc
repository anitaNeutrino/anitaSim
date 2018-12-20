#include "TriggerState.h"
#include "AnitaSimSettings.h"

anitaSim::TriggerState::TriggerState(const Settings* s)
  : antennaClump(s ? s->antennaclump : 0)
{
  if(s){
    ///@todo is +1 necessary?
    for(auto& vec : L1){
      vec.resize(s->NLAYERS, 0);
    }
    for(auto& vec : L2){
      vec.resize(s->NLAYERS, 0);
    }

    const int numPhi = 16; ///@todo figure me out from settings!
    for(int pol=0; pol <  Anita::NPOL; pol++){
      for(int layer=0; layer < s->NLAYERS; layer++){
	location.at(pol).at(layer).resize(numPhi,  0);
      }
      locationNadirOnly.at(pol).resize(numPhi, 0);
      
    }    
  }
  
  
}
