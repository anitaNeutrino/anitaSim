#include "AnitaSimSettings.h"
#include "CommandLineOptions.h"
#include "EventGenerator.h"
#include "AnitaPayload.h"

int main(int argc,  char **argv) {

  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------
  anitaSim::Settings settings;
  icemc::CommandLineOptions clOpts(argc, argv, settings);

  if(clOpts.are_good){
    auto anita = std::make_shared<anitaSim::AnitaPayload>(&settings);    
    icemc::EventGenerator uhen(&settings);
    uhen.generate(*anita.get());
  }
  
  return 0;
}
