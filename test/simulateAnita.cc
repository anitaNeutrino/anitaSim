#include "AnitaSimSettings.h"
#include "CommandLineOptions.h"
#include "EventGenerator.h"
#include "AnitaPayload.h"

#include <chrono>

int main(int argc,  char **argv) {

  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------
  auto starttime = std::chrono::system_clock::now();
  anitaSim::Settings settings;
  icemc::CommandLineOptions clOpts(argc, argv, settings);

  if(clOpts.are_good){
    auto anita = std::make_shared<anitaSim::AnitaPayload>(&settings);    
    icemc::EventGenerator uhen(&settings);
    uhen.generate(*anita.get());
  }

  auto endtime = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime);
  int elapsed_seconds = elapsed.count();
  std::cout << "Elapsed runtime is " << (int)(elapsed_seconds/60) << ":" << (elapsed_seconds%60) << " minutes" <<  std::endl;
  return 0;
}
