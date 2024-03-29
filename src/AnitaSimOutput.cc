#include "AnitaSimOutput.h"
#include "AnitaSimSettings.h"
#include "EventGenerator.h"
#include "RootOutput.h"
#include "Report.h"
#include "EnvironmentVariable.h"
#include "AnitaPayload.h"
#include "RayTracer.h"
#include "SurfaceScreen.h"
#include "Tools.h"

#include "TFile.h"

#ifdef ANITA_UTIL_EXISTS
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"
#include "Adu5Pat.h"
#include "FFTtools.h"

#ifdef ANITA3_EVENTREADER
#include "TruthAnitaEvent.h"
#endif

#endif


anitaSim::AnitaSimOutput::AnitaSimOutput(const AnitaPayload* detector, const Settings* settings, const char* outputDir, int run)
  : fDetector(detector), fSettings(settings),
    fRun(run), fOutputDir(TString::Format("%s/run%d/", outputDir, run)), 
    fEvent(nullptr), fHeader(nullptr), fGps(nullptr)
{
  TString mkdirCommand("mkdir -p ");
  mkdirCommand += fOutputDir;
  system(mkdirCommand.Data());

  initRootifiedAnitaDataFiles();
}





anitaSim::AnitaSimOutput::~AnitaSimOutput(){

  // write, close and delete all non-NULL member files.
  for(auto f : {fHeadFile, fGpsFile, fEventFile, fTruthFile}){
    if(f){
      f->Write();
      f->Close();
      delete f;
    }
  }
}



void anitaSim::AnitaSimOutput::initRootifiedAnitaDataFiles(){

#ifdef ANITA_UTIL_EXISTS

  TString eventFileName = fOutputDir + TString::Format("/SimulatedAnitaEventFile%d.root", fRun);
  fEventFile = new TFile(eventFileName, "RECREATE");

  icemc::RootOutput::initTree(&eventTree, "eventTree", "eventTree", fEventFile);
  eventTree.Branch("event",             &fEvent           );
  eventTree.Branch("run",               fRun);//,   "run/I"   );  
  // eventTree.Branch("weight",            &uhen->fNeutrinoPath->weight,   "weight/D"); ///@todo restore weight here 

  TString headFileName = fOutputDir + TString::Format("/SimulatedAnitaHeadFile%d.root", fRun);
  fHeadFile = new TFile(headFileName, "RECREATE");
  
  icemc::RootOutput::initTree(&headTree, "headTree", "headTree", fHeadFile);
  headTree.Branch("header",  &fHeader           );
  // headTree.Branch("weight",  &uhen->fNeutrinoPath->weight,      "weight/D"); ///@todo restore weight here

  TString gpsFileName = fOutputDir + TString::Format("/SimulatedAnitaGpsFile%d.root", fRun);
  fGpsFile = new TFile(gpsFileName, "RECREATE");

  icemc::RootOutput::initTree(&adu5PatTree, "adu5PatTree", "adu5PatTree", fGpsFile);
  adu5PatTree.Branch("pat",          &fGps);
  adu5PatTree.Branch("eventNumber",  &fEventNumber,  "eventNumber/I");
  // adu5PatTree.Branch("weight",       &uhen->fNeutrinoPath->weight,       "weight/D"     ); ///@todo restore weight here
  
#ifdef ANITA3_EVENTREADER
  // Set AnitaVersion so that the right payload geometry is used
  AnitaVersion::set(fSettings->ANITAVERSION);
  
  TString truthFileName = fOutputDir + TString::Format("/SimulatedAnitaTruthFile%d.root",  fRun);
  fTruthFile = new TFile(truthFileName, "RECREATE");
  TNamed* ss = fSettings->makeRootSaveableSettings();
  ss->Write();
  delete ss;
  ss = NULL;

  TNamed gitversion("GitVersion", icemc::EnvironmentVariable::ICEMC_VERSION(fOutputDir).c_str());
  gitversion.Write();

  TNamed nnu("NNU", TString::Format("%d", fSettings->NNU).Data());
  nnu.Write();

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  TNamed starttime("StartTime", asctime(timeinfo));
  starttime.Write();

  icemc::report() << icemc::severity::info << "ICEMC GIT Repository Version: " <<  gitversion.GetTitle() << std::endl;

  icemc::RootOutput::initTree(&triggerSettingsTree, "triggerSettingsTree", "Trigger settings", fTruthFile);
  // Anita* anita2 = const_cast<Anita*>(uhen->anita1);
  AnitaPayload* anita2 = const_cast<AnitaPayload*>(fDetector);
  
  triggerSettingsTree.Branch("dioderms",  anita2->bwslice_dioderms_fullband_allchan,  "dioderms[2][48][6]/D" );
  triggerSettingsTree.Branch("diodemean", anita2->bwslice_diodemean_fullband_allchan, "diodemean[2][48][6]/D");
  triggerSettingsTree.Fill();
  
  icemc::RootOutput::initTree(&truthTree, "truthAnitaTree", "Truth Anita Tree", fTruthFile);
  truthTree.Branch("truth", &fTruth);

  //icemc::RootOutput::initTree(&summaryTree, "summaryAnitaTree", "Summary Anita Tree", fTruthFile);
  //summaryTree.Branch("summary", &fSummary);

#else
  Log() << severity::warning << "Need ANITA EventReader version at least 3 to produce Truth output." << std::endl;    
#endif

#else
  Log() << severity::warning << "Can't generate ROOTified output without satisfying ANITA_UTIL_EXISTS at compile time" << std::endl;
#endif
}







void anitaSim::AnitaSimOutput::fillRootifiedAnitaDataTrees(const icemc::Event& icemcEvent){
  
#ifdef ANITA_UTIL_EXISTS
  AnitaGeomTool* geom = AnitaGeomTool::Instance();

  const AnitaPayload* bn1 = fDetector;
  const AnitaPayload* anita1 = fDetector;

  fEvent  = new UsefulAnitaEvent();
  fHeader = new RawAnitaHeader();
  Adu5Pat gps = fDetector->pat();
  fGps    = &gps;
  fGps->run = fRun;

  memset(fEvent->fNumPoints, 0, sizeof(fEvent->fNumPoints) );
  memset(fEvent->fVolts,     0, sizeof(fEvent->fVolts)     );
  memset(fEvent->fTimes,     0, sizeof(fEvent->fTimes)     );

  int fNumPoints = NUM_SAMP;
  for (int ichan=0; ichan<NUM_DIGITZED_CHANNELS; ichan++){
    fEvent->fNumPoints[ichan] = fNumPoints;

    for (int j = 0; j < fNumPoints; j++) {
      // convert seconds to nanoseconds
      fEvent->fTimes[ichan][j] = j * anita1->TIMESTEP * 1.0E9;
    }
  }
  fEvent->fRFSpike = 0;// when we get as far as simulating this, we're doing well

  for (int iant = 0; iant < fSettings->NANTENNAS; iant++){
    int UsefulChanIndexH = geom->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
    int UsefulChanIndexV = geom->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);
    fEvent->fNumPoints[UsefulChanIndexV] = fNumPoints;
    fEvent->fNumPoints[UsefulChanIndexH] = fNumPoints;
    fEvent->chanId[UsefulChanIndexV] = UsefulChanIndexV;
    fEvent->chanId[UsefulChanIndexH] = UsefulChanIndexH;

    const int offset = (fDetector->fVoltsRX.channelsH.at(iant).rfcm_lab_all.size() - fNumPoints)/2; ///@todo find a better way to do this...
    for (int j = 0; j < fNumPoints; j++) {
      // convert seconds to nanoseconds
      fEvent->fTimes[UsefulChanIndexV][j] = j * anita1->TIMESTEP * 1.0E9;
      fEvent->fTimes[UsefulChanIndexH][j] = j * anita1->TIMESTEP * 1.0E9;

      const double voltsToMilliVolts = 1000; // volts to millivolts
      // fEvent->fVolts[UsefulChanIndexH][j] = fDetector->fVoltsRX.rfcm_lab_h_all[iant][j+offset]*voltsToMilliVolts;
      // fEvent->fVolts[UsefulChanIndexV][j] = fDetector->fVoltsRX.rfcm_lab_e_all[iant][j+offset]*voltsToMilliVolts;
      fEvent->fVolts[UsefulChanIndexH][j] = fDetector->fVoltsRX.channelsH.at(iant).rfcm_lab_all[j+offset]*voltsToMilliVolts;
      fEvent->fVolts[UsefulChanIndexV][j] = fDetector->fVoltsRX.channelsV.at(iant).rfcm_lab_all[j+offset]*voltsToMilliVolts;
      
      fEvent->fCapacitorNum[UsefulChanIndexH][j] = j;
      fEvent->fCapacitorNum[UsefulChanIndexV][j] = j;
    }//end int j
  }// end int iant

  fEventNumber = icemcEvent.loop.eventNumber;
  fEvent->eventNumber = icemcEvent.loop.eventNumber;
  fHeader->eventNumber = icemcEvent.loop.eventNumber;
  fHeader->surfSlipFlag = 0;
  fHeader->errorFlag = 0;

  if (fSettings->MINBIAS==1){
    fHeader->trigType = 8; // soft-trigger
  }
  else{
    fHeader->trigType = 1; // RF trigger
  }

  fHeader->run = icemcEvent.loop.run;
  // put the vpol only as a placeholder - these are only used in Anita-2 anyway

  fHeader->upperL1TrigPattern = fDetector->fTriggerState.L1[0][0];
  fHeader->lowerL1TrigPattern = fDetector->fTriggerState.L1[0][1];
  fHeader->nadirL1TrigPattern = fDetector->fTriggerState.L1[0][2];
  fHeader->upperL2TrigPattern = fDetector->fTriggerState.L2[0][0];
  fHeader->lowerL2TrigPattern = fDetector->fTriggerState.L2[0][1];
  fHeader->nadirL2TrigPattern = fDetector->fTriggerState.L2[0][2];

  // fHeader->upperL1TrigPattern = fDetector->fL1trig[0][0];
  // fHeader->lowerL1TrigPattern = fDetector->fL1trig[0][1];
  // fHeader->nadirL1TrigPattern = fDetector->fL1trig[0][2];
  // fHeader->upperL2TrigPattern = fDetector->fL2trig[0][0];
  // fHeader->lowerL2TrigPattern = fDetector->fL2trig[0][1];
  // fHeader->nadirL2TrigPattern = fDetector->fL2trig[0][2];
  
  if (fSettings->WHICH<Payload::Anita3){
    fHeader->phiTrigMask  = (short) anita1->phiTrigMask;
    // fHeader->l3TrigPattern = (short) fDetector->fL3trig[0];
    fHeader->l3TrigPattern = (short) fDetector->fTriggerState.L3[0];    
    // fHeader->l3TrigPattern = (short) uhen->l3trig[0];        
  }

  fHeader->calibStatus = 31;
  fHeader->realTime = bn1->getRealTime();//realTime_flightdata;
  fHeader->triggerTime = bn1->getRealTime(); //realTime_flightdata;

#ifdef ANITA3_EVENTREADER
  if (fSettings->WHICH==Payload::Anita3 || fSettings->WHICH==Payload::Anita4) {
    fHeader->setTrigPattern((short) fDetector->fTriggerState.L3.at(0), AnitaPol::kVertical);
    fHeader->setTrigPattern((short) fDetector->fTriggerState.L3.at(1), AnitaPol::kHorizontal);
    fHeader->setMask( (short) anita1->l1TrigMask,  (short) anita1->phiTrigMask,  AnitaPol::kVertical);
    fHeader->setMask( (short) anita1->l1TrigMaskH, (short) anita1->phiTrigMaskH, AnitaPol::kHorizontal);
  }

  fTruth                   = new TruthAnitaEvent();
  fTruth->eventNumber      = icemcEvent.loop.eventNumber;
  fTruth->realTime         = bn1->getRealTime(); //realTime_flightdata;
  fTruth->run              = icemcEvent.loop.run; //clOpts.run_no;

  //@todo URGENT RESTORE these parameters to the truth tree! FIX ME!  
  fTruth->nuMom            = icemcEvent.neutrino.energy.in(icemc::Energy::Unit::eV);
  fTruth->nu_pdg           = icemcEvent.neutrino.pdgCode();
  // fTruth->e_component      = fDetector->e_component;
  // fTruth->h_component      = fDetector->h_component;
  // fTruth->n_component      = fDetector->n_component;
  // fTruth->e_component_k    = fDetector->e_component_kvector;
  // fTruth->h_component_k    = fDetector->h_component_kvector;
  // fTruth->n_component_k    = fDetector->n_component_kvector;
  fTruth->sourceLon        = icemcEvent.interaction.position.Longitude();
  fTruth->sourceLat        = icemcEvent.interaction.position.Latitude();
  fTruth->sourceAlt        = icemcEvent.interaction.position.Altitude();
  fTruth->weight           = icemcEvent.loop.phaseWeight();
  TVector3 n_bn = bn1->position().Unit();
  TVector3 n_int_pos = icemcEvent.interaction.position.Unit();
  for (int i=0;i<3;i++){
    fTruth->balloonPos[i]  = bn1->position()[i];
    fTruth->balloonDir[i]  = n_bn[i];
    fTruth->nuPos[i] = icemcEvent.interaction.position[i];
    fTruth->nuDir[i] = n_int_pos[i];
  }
  // for (int i=0;i<5;i++){
  //   for (int j=0;j<3;j++){
  //     fTruth->rfExitNor[i][j] = ray1->n_exit2bn[i][j];
  //     fTruth->rfExitPos[i][j] = ray1->rfexit[i][j];
  //   }
  // }
  // for (int i=0;i<48;i++){
  // @todo URGENT RESTORE THESE PARAMETERS!
  // fTruth->hitangle_e[i]  = uhen->hitangle_e_all[i];
  // fTruth->hitangle_h[i]  = uhen->hitangle_h_all[i];
  // }
  // if(!settings1.ROUGHNESS){
  //   for (int i=0;i<Anita::NFREQ;i++){
  //     fTruth->vmmhz[i]       = panel1->GetVmmhz_freq(i);
  //   }
  // }

	    
  memset(fTruth->SNRAtTrigger,       0, sizeof(fTruth->SNRAtTrigger)       );
  memset(fTruth->fSignalAtTrigger,   0, sizeof(fTruth->fSignalAtTrigger)   );
  memset(fTruth->fNoiseAtTrigger,    0, sizeof(fTruth->fNoiseAtTrigger)    );
  memset(fTruth->SNRAtDigitizer,     0, sizeof(fTruth->SNRAtDigitizer)     );
  memset(fTruth->thresholds,         0, sizeof(fTruth->thresholds)         );
  memset(fTruth->fDiodeOutput,       0, sizeof(fTruth->fDiodeOutput)       );
	    
  fTruth->maxSNRAtTriggerV=0;
  fTruth->maxSNRAtTriggerH=0;
  fTruth->maxSNRAtDigitizerV=0;
  fTruth->maxSNRAtDigitizerH=0;

  for (int iant = 0; iant < fSettings->NANTENNAS; iant++){
    int UsefulChanIndexH = geom->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
    int UsefulChanIndexV = geom->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);

    fTruth->SNRAtTrigger[UsefulChanIndexV] = icemc::Tools::calculateSNR(fDetector->fVoltsRX.channelsV.at(iant).justSignalTrig, fDetector->fVoltsRX.channelsV.at(iant).justNoiseTrig);
    fTruth->SNRAtTrigger[UsefulChanIndexH] = icemc::Tools::calculateSNR(fDetector->fVoltsRX.channelsH.at(iant).justSignalTrig, fDetector->fVoltsRX.channelsH.at(iant).justNoiseTrig);
	      
    if (fTruth->SNRAtTrigger[UsefulChanIndexV]>fTruth->maxSNRAtTriggerV) fTruth->maxSNRAtTriggerV=fTruth->SNRAtTrigger[UsefulChanIndexV];
    if (fTruth->SNRAtTrigger[UsefulChanIndexH]>fTruth->maxSNRAtTriggerH) fTruth->maxSNRAtTriggerH=fTruth->SNRAtTrigger[UsefulChanIndexH];

    fTruth->SNRAtDigitizer[UsefulChanIndexV] = icemc::Tools::calculateSNR(fDetector->fVoltsRX.channelsV.at(iant).justSignalDig, fDetector->fVoltsRX.channelsV.at(iant).justNoiseDig);
    fTruth->SNRAtDigitizer[UsefulChanIndexH] = icemc::Tools::calculateSNR(fDetector->fVoltsRX.channelsH.at(iant).justSignalDig, fDetector->fVoltsRX.channelsH.at(iant).justNoiseDig);

    if (fTruth->SNRAtDigitizer[UsefulChanIndexV]>fTruth->maxSNRAtDigitizerV) fTruth->maxSNRAtDigitizerV=fTruth->SNRAtDigitizer[UsefulChanIndexV];
    if (fTruth->SNRAtDigitizer[UsefulChanIndexH]>fTruth->maxSNRAtDigitizerH) fTruth->maxSNRAtDigitizerH=fTruth->SNRAtDigitizer[UsefulChanIndexH];

    fTruth->thresholds[UsefulChanIndexV] = fDetector->fThresholdsAnt[iant][0][4];
    fTruth->thresholds[UsefulChanIndexH] = fDetector->fThresholdsAnt[iant][1][4];
    int irx = iant;
    if (iant<16){
      if (iant%2) irx = iant/2;
      else        irx = iant/2 + 1;
    }
	      
    for (int j = 0; j < fNumPoints; j++) {
      fTruth->fTimes[UsefulChanIndexV][j]             = j * anita1->TIMESTEP * 1.0E9;
      fTruth->fTimes[UsefulChanIndexH][j]             = j * anita1->TIMESTEP * 1.0E9;

      ///@todo replace this 128 with offset as in event->fVolts filling loop
      const int offset = 128;
      fTruth->fSignalAtTrigger[UsefulChanIndexV][j]   = fDetector->fVoltsRX.channelsV.at(iant).justSignalTrig[j+offset]*1000;
      fTruth->fSignalAtTrigger[UsefulChanIndexH][j]   = fDetector->fVoltsRX.channelsH.at(iant).justSignalTrig[j+offset]*1000;
      fTruth->fNoiseAtTrigger[UsefulChanIndexV][j]    = fDetector->fVoltsRX.channelsV.at(iant).justNoiseTrig[j+offset]*1000;
      fTruth->fNoiseAtTrigger[UsefulChanIndexH][j]    = fDetector->fVoltsRX.channelsH.at(iant).justNoiseTrig[j+offset]*1000;
      fTruth->fSignalAtDigitizer[UsefulChanIndexV][j] = fDetector->fVoltsRX.channelsV.at(iant).justSignalDig[j+offset]*1000;
      fTruth->fSignalAtDigitizer[UsefulChanIndexH][j] = fDetector->fVoltsRX.channelsH.at(iant).justSignalDig[j+offset]*1000;
      fTruth->fNoiseAtDigitizer[UsefulChanIndexV][j]  = fDetector->fVoltsRX.channelsV.at(iant).justNoiseDig[j+offset]*1000;
      fTruth->fNoiseAtDigitizer[UsefulChanIndexH][j]  = fDetector->fVoltsRX.channelsH.at(iant).justNoiseDig[j+offset]*1000;
      fTruth->fDiodeOutput[UsefulChanIndexV][j]       = anita1->timedomain_output_allantennas[0][irx][j];
      fTruth->fDiodeOutput[UsefulChanIndexH][j]       = anita1->timedomain_output_allantennas[1][irx][j];
    }//end int j
  }// end int iant

  truthTree.Fill();
  delete fTruth;
  fTruth = nullptr;
#endif

  headTree.Fill();
  eventTree.Fill();
  adu5PatTree.Fill();

  std::cout << "HeadTree has " << headTree.GetEntries() << " entries " << std::endl;

  delete fEvent;
  delete fHeader;
  // delete fGps; // now on stack, don't delete!

  fEvent = nullptr;
  fHeader = nullptr;
  fGps = nullptr;
  
#endif

}
