#include "ANITA.h"
#include "AnitaSimSettings.h"
#include "GlobalTrigger.h"
#include "ChanTrigger.h"
#include "Tools.h"
#include "AnitaSimSettings.h"
#include "Report.h"
#include "Constants.h"
#include "RayTracer.h"
#include "VoltsRX.h"
#include "Seavey.h"
#include "RootOutput.h"
#include <memory>

#include "TFile.h" ///@todo remove after done debugging

anitaSim::ANITA::ANITA(const Settings* settings)
// anitaSim::ANITA::ANITA(const Settings* settings, const RootOutput* ro ) 
  : FlightDataManager(settings), Anita(settings, settings->getOutputDir(), this),
    fSettings(settings),
    fNumRX(settings ? settings->NANTENNAS : 48),
    fVoltsRX(settings ? settings->NANTENNAS : 0),
    fTriggerState(settings),
    fAnitaOutput(this, settings, settings->getOutputDir(), settings->getRun())
{
  initSeaveys(settings, this);

  if(fSettings->WHICH == Payload::Anita1Simple ||
     fSettings->WHICH == Payload::Anita1){
    SetDiffraction(); // for the upper ring
  }

  saveGainsPlot(std::string(fSettings->getOutputDir())+"/gains.eps");
  
  // sets position of balloon and related quantities
  // SetDefaultBalloonPosition();
  // SetNoise(fSettings, this, antarctica);  
}


anitaSim::ANITA::~ANITA(){

  
}


bool anitaSim::ANITA::chanceInHell(const icemc::PropagatingSignal& signal){
  ///@todo do something much, much cleverer here...
  /// it's not even clear this is anywhere near the correct value
  if(signal.maxEField() > 1e-5){
    return true;
  }
  else {
    return false;
  }
  return true;
}


const Geoid::Position& anitaSim::ANITA::getPosition(double time){
  /**
   * Here I'm using TMath's QuietNaN() as a default argument.
   * That might be a bad idea, but for now...
   */

  if(TMath::IsNaN(time) && TMath::IsNaN(fLastPositionTime)){
    time = getStartTime();
  }

  if(!TMath::IsNaN(time) && time != fLastPositionTime){

    getBalloonPosition(time, this);

    for(auto& s : fSeaveys){
      s.updatePosition(FlightDataManager::position(),
		       FlightDataManager::getHeading(),
		       FlightDataManager::getPitch(),
		       FlightDataManager::getRoll());
    }

    fLastPositionTime = time;
  }
  return FlightDataManager::position();
}


TVector3 anitaSim::ANITA::getPositionRX(Int_t rx) const {
  return fSeaveys.at(rx).getPosition(Seavey::Pol::V);
}



void anitaSim::ANITA::initSeaveys(const Settings *settings1, const Anita *anita1) {

  for(int rx = 0; rx < getNumRX(); rx++){
    int ilayer = -1;
    int ifold = -1;
    getLayerFoldFromTriggerRX(rx, ilayer, ifold);
    // std::cout << "Seaveys: rx = " << rx << ", ilayer = " << ilayer << ", ifold = " << ifold << ",  rxTrigger = " << this->GetRxTriggerNumbering(ilayer, ifold) << std::endl;
    

    TVector3 n_eplane;
    TVector3 n_hplane;
    TVector3 n_normal;
    
    if(settings1->WHICH==Payload::Anita1 ||
       settings1->WHICH==Payload::Anita2 ||
       settings1->WHICH==Payload::Anita3 ||
       settings1->WHICH==Payload::Anita4) { /// @todo presumably this is also correct for ANITA-4

      n_eplane = icemc::constants::const_z;
      n_eplane.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
      n_hplane = -icemc::constants::const_y;
      n_hplane.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
      n_normal = icemc::constants::const_x;
      n_normal.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
    }
    else {
      n_eplane = icemc::constants::const_z;
      n_eplane.RotateY(anita1->THETA_ZENITH[ilayer] - icemc::constants::PI/2);
      n_hplane = (-icemc::constants::const_y);
      n_hplane.RotateY(anita1->THETA_ZENITH[ilayer] - icemc::constants::PI/2);
      n_normal = icemc::constants::const_x;
      n_normal.RotateY(anita1->THETA_ZENITH[ilayer] - icemc::constants::PI/2);
    }

    double phi = 0;
    // rotate about z axis for phi
    if (settings1->CYLINDRICALSYMMETRY==1) {
      phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*icemc::constants::PI + anita1->PHI_OFFSET[ilayer]; // + phi_spin;
    }
    else{
      //phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer] +phi_spin;
      phi=anita1->PHI_EACHLAYER[ilayer][ifold];
    }

    n_eplane.RotateZ(phi);
    n_hplane.RotateZ(phi);
    n_normal.RotateZ(phi);
    
    TVector3 positionH;
    TVector3 positionV;    
    for(auto pol : {Seavey::Pol::V, Seavey::Pol::H}){
      int ipol = static_cast<int>(pol);
      TVector3& seaveyPayloadPos = pol == Seavey::Pol::V ? positionV : positionH;

      double phi = 0;
      if (settings1->WHICH == Payload::Anita1 ||
	  settings1->WHICH == Payload::Anita2 ||
	  settings1->WHICH == Payload::Anita3 ||
	  settings1->WHICH == Payload::Anita4 ){
	seaveyPayloadPos = anita1->ANTENNA_POSITION_START[ipol][ilayer][ifold];
      }
      else {
	if (settings1->CYLINDRICALSYMMETRY==1){ // for timing code
	  // phi is 0 for antenna 0 (0-31) and antenna 16 (0-31)
	  // antenna 1 (1-32) and antenna 18 (1-32)
	  phi = (double) ifold / (double) anita1->NRX_PHI[ilayer] * 2 * icemc::constants::PI + anita1->PHI_OFFSET[ilayer];
	}
	else{
	  phi = anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer];
	}
	seaveyPayloadPos = TVector3(anita1->RRX[ilayer]*cos(phi) + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer]),
				    anita1->RRX[ilayer]*sin(phi) + anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer]),
				    anita1->LAYER_VPOSITION[ilayer]);
      }
     
    } 
    fSeaveys.emplace_back(Seavey(positionV, positionH,  n_eplane,  n_hplane, n_normal, settings1));
  }
}




void anitaSim::ANITA::addSignalToRX(const icemc::PropagatingSignal& signal, int rx, int inu){

  int ifold, ilayer;
  getLayerFoldFromTriggerRX(rx, ilayer, ifold);
  
  static bool firstTime = true;
  if(inu == 522 && firstTime){
    for(int i=0; i < fSeaveys.size(); i++){
      // if(true || i==34){
      if(true || i==41){
  	fSeaveys.at(i).setDebug(true);
      }
    }
    if(rx==fSeaveys.size()-1){
      firstTime = false;
    }
  }
  else{
    for(int i=0; i < fSeaveys.size(); i++){
      fSeaveys.at(i).setDebug(false);
    }
  }

  if(rx >= 0 && rx < fSeaveys.size()){
    // actually we add the signal to the new Seavey class
    fSeaveys.at(rx).addSignal(signal);
  }
}



// bool anitaSim::ANITA::applyTrigger(const std::vector<TGraph>& pureSignalVoltageTimeGraphs, const TVector& poyntingVector, const TVector& polarizationVector){
bool anitaSim::ANITA::applyTrigger(int inu){
  
  //////////////////////////////////////
  //       EVALUATE GLOBAL TRIGGER    //
  //          FOR VPOL AND HPOL       //
  //////////////////////////////////////

  fVoltsRX.reset();

  auto globalTrigger = std::make_shared<GlobalTrigger>(fSettings, dynamic_cast<Anita*>(this));

  // int discones_passing = 0;  // number of discones that pass

  for(int pol=0;  pol < NPOL; pol++){
    for(int ant=0; ant < nAnt; ant++){
      fJustNoiseTrig[pol][ant].resize(Anita::HALFNFOUR, 0);
      fJustSignalTrig[pol][ant].resize(Anita::HALFNFOUR, 0);
      fJustNoiseDig[pol][ant].resize(Anita::HALFNFOUR, 0);
      fJustSignalDig[pol][ant].resize(Anita::HALFNFOUR, 0);
    }
  }

  // start looping over antennnas.
  // ilayer loops through vertical layers

  globalTrigger->volts_rx_rfcm_trigger.assign(16,  std::vector <std::vector <double> >(3,  std::vector <double>(0)));

  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  
  for (int antNum=0; antNum < getNumRX(); antNum++) { // loop over layers on the payload
      
      ChanTrigger ct(this);

      // int antNum = this->GetRxTriggerNumbering(ilayer, ifold);
      int ilayer, ifold;
      getLayerFoldFromTriggerRX(antNum, ilayer, ifold);
      // int antNum = this->GetRx(ilayer, ifold);
      ct.readInSeavey(fSettings,  &fSeaveys.at(antNum), antNum, this);

      // this->GetAntennaOrientation(fSettings,  this,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
      // ct.ApplyAntennaGain(fSettings, this, fScreenPtrIDontOwn, antNum, n_eplane, n_hplane, n_normal, inu);

      // std::cout << antNum << "\t" << count_rx << std::endl;

      ct.TriggerPath(fSettings, this, antNum, this);
      ct.DigitizerPath(fSettings, this, antNum);
      ct.TimeShiftAndSignalFluct(fSettings, this, antNum,
				 fVoltsRX.rfcm_lab_e_all.at(antNum).data(),
				 fVoltsRX.rfcm_lab_h_all.at(antNum).data());
      ct.saveTriggerWaveforms(&fJustSignalTrig[0][antNum][0], &fJustSignalTrig[1][antNum][0], &fJustNoiseTrig[0][antNum][0], &fJustNoiseTrig[1][antNum][0]);
      ct.saveDigitizerWaveforms(&fJustSignalDig[0][antNum][0], &fJustSignalDig[1][antNum][0], &fJustNoiseDig[0][antNum][0], &fJustNoiseDig[1][antNum][0]);

      if (fSettings->SCALEDOWNLCPRX1){
	globalTrigger->volts[0][ilayer][0] = globalTrigger->volts[0][ilayer][0]/sqrt(2.);
      }

      if (fSettings->RCPRX2ZERO){
	globalTrigger->volts[1][ilayer][1]=0.;
      }

      if (fSettings->LCPRX2ZERO){
	globalTrigger->volts[0][ilayer][1]=0.;
      }

      if (fSettings->SIGNAL_FLUCT) {
	if (fSettings->WHICH==Payload::AnitaLite) {
	  globalTrigger->volts[ilayer][ifold][0]+=gRandom->Gaus(0., this->VNOISE_ANITALITE[ifold]);
	  globalTrigger->volts[ilayer][ifold][1]+=gRandom->Gaus(0., this->VNOISE_ANITALITE[ifold]);
	} //else
      } //if adding noise

      ct.WhichBandsPass(fSettings, this, globalTrigger.get(), this, ilayer, ifold, fThresholdsAnt[antNum]);

  //   } //loop through the phi-fold antennas
  // }  //loop through the layers of antennas
    }

  int count_pass = 0;
  // globalTrigger->PassesTrigger(fSettings, this, discones_passing, 2, fL3trig, fL2trig, fL1trig, fSettings->antennaclump, loctrig, loctrig_nadironly, inu, thispasses);
  const int triggerMode = 2;
  globalTrigger->PassesTrigger(fSettings, this, triggerMode, fTriggerState);  

  ///////////////////////////////////////
  //       Require that it passes      //
  //            global trigger         //
  ///////////////////////////////////////
  // for Anita-lite,  Anita Hill, just L1 requirement on 2 antennas. This option is currently disabled
  // Save events that generate an RF trigger or that are part of the min bias sample
  // Minimum bias sample: save all events that we could see at the payload
  // Independentely from the fact that they generated an RF trigger

  bool eventPassesTrigger = false;
  if ( (fTriggerState.passes.at(0) > 0 && this->pol_allowed[0]==1)
       || (fTriggerState.passes.at(1) > 0 && this->pol_allowed[1]==1)
       || (fSettings->TRIGTYPE==0 && count_pass>=fSettings->NFOLD)
       || (fSettings->MINBIAS==1)){
    eventPassesTrigger = true;
  }

  return eventPassesTrigger;
}


void anitaSim::ANITA::write(const icemc::Event& event) {
  fAnitaOutput.fillRootifiedAnitaDataTrees(event);
}




double anitaSim::ANITA::GetAverageVoltageFromAntennasHit(const Settings *settings1, int *nchannels_perrx_triggered, const double *voltagearray, double& volts_rx_sum) const {
  double sum=0;
  int count_hitantennas=0;
  for (int i=0;i<settings1->NANTENNAS;i++) {
    if (nchannels_perrx_triggered[i]>=3) {
      sum+=voltagearray[i];
      count_hitantennas++;
    } //if
  } //for
  volts_rx_sum = sum;
  sum = sum/(double)count_hitantennas;
  return sum;
}
//end GetAverageVoltageFromAntennasHit()
