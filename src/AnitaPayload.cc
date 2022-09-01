#include "AnitaPayload.h"
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

anitaSim::AnitaPayload::AnitaPayload(const Settings* settings)
// anitaSim::ANITA::ANITA(const Settings* settings, const RootOutput* ro ) 
  : FlightDataManager(settings), Anita(settings, settings->getOutputDir(), this),
    fSettings(settings),
    fNumRX(settings ? settings->NANTENNAS : 48),
    fVoltsRX(settings ? settings->NANTENNAS : 0),
    fTriggerState(settings),
    fTriggerStateNoise(settings),
    fAnitaOutput(this, settings, settings->getOutputDir(), settings->getRun())
{
  initSeaveys();

  if(fSettings->WHICH == Payload::Anita1Simple ||
     fSettings->WHICH == Payload::Anita1){
    SetDiffraction(); // for the upper ring
  }

  fThresholdsAnt.resize(getNumRX());

  
  // sets position of balloon and related quantities
  // SetDefaultBalloonPosition();
  // SetNoise(fSettings, this, antarctica);  
}


anitaSim::AnitaPayload::~AnitaPayload(){
    
}


bool anitaSim::AnitaPayload::chanceInHell(const icemc::PropagatingSignal& signal){
  ///@todo do something much, much cleverer here...
  /// it's not even clear this is anywhere near the correct value
  //return (signal.maxEField() > signalThreshold());
  return (signal.maxEField() > 1e-18); // This value seemed good? Needs further verification
}

double anitaSim::AnitaPayload::signalThreshold() const {
  //@todo replace rough values with actual implementation
  //double vnoise = anitaSim::VNOISE[0]; //VNOISE[0]
  double maxthreshold = 2.3; //Need to implement
  double heff_max = 0.62639; // max effective height = 0.62639, need to implement
  //double bwmin = anitaSim::bwmin; // Minimum width of bandwidth slice

  double thresh = (Anita::VNOISE[0]/10)*maxthreshold/(heff_max*(Anita::bwmin/1.E6));
  
  return thresh;
}

const Geoid::Position& anitaSim::AnitaPayload::getPosition(double time){
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


TVector3 anitaSim::AnitaPayload::getPositionRX(Int_t rx) const {
  return fSeaveys.at(rx).getPosition(Seavey::Pol::V);
}




void anitaSim::AnitaPayload::initSeaveys() {

  for(int rx = 0; rx < getNumRX(); rx++){
    int ilayer = -1;
    int ifold = -1;
    getLayerFoldFromTriggerRX(rx, ilayer, ifold);
    // std::cout << "Seaveys: rx = " << rx << ", ilayer = " << ilayer << ", ifold = " << ifold << ",  rxTrigger = " << this->GetRxTriggerNumbering(ilayer, ifold) << std::endl;
    

    TVector3 n_eplane;
    TVector3 n_hplane;
    TVector3 n_normal;
    
    if(fSettings->WHICH==Payload::Anita1 ||
       fSettings->WHICH==Payload::Anita2 ||
       fSettings->WHICH==Payload::Anita3 ||
       fSettings->WHICH==Payload::Anita4) { /// @todo presumably this is also correct for ANITA-4

      n_eplane = icemc::constants::const_z;
      n_eplane.RotateY(ANTENNA_DOWN[ilayer][ifold]);
      n_hplane = -icemc::constants::const_y;
      n_hplane.RotateY(ANTENNA_DOWN[ilayer][ifold]);
      n_normal = icemc::constants::const_x;
      n_normal.RotateY(ANTENNA_DOWN[ilayer][ifold]);
    }
    else {
      n_eplane = icemc::constants::const_z;
      n_eplane.RotateY(THETA_ZENITH[ilayer] - icemc::constants::PI/2);
      n_hplane = (-icemc::constants::const_y);
      n_hplane.RotateY(THETA_ZENITH[ilayer] - icemc::constants::PI/2);
      n_normal = icemc::constants::const_x;
      n_normal.RotateY(THETA_ZENITH[ilayer] - icemc::constants::PI/2);
    }

    double phi = 0;
    // rotate about z axis for phi
    if (fSettings->CYLINDRICALSYMMETRY==1) {
      phi=(double)ifold/(double)NRX_PHI[ilayer]*2*icemc::constants::PI + PHI_OFFSET[ilayer]; // + phi_spin;
    }
    else{
      //phi=PHI_EACHLAYER[ilayer][ifold] + PHI_OFFSET[ilayer] +phi_spin;
      phi=PHI_EACHLAYER[ilayer][ifold];
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
      if (fSettings->WHICH == Payload::Anita1 ||
	  fSettings->WHICH == Payload::Anita2 ||
	  fSettings->WHICH == Payload::Anita3 ||
	  fSettings->WHICH == Payload::Anita4 ){
	seaveyPayloadPos = ANTENNA_POSITION_START[ipol][ilayer][ifold];
      }
      else {
	if (fSettings->CYLINDRICALSYMMETRY==1){ // for timing code
	  // phi is 0 for antenna 0 (0-31) and antenna 16 (0-31)
	  // antenna 1 (1-32) and antenna 18 (1-32)
	  phi = (double) ifold / (double) NRX_PHI[ilayer] * 2 * icemc::constants::PI + PHI_OFFSET[ilayer];
	}
	else{
	  phi = PHI_EACHLAYER[ilayer][ifold] + PHI_OFFSET[ilayer];
	}
	seaveyPayloadPos = TVector3(RRX[ilayer]*cos(phi) + LAYER_HPOSITION[ilayer]*cos(LAYER_PHIPOSITION[ilayer]),
				    RRX[ilayer]*sin(phi) + LAYER_HPOSITION[ilayer]*sin(LAYER_PHIPOSITION[ilayer]),
				    LAYER_VPOSITION[ilayer]);
      }
    }
    fSeaveys.emplace_back(Seavey(positionV, positionH,  n_eplane,  n_hplane, n_normal, fSettings));

    if(fSeaveys.size()==37){
      fSeaveys.back().setDebug(true);
    }
  }
}




void anitaSim::AnitaPayload::addSignalToRX(const icemc::PropagatingSignal& signal, int rx, int inu){
  if(rx >= 0 && rx < fSeaveys.size()){
    fSeaveys.at(rx).addSignal(signal);
  }
  else{
    icemc::report() << icemc::severity::error << "rx " << rx << " is out of range (" << 0 << ", " << fSeaveys.size() << ")." << std::endl;
  }
}


bool anitaSim::AnitaPayload::applyTrigger(int inu){
  
  //////////////////////////////////////
  //       EVALUATE GLOBAL TRIGGER    //
  //          FOR VPOL AND HPOL       //
  //////////////////////////////////////

  fVoltsRX.reset(); /// All waveforms downstream of antennas (e.g. signal/ditizger path) set to zero

  auto globalTrigger = std::make_shared<GlobalTrigger>(fSettings, dynamic_cast<Anita*>(this));

  for (int antNum=0; antNum < getNumRX(); antNum++) { // loop over layers on the payload

    ChanTrigger ct(fSettings, this);
    ct.readInSeavey(&fSeaveys.at(antNum), antNum, this);

    ct.TriggerPath(this, antNum);
    ct.DigitizerPath(this, antNum);
    ct.SaveLabradorWaveforms(fVoltsRX.channelsV.at(antNum).rfcm_lab_all.data(),
			     fVoltsRX.channelsH.at(antNum).rfcm_lab_all.data());
    ct.saveTriggerWaveforms(fVoltsRX.channelsV.at(antNum).justSignalTrig.data(),
			    fVoltsRX.channelsH.at(antNum).justSignalTrig.data(),
			    fVoltsRX.channelsV.at(antNum).justNoiseTrig.data(),
			    fVoltsRX.channelsH.at(antNum).justNoiseTrig.data());
    
    ct.saveDigitizerWaveforms(fVoltsRX.channelsV.at(antNum).justSignalDig.data(),
			      fVoltsRX.channelsH.at(antNum).justSignalDig.data(),
			      fVoltsRX.channelsV.at(antNum).justNoiseDig.data(),
			      fVoltsRX.channelsH.at(antNum).justNoiseDig.data());


    int ilayer, ifold;
    getLayerFoldFromTriggerRX(antNum, ilayer, ifold);
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

    ct.WhichBandsPass(this, globalTrigger.get(), this, ilayer, ifold, fThresholdsAnt[antNum]);
  }

  // globalTrigger->PassesTrigger(fSettings, this, discones_passing, 2, fL3trig, fL2trig, fL1trig, fSettings->antennaclump, loctrig, loctrig_nadironly, inu, thispasses);
  const int triggerMode = 2;
  globalTrigger->PassesTrigger(this, triggerMode, fTriggerStateNoise, true);
  globalTrigger->PassesTrigger(this, triggerMode, fTriggerState);
  
  ///////////////////////////////////////
  //       Require that it passes      //
  //            global trigger         //
  ///////////////////////////////////////
  // for Anita-lite,  Anita Hill, just L1 requirement on 2 antennas. This option is currently disabled
  // Save events that generate an RF trigger or that are part of the min bias sample
  // Minimum bias sample: save all events that we could see at the payload
  // Independentely from the fact that they generated an RF trigger

  bool eventPassesTrigger = false;
  if ( (fTriggerState.passes.at(0) > 0 && fTriggerStateNoise.passes.at(0) == 0 && this->pol_allowed[0]==1)
       || (fTriggerState.passes.at(1) > 0 && fTriggerStateNoise.passes.at(1) == 0 && this->pol_allowed[1]==1)
       || (fSettings->MINBIAS==1)){
    eventPassesTrigger = true;
  }
  else if (fTriggerStateNoise.passes.at(0) > 0 || fTriggerStateNoise.passes.at(1) > 0)
    std::cout << "~~~ Event " << inu << " triggered on noise!\n";
  
  return eventPassesTrigger;
}


void anitaSim::AnitaPayload::write(const icemc::Event& event) {
  fAnitaOutput.fillRootifiedAnitaDataTrees(event);
}




