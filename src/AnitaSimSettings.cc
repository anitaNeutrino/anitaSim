#include "AnitaSimSettings.h"
#include "Report.h"
#include "Constants.h"

std::ostream& operator<<(std::ostream& os, const anitaSim::Payload& which){
  switch(which){
  case anitaSim::Payload::AnitaLite:
    return os << "Payload::AnitaLite";
  case anitaSim::Payload::Ross:
    return os << "Payload::Ross";
  case anitaSim::Payload::Anita1Simple:
    return os << "Payload::Anita1Simple";
  case anitaSim::Payload::Custom:
    return os << "Payload::Custom";
  case anitaSim::Payload::AnitaHill:
    return os << "Payload::AnitaHill";
  case anitaSim::Payload::SLAC:
    return os << "Payload::SLAC";
  case anitaSim::Payload::Anita1:
    return os << "Payload::Anita1";
  case anitaSim::Payload::EeVEX:
    return os << "Payload::EeVEX";
  case anitaSim::Payload::Anita2:
    return os << "Payload::Anita2";
  case anitaSim::Payload::Anita3:
    return os << "Payload::Anita3";
  case anitaSim::Payload::Anita4:
    return os << "Payload::Anita4";
  case anitaSim::Payload::Satellite:
    return os << "Payload::Satellite";
  default:
    return os << "Unknown Payload";
  }
}


anitaSim::Settings::Settings(): icemc::Settings() {
  Initialize();
}


void anitaSim::Settings::Initialize(){
  NDISCONES_PASS=3;
  FREQ_LOW_SEAVEYS=200.E6;
  FREQ_HIGH_SEAVEYS=1200.E6;
  BW_SEAVEYS=FREQ_HIGH_SEAVEYS-FREQ_LOW_SEAVEYS;

  TUFFSON=0;
  ADDCW=0;
}

void anitaSim::Settings::ReadInputs(const char* fileName , std::ofstream &foutput){
  icemc::Settings::ReadInputs(fileName, foutput);


  int whichInt = 0;
  getSetting("Which payload", whichInt);
  WHICH = static_cast<anitaSim::Payload>(whichInt);
  getSetting("Antenna layers", NLAYERS);

  switch(WHICH){
  case anitaSim::Payload::Anita3:
    ANITAVERSION = 3;
    break;
  case anitaSim::Payload::Anita4:
    ANITAVERSION = 4;
    break;
  default:
    ANITAVERSION = 0;
    break;
  }

  setNrxPhiAndNantennasFromWhich();
  
  if(((WHICH==anitaSim::Payload::Ross || WHICH==anitaSim::Payload::Anita1) && NLAYERS!=4) || (WHICH==anitaSim::Payload::AnitaLite && NLAYERS!=1) || (WHICH==anitaSim::Payload::EeVEX && NLAYERS!=1)){
    std::cout << "Non-default setting: WHICH = " << WHICH << " and NLAYERS= " << NLAYERS << std::endl;
  }  

  getSetting("Inclination top three layers", INCLINE_TOPTHREE);

  if(INCLINE_TOPTHREE!=10){
    std::cout << "Non-default setting: INCLINE_TOPTHREE= " << INCLINE_TOPTHREE << std::endl;
  }

  getSetting("Inclination fourth layer", INCLINE_NADIR);

  if(INCLINE_NADIR!=10){
    std::cout << "Non-default setting: INCLINE_NADIR= " << INCLINE_NADIR << std::endl;
  }

  int whichPathInt = 0;
  getSetting("Flight path", whichPathInt);
  WHICHPATH = static_cast<FlightPath>(whichPathInt);

  if((WHICH==anitaSim::Payload::AnitaLite    && WHICHPATH!=FlightPath::AnitaLite) ||
     (WHICH==anitaSim::Payload::Anita1Simple && WHICHPATH!=FlightPath::Anita1) || 
     (WHICH==anitaSim::Payload::Anita1       && WHICHPATH!=FlightPath::Anita1) ||
     (WHICH==anitaSim::Payload::Anita2       && WHICHPATH!=FlightPath::Anita2) ||
     (WHICH==anitaSim::Payload::Anita3       && WHICHPATH!=FlightPath::Anita3) ||
     (WHICH==anitaSim::Payload::Anita4       && WHICHPATH!=FlightPath::Anita4)){
    std::cout << "Non-default setting: WHICHPATH = " << WHICHPATH
	      << " and WHICH = " << WHICH << std::endl;
  }

  getSetting("Balloon latitude", BN_LATITUDE);

  getSetting("Balloon longitude", BN_LATITUDE);

  if((BN_LONGITUDE!=999 || BN_LATITUDE!=999) && WHICHPATH==FlightPath::FixedPosition){
    std::cout << "BN_LATITUDE: "<< BN_LATITUDE << ", BN_LONGITUDE: " << BN_LONGITUDE << std::endl;
  }

  if (BN_LONGITUDE>180. && BN_LONGITUDE!=999){
    std::cout << "Entered balloon longitude wrong!  Should be between -180 and 180 degrees." << std::endl;
  }
  getSetting("Balloon orientation", RANDOMIZE_BN_ORIENTATION);


  if (RANDOMIZE_BN_ORIENTATION==1 &&
      (WHICHPATH==FlightPath::AnitaLite ||
       WHICHPATH==FlightPath::Anita1 ||
       WHICHPATH==FlightPath::Anita2 ||
       WHICHPATH==FlightPath::Anita3 ||
       WHICHPATH==FlightPath::Anita4)){
    std::cout << "Warning:: Strangely you asked for a real flight path but a randomized balloon orientation.  WILL BE OVERRIDDEN." << std::endl;
  }

  getSetting("Balloon altitude", BN_ALTITUDE);


  // whether to use constant gains as entered in GetBeamWidths (0) or to use Ped's measurements as entered in ReadGains (1)
  // GAINS is actually an int, not a double...
  getSetting("Gain setting", GAINS);



  getSetting("Trigger scheme", TRIGGERSCHEME);

  // Ben S, leaving this here... should it get its own config file entry?
  TRIGTYPE=1; // ANITA, not anita-lite.  But the anita-lite code back in later

  getSetting("Band thresholds", tempThresholds);

  getSetting("Banding", BANDING);


  if(BANDING !=0 && BANDING!= 1 && BANDING!=2 && BANDING!=4) {
    std::cout << "Banding should be set to 0 (Anita 1), 1 (custum), 2 (Anita 2), 3 (Satellite) or 4 (Anita 3)."
	      << std::endl;
    exit(1);
  }

  if ((TRIGGERSCHEME==0 || TRIGGERSCHEME==1) && BANDING!=1) {
    std::cout << "Frequency domain trigger schemes can only be used with user-set sub-bands." << std::endl;
    exit(1);
  }

  if (TRIGGERSCHEME==2 && BANDING==1) {
    std::cout << "Time domain trigger scheme only works with Anita 1, Anita 2 or Anita 3 banding data, you can't set your own bands." << std::endl;
    exit(1);
  }


  getSetting("Lower band edges", bandLowEdgesMHz);
  getSetting("Upper band edges", bandHighEdgesMHz);
  getSetting("Required bands", requiredBands);
  getSetting("Allowed bands", allowedBands);
  getSetting("Number of bands", NBANDS);
  getSetting("Percent bandwidth", PERCENTBW);

  getSetting("Notch filter limits", notchFilterLimitsMHz);

  getSetting("Num antenna channels for L1 trigger", trigRequirements[0]);
  getSetting("Num L1 hits to pass L2", trigRequirements[1]);
  getSetting("Num antenna for L2 trigger", antennaclump);
  getSetting("Require centre antenna", REQUIRE_CENTRE);
  getSetting("L3 trigger requirement", trigRequirements[2]);
  getSetting("LCP/RCP or V/H", LCPRCP);
  getSetting("Channels required polarization", channelRequirePol);
  getSetting("Channels allowed polarization", channelAllowedPol);

  getSetting("Exta antennas in trigger", DISCONES);
  getSetting("Nadir only trigger", INCLUDE_NADIRONLY);

  if (INCLUDE_NADIRONLY!=0){
    std::cout << "Non-default setting:  INCLUDE_NADIRONLY = " << INCLUDE_NADIRONLY << std::endl;
  }

  getSetting("ANITA-1 channel masking", CHMASKING);

  if (WHICHPATH!=FlightPath::Anita1 && CHMASKING==1) {
    std::cout << "Cannot include masking for flights other than the ANITA-1 flight." << std::endl;
    std::cout << "For the ANITA-3 channel masking, it is implemented together with the phi masking and it's turned on whenever the PHIMASKING is ON." << std::endl;
    std::cout << "CHMASKING set to 0." << std::endl;
    CHMASKING=0;
  }

  getSetting("ANITA-2 channel masking", PHIMASKING);



  getSetting("Scale down LCP voltage 1st ant", SCALEDOWNLCPRX1);
  getSetting("Scale down E pol 1st ant", SCALEDOWNEPOLRX1);
  getSetting("Scale down H pol 1st ant", SCALEDOWNHPOLRX1);
  getSetting("Scale down E pol 2nd ant", SCALEDOWNEPOLRX2);
  getSetting("E pol scale down factor 2nd ant", SCALEFACTOREPOLRX2);
  getSetting("H pol scale down factor 2nd ant", SCALEDOWNHPOLRX2);
  getSetting("E pol 2nd ant dead", EPOLRX2ZERO);
  getSetting("H pol 2nd ant dead", HPOLRX2ZERO);
  getSetting("RCP 2nd ant dead", RCPRX2ZERO);
  getSetting("LCP 2nd ant dead", LCPRX2ZERO);

  if (WHICH==anitaSim::Payload::AnitaLite && !(SCALEDOWNEPOLRX1==1 && RCPRX2ZERO==1)){
    std::cout << "Non-default setting:  WHICH= " << WHICH << " and EPOLRX2ZERO= " << EPOLRX2ZERO << std::endl;
  }



  getSetting("Loop over boresights", BORESIGHTS);
  if (BORESIGHTS==0){
    std::cout << "Warning!  Non-standard parameter setting.  BORESIGHTS = " << BORESIGHTS << std::endl;
  }
  getSetting("Max interaction distance", MAXHORIZON);
  
  getSetting("Theta resolution", SIGMA_THETA);
  if (SIGMA_THETA==1){
    std::cout << "Non-default setting:  SIGMA_THETA = 1" << std::endl;
  }
  SIGMA_THETA*=TMath::DegToRad(); // immediately convert degrees to radians
  getSetting("Low frequency", FREQ_LOW);

  // if (FREQ_LOW_SEAVEYS>anita1->FREQ_LOW){
  //   FREQ_LOW_SEAVEYS=anita1->FREQ_LOW;
  // }
  getSetting("High frequency", FREQ_HIGH);
  BW = FREQ_HIGH - FREQ_LOW; // total bandwidth of simulation






  getSetting("SLAC run", SLAC);
  // get positions of the anita payload during the slac test
  if (SLAC){
    icemc::report() << icemc::severity::error << "SLAC is no longer supported!" << std::endl;
    exit(1);
  }




  // this rotates surface slope 10 degrees away from south pole
  if (SLAC==1){
    std::cout << "Warning!  Non-standard parameter setting.  SLAC = " << SLAC << std::endl;
  }
  if (SLAC) {
    foutput << "!!!!SLAC setting causes some settings to be superseded:" << std::endl;
    FIRN=0; // no firn
    foutput << "FIRN=0" << std::endl;
    SLOPEY=0; // slopeyness off
    foutput << "SLOPEY=0" << std::endl;
    BORESIGHTS=1; // loop over boresights
    foutput << "BORESIGHTS=1" << std::endl;
    BN_ALTITUDE=4.22/0.3; // balloon altitude in ft.!!
    foutput << "BN_ALTITUDE=4.22/0.3" << std::endl;
    RANDOMIZE_BN_ORIENTATION=0; // don't randomize the balloon orientation
    foutput << "RANDOMIZE_BN_ORIENTATION=0" << std::endl;
    SKIPCUTS=1; // don't make chance in hell cuts
    foutput << "SKIPCUTS=1" << std::endl;
    SLACSLOPE=5.8; // slope of the ice in degrees
    foutput << "SLACSLOPE=5.8" << std::endl;
    SLACICELENGTH=5.02; // length of the block of ice
    foutput << "SLACICELENGTH=5.02" << std::endl;
  }
  getSetting("SLAC horizontal distance", SLAC_HORIZDIST);
  getSetting("SLAC ice slope", SLACSLOPE);
  getSetting("SLAC block length", SLACICELENGTH);
  getSetting("SLAC interaction depth", SLAC_HORIZ_DEPTH);
  SLAC_DEPTH=tan(SLACSLOPE*icemc::constants::RADDEG)*(SLACICELENGTH-SLAC_HORIZ_DEPTH) // height from lowest point of ice
    +21.375*icemc::constants::CMINCH/100.; // height from beam to lowest point of ice








  getSetting("Coherent power threshold", COHERENT_THRESHOLD );



  // default values are 0
  APPLYIMPULSERESPONSEDIGITIZER=0;
  APPLYIMPULSERESPONSETRIGGER=0;
  USETIMEDEPENDENTTHRESHOLDS=0;
  USEDEADTIME=0;
  getSetting("Digitizer path impulse response", APPLYIMPULSERESPONSEDIGITIZER);
  std::cout << "Apply impulse response to digitizer path: " << APPLYIMPULSERESPONSEDIGITIZER << std::endl;
  getSetting("Trigger path impulse response", APPLYIMPULSERESPONSETRIGGER);
  std::cout << "Apply impulse response to trigger path: " << APPLYIMPULSERESPONSETRIGGER << std::endl;

#ifdef ANITA_UTIL_EXISTS
  if ( (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER) && WHICH!=anitaSim::Payload::Anita2 && WHICH!=anitaSim::Payload::Anita3 && WHICH!=anitaSim::Payload::Anita4) {
    std::cout << "Signal chain impulse response is only available for " << anitaSim::Payload::Anita2  << ", " << anitaSim::Payload::Anita3 << ", and " << anitaSim::Payload::Anita4 << std::endl;
    exit(1);
  }
#endif
#ifndef ANITA_UTIL_EXISTS
  if (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER){
    std::cout << "Signal chain impulse response can only be applied when the Anita tools are sourced." << std::endl;
    exit(1);
  }
#endif
  getSetting("Time dependent thresholds", USETIMEDEPENDENTTHRESHOLDS);
  std::cout << "Use time-dependent thresholds: " << USETIMEDEPENDENTTHRESHOLDS << std::endl;
  getSetting("Dead time", USEDEADTIME);
  std::cout << "Use dead time from flight: " << USEDEADTIME << std::endl;
  
  if ( (USETIMEDEPENDENTTHRESHOLDS || USEDEADTIME) && WHICH!=anitaSim::Payload::Anita3 && WHICH!=anitaSim::Payload::Anita4) {
    std::cout << "Time-dependent thresholds are only available for " << anitaSim::Payload::Anita3 << " and " << anitaSim::Payload::Anita4 << std::endl;
    exit(1);
  }


  getSetting("Digitizer noise from flight", NOISEFROMFLIGHTDIGITIZER);
  std::cout << "Use noise from flight for digitizer path: " << NOISEFROMFLIGHTDIGITIZER << std::endl;

  getSetting("Trigger noise from flight", NOISEFROMFLIGHTTRIGGER);
  std::cout << "Use noise from flight for trigger path: " << NOISEFROMFLIGHTTRIGGER << std::endl;

#ifdef ANITA3_EVENTREADER
  if ( (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER) && (WHICH!=anitaSim::Payload::Anita3 && WHICH!=anitaSim::Payload::Anita4)) {
    std::cout << "Noise from flight only available for " << anitaSim::Payload::Anita3 << " and " << anitaSim::Payload::Anita4 << std::endl;
    exit(1);
  }
  if (!APPLYIMPULSERESPONSETRIGGER && NOISEFROMFLIGHTTRIGGER ){
    std::cout << "Noise from flight can only be applied to trigger path if impulse reponse is also used " << std::endl;
    exit(1);
  }
#endif

#ifndef ANITA_UTIL_EXISTS
  if (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER){
    std::cout << "Noise from flight can only be applied when the Anita tools are sourced." << std::endl;
    exit(1);
  }
#endif




  getSetting("Efficiency scan",      TRIGGEREFFSCAN                       );
  if (TRIGGEREFFSCAN==1){
    getSetting("Central phi-sector",   trigEffScanPhi                     );
    getSetting("Apply pulse at surf",  TRIGGEREFFSCAPULSE                 );
    getSetting("Off-axis attenuation", efficiencyScanOffAxisAttenuations  );
    getSetting("Rings used",           efficiencyScanRingsUsed            );
    getSetting("Phi-sectors delays",   efficiencyScanPhiSectorDelay       );
    getSetting("Ring delays",          efficiencyScanRingDelay            );
    getSetting("Ring delays to phi",   efficiencyScanApplyRingDelay       );
  }
  
  getSetting("Simulate TUFFs", TUFFSON);
  getSetting("Which TUFFs are on", whichTUFFsON);
  if (TUFFSON) std::cout << "TUFFs are simulated " << std::endl;
  
  getSetting("Add CW", ADDCW);
  if(ADDCW) std::cout << "Adding CW " << std::endl;




  
};




void anitaSim::Settings::setNrxPhiAndNantennasFromWhich(){
  if (WHICH==anitaSim::Payload::AnitaLite) { // anita-lite
    NFOLD=3;
    CYLINDRICALSYMMETRY=0;
    NRX_PHI[0]=2;
  }
  else if (WHICH==anitaSim::Payload::Ross) {
    CYLINDRICALSYMMETRY=1;
    NRX_PHI[0]=5;
    NRX_PHI[1]=5;
    NRX_PHI[2]=5;
    NRX_PHI[3]=5;
    NRX_PHI[4]=4;
  }
  else if (WHICH==anitaSim::Payload::Anita1Simple) {
    CYLINDRICALSYMMETRY=1;
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
    NRX_PHI[3]=8;
  }
  else if (WHICH==anitaSim::Payload::Custom) {
    std::cout << "Is this configuration cylindrically symmetric? Yes(1) or No(0)\n";
    std::cin >> CYLINDRICALSYMMETRY;
    std::cout << "How many layers?\n";
    std::cin >> NLAYERS;
    for (int i=0;i<NLAYERS;i++) {
      std::cout << "How many antennas in the " << i << "th layer?\n";
      std::cin >> NRX_PHI[i];
    }
    std::cout << "How many polarizations must pass a voltage threshold?\n";
    std::cin >> NFOLD;

  } //else if (custom payload)
  else if (WHICH==anitaSim::Payload::AnitaHill) {// anita hill
    if (NLAYERS!=2){
      std::cout << "Warning!!! Did not enter the right number of layers in the input file.  For Anita Hill, it's 2." << std::endl;
    }
    CYLINDRICALSYMMETRY=1;
    NRX_PHI[0]=1; // this is how many antennas we have in phi on each "layer"
    NRX_PHI[1]=1; // for anita hill, we are calling each station a different "layer"
  }
  else if(WHICH==anitaSim::Payload::Anita1) { // Kurt's measurements for the first flight in elog 345
    //NFOLD=8;
    CYLINDRICALSYMMETRY=0;
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
  }
  else if (WHICH==anitaSim::Payload::EeVEX) {
    CYLINDRICALSYMMETRY=1;
    NRX_PHI[0]=360;
    NRX_PHI[1]=360;
    NRX_PHI[2]=360;
    NRX_PHI[3]=360;
    NRX_PHI[4]=360;
  }
  else if (WHICH==anitaSim::Payload::Anita2) {
    CYLINDRICALSYMMETRY=0;
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
    NRX_PHI[3]=8;
  }
  else if (WHICH==anitaSim::Payload::Anita3 || WHICH==anitaSim::Payload::Anita4) { // ANITA-3 and ANITA-4
    CYLINDRICALSYMMETRY=0;
    //these are physical layers
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
    NRX_PHI[3]=16;
  }
  else if (WHICH==anitaSim::Payload::Satellite) { // satellite
    CYLINDRICALSYMMETRY=1;
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;

  } //else if (satellite)

  NANTENNAS=0;
  for (int i=0;i<NLAYERS;i++){
    NANTENNAS+=NRX_PHI[i];
  }
}




void anitaSim::Settings::ApplyInputs(Anita* anita1) const {
  
   //When you look at the Anita payload there are 4 layers, with 8,8,16 and 8 antennas each.  But in the trigger, the top two become one layer of 16 antennas. 
  if (WHICH==Payload::Anita1Simple || WHICH==Payload::Anita1 || WHICH==Payload::Anita2 || WHICH==Payload::Anita3 || WHICH==Payload::Anita4){
    anita1->NTRIGGERLAYERS = NLAYERS - 1;
  }
  else{
    anita1->NTRIGGERLAYERS=NLAYERS;
  }


  anita1->INCLINE_TOPTHREE=INCLINE_TOPTHREE;
  anita1->INCLINE_NADIR=INCLINE_NADIR;

  if(WHICHPATH==FlightPath::AnitaLite){
    anita1->LIVETIME=45.*24.*3600.*0.75; // 45 days for anita-lite
  }
  else if (WHICHPATH==FlightPath::FixedPosition){
    anita1->LIVETIME=6.02*24.*3600.; // anita-lite
  } else if (WHICHPATH==FlightPath::Anita1){
    // kim's livetime for anita
    anita1->LIVETIME=17.*24.*3600.; // for anita, take 34.78 days * 0.75 efficiency
  }
  else if (WHICHPATH==FlightPath::Anita2){
    anita1->LIVETIME=28.5*24*3600;  // Anita-2 livetime taken from paper
  }
  else if (WHICHPATH==FlightPath::Anita3){
    anita1->LIVETIME=17.4*24*3600;  // Anita-3 livetime taken from Ben Strutt's thesis (elog note 698)
  }
  else {
    ///@todo warn here?
    anita1->LIVETIME=14.*24.*3600.; // otherwise use 2 weeks by default
  }
  
  if (WHICH==Payload::EeVEX){
    // EeVEX
    anita1->LIVETIME=100.*24.*3600.; // ultra-long duration balloon flight of 100 days
  }

  // bn1->BN_LATITUDE              = BN_LATITUDE;
  // bn1->BN_LONGITUDE             = BN_LONGITUDE;
  // bn1->BN_ALTITUDE              = BN_ALTITUDE;
  // bn1->RANDOMIZE_BN_ORIENTATION = RANDOMIZE_BN_ORIENTATION;
  // bn1->MAXHORIZON               = MAXHORIZON;

  anita1->GAINS   = GAINS;
  anita1->BANDING = BANDING;
  anita1->NBANDS  = NBANDS;
  anita1->PERCENTBW = PERCENTBW;
  anita1->SIGMA_THETA = SIGMA_THETA;
  anita1->FREQ_LOW    = FREQ_LOW;
  anita1->FREQ_HIGH   = FREQ_HIGH;

  for (unsigned int i=0;i<tempThresholds.size();i++) {
    anita1->bwslice_thresholds.at(i) = tempThresholds.at(i);
  }


  for (unsigned int i=0; i < bandLowEdgesMHz.size(); i++) {
    anita1->bwslice_min[i] = 1e6*bandLowEdgesMHz.at(i);
    anita1->bwslice_max[i] = 1e6*bandHighEdgesMHz.at(i);
    anita1->bwslice_center[i] = 0.5*(anita1->bwslice_min[i] + anita1->bwslice_max[i]);
    anita1->bwslice_width[i] = anita1->bwslice_max[i] - anita1->bwslice_min[i];
  }

for(unsigned int i=0; i < requiredBands.size(); i++){
    anita1->bwslice_required[i] = requiredBands.at(i);
  }

  for(unsigned int i=0; i < allowedBands.size(); i++){
    anita1->bwslice_allowed[i] = allowedBands.at(i);
  }

  anita1->bwmin=1.E10;
  if (anita1->BANDING!=1){
    anita1->bwmin=200.E6;
  }

  const int numBands = 5;
  for (int i=0; i<numBands; i++) {
    if (anita1->BANDING==1) {
      if ((anita1->bwslice_max[i] - anita1->bwslice_min[i]) < anita1->bwmin && anita1->bwslice_allowed[i] == 1){

	anita1->bwmin = anita1->bwslice_max[i] - anita1->bwslice_min[i];
      }
    }
  }


  anita1->NOTCH_MIN = 1e6*notchFilterLimitsMHz.at(0);
  anita1->NOTCH_MAX = 1e6*notchFilterLimitsMHz.at(1);
  
  if (anita1->NOTCH_MIN>anita1->NOTCH_MAX) {
    std::cout << "Min of notch filter is greater than max. Try again." << std::endl;
  }
  if (anita1->NOTCH_MIN!=0 || anita1->NOTCH_MAX!=0){
    std::cout << "Applying a notch filter from " << anita1->NOTCH_MIN << " Hz to "
	      << anita1->NOTCH_MAX << " Hz" << std::endl;
  }	      
  anita1->trigRequirements[0] = trigRequirements[0];
  anita1->trigRequirements[1] = trigRequirements[1];	
  anita1->trigRequirements[2] = trigRequirements[2];	
  anita1->REQUIRE_CENTRE = REQUIRE_CENTRE;		
 for(unsigned int i=0; i < channelRequirePol.size(); i++){
    anita1->pol_required[i] = channelRequirePol.at(i); 
  }
 for(unsigned int i=0; i < channelAllowedPol.size(); i++){
    anita1->pol_allowed[i] = channelAllowedPol.at(i);
  }



  // askFreqGen->SetLPM(useLPM);
  // if (askFreqGen->GetLPM()!=1){
  //   std::cout << "Non-default setting:  LPM= " << askFreqGen->GetLPM() << std::endl;
  // }
  // askFreqGen->SetParameterization(askaryanParameterization);
  // askFreqGen->SetJaime_Factor(jamieFactor);
  // askFreqGen->SetMedium(medium);


  // sec1->SECONDARIES=SECONDARIES;
  // sec1->TAUDECAY=TAUDECAY;

  anita1->trigEffScanPhi = trigEffScanPhi;
  for (unsigned int i=0; i < efficiencyScanOffAxisAttenuations.size(); i++){
    anita1->trigEffScanAtt[i] = efficiencyScanOffAxisAttenuations.at(i);
  }

  for (unsigned int i=0; i < efficiencyScanPhiSectorDelay.size(); i++){
    //convert ns to s
    anita1->trigEffScanPhiDelay[i] = efficiencyScanPhiSectorDelay.at(i)*1e-9;
  }

  for (unsigned int i=0; i < efficiencyScanRingDelay.size(); i++){
    //convert ns to s
    anita1->trigEffScanRingDelay[i] = efficiencyScanRingDelay.at(i)*1e-9;
  }

  for (unsigned int i=0; i < efficiencyScanApplyRingDelay.size(); i++){
    anita1->trigEffScanApplyRingDelay[i] = efficiencyScanApplyRingDelay.at(i);
  }

  for (unsigned int i=0; i < efficiencyScanRingsUsed.size(); i++){
    anita1->trigEffScanRingsUsed[i] = efficiencyScanRingsUsed.at(i);
  }

  
  
  if (TRIGGEREFFSCAN){
    std::cout << "Let's do a trigger efficiency scan!" << std::endl;
    std::cout << "Apply pulse at AMPA (0) or SURF (1) : " << TRIGGEREFFSCAPULSE << std::endl;
    std::cout << "Central phi sector is " << anita1->trigEffScanPhi << std::endl;
    std::cout << "Attenuations are ";
    for (int i=0;i<5;i++) std::cout << anita1->trigEffScanAtt[i] << " ";
    std::cout << std::endl;
    std::cout << "Phi sector delays are ";
    for (int i=0;i<5;i++) std::cout << anita1->trigEffScanPhiDelay[i] << " ";
    std::cout << std::endl;
    std::cout << "The rings used in this scan are ";
    for (int i=0;i<3;i++) std::cout << anita1->trigEffScanRingsUsed[i] << " ";
    std::cout << std::endl;
    std::cout << "Ring delays are applie to : ";
    for (int i=0;i<5;i++) std::cout << anita1->trigEffScanApplyRingDelay[i] << " ";
    std::cout << std::endl;
    std::cout << "Ring delays are for T-M, M-B, T-B : ";
    for (int i=0;i<3;i++) std::cout << anita1->trigEffScanRingDelay[i] << " ";
    std::cout << std::endl;
  }
  


  for (unsigned int i=0; i< whichTUFFsON.size(); i++){
    anita1->TUFFstatus[i] = whichTUFFsON.at(i);
  }

  if (TUFFSON){
    std::cout << "The TUFFs are ON for the whole flight!" << std::endl;
    std::cout << "Notch 0 status " << anita1->TUFFstatus[0] << std::endl;
    std::cout << "Notch 1 status " << anita1->TUFFstatus[1] << std::endl;
    std::cout << "Notch 2 status " << anita1->TUFFstatus[2] << std::endl;
  }


}
