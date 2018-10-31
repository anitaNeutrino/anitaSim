#include "PayloadGeometry.h"
#include "AnitaSimSettings.h"
#include "Constants.h"
#include "EnvironmentVariable.h"
#include "Report.h"


#ifdef ANITA_UTIL_EXISTS
#include "FFTtools.h"
#include "AnitaEventCalibrator.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#endif

anitaSim::PayloadGeometry::PayloadGeometry(const Settings* settings) : fSettings(settings){
  INCLINE_TOPTHREE = 10.; // cant angle of top three layers of antennas
  INCLINE_NADIR = 10.; // cant angle of nadir (bottom) layer of antennas
  
  GetPayload();
}

int anitaSim::PayloadGeometry::GetRxTriggerNumbering(int ilayer, int ifold) const { // get antenna number based on which layer and position it is
  // make the top trigger layer count 1-16 left to right
  if (ilayer==0){
    //cout << "ilayer, ifold, getrx are " << ilayer << "\t" << ifold << "\t" << 2*ifold+ilayer << "\n";
    return 2*ifold;
  }
  else if(ilayer==1) {
    return 2*ifold+1;
  }
  else {
    //cout << "ilayer, ifold, getrx are " << ilayer << "\t" << ifold << "\t" << GetRx(ilayer,ifold) << "\n";
    // return GetRx(ilayer,ifold);

    int irx=0;
    for (int i=0;i<ilayer;i++) {
      irx+=NRX_PHI[i];
    }
    irx+=ifold;
    return irx;
    
  }
}


void anitaSim::PayloadGeometry::getLayerFoldFromTriggerRX(int rx, int& ilayer, int& ifold) const {
  int antNum = rx; 
  ///@todo Do something smarter, this is NOT how to do things...
  ilayer = -1;
  ifold = -1;
  for (int ilayerTemp=0 ;ilayerTemp < fSettings->NLAYERS; ilayerTemp++) { // loop over layers on the payload
    for (int ifoldTemp=0;ifoldTemp<this->NRX_PHI[ilayerTemp];ifoldTemp++) { // ifold loops over phi
      Int_t antNum2 = this->GetRxTriggerNumbering(ilayerTemp, ifoldTemp);
      if(antNum==antNum2){
	ilayer = ilayerTemp;
	ifold = ifoldTemp;
	break;
      }
    }
    if(ilayer > -1){
      break;
    }
  }
}




void anitaSim::PayloadGeometry::GetPayload(){

  std::string ICEMC_DATA_DIR(icemc::EnvironmentVariable::ICEMC_SRC_DIR());
  ICEMC_DATA_DIR += "/data";
  
  // anita-lite payload
  // see comments next to variable definitions
  ///@todo restore me
  // double temp_eachrx[anitaSim::Anita::NPHI_MAX]; // temperature of each antenna (for the anita-lite configuration)
    
  const double gps_offset = atan2(-0.7042,0.71), MINCH = 0.0254, phase_center = 0.17;
  // const double phase_center_anita2=0.17;
  const double phase_center_anita2_analysis=.2;
  //const double gps_offset_anita2=atan2(0.89,-0.29);
  const double gps_offset_anita2=atan2(-0.7085,0.7056); // from elog 473
  const double gps_offset_anita3= 45*icemc::constants::RADDEG; // Linda: 45 degrees from EventReader
  const double phase_center_anita3=0.20; // Linda: phase-centers are around 20 cm inwards of antennas face-end

  number_all_antennas = fSettings->NANTENNAS;
  
  for(int layer=0; layer < NLAYERS_MAX; layer++){
    NRX_PHI[layer] = fSettings->NRX_PHI[layer];
  }

  if (fSettings->WHICH==Payload::AnitaLite) { // anita-lite
		
    // fSettings->NFOLD=3;
		
    // fSettings->CYLINDRICALSYMMETRY=0;
    PHI_EACHLAYER[0][0]=0.;
    PHI_EACHLAYER[0][1]=22.5*icemc::constants::RADDEG;
    // NRX_PHI[0]=2;
    PHI_OFFSET[0]=0;
    THETA_ZENITH[0] = icemc::constants::PI/2+10.*icemc::constants::RADDEG;
    LAYER_VPOSITION[0]=0.; // vertical separation between layers
    LAYER_HPOSITION[0]=0.;  // position of layer relative to center axis in horizontal plane
    LAYER_PHIPOSITION[0]=0.; // phi of position of layer
		
    for (int i=0;i<5;i++) {
      RRX[i]=3.006;
    }

    ///@todo re-enable VNOISE_ANITALITE
    // for (int i=0;i<NRX_PHI[0];i++) {
    //   VNOISE_ANITALITE[i] = ChanTrigger::GetNoise(fSettings,
    // 						  bn1->getAltitude(),
    // 						  //@todo fix this, replaced surfaceUnderBalloon with 0
    // 						  // bn1->getSurfaceUnderBalloon(),
    // 						  0,
						  
    // 						  THETA_ZENITH[i],fSettings->BW_SEAVEYS,temp_eachrx[i]);
    // }
  } //if (ANITA-lite)
    
    //Ross Payload
  else if (fSettings->WHICH==Payload::Ross) {
    //fSettings->NFOLD=8;
    // fSettings->CYLINDRICALSYMMETRY=1;
		
    // NRX_PHI[0]=5;
    // NRX_PHI[1]=5;
    // NRX_PHI[2]=5;
    // NRX_PHI[3]=5;
    // NRX_PHI[4]=4;
		
    PHI_OFFSET[0]=0;
    PHI_OFFSET[1]=2*icemc::constants::PI/(double)NRX_PHI[1]/2;
    PHI_OFFSET[2]=0;
    PHI_OFFSET[3]=2*icemc::constants::PI/(double)NRX_PHI[3]/2;
    PHI_OFFSET[4]=0;
		
    THETA_ZENITH[0]=icemc::constants::PI/2+10.*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+10.*icemc::constants::RADDEG;
    THETA_ZENITH[2]=icemc::constants::PI/2+10.*icemc::constants::RADDEG;
    THETA_ZENITH[3]=icemc::constants::PI/2+10.*icemc::constants::RADDEG;
    THETA_ZENITH[4]=3*icemc::constants::PI/4;
		
    // anita proposal "says that the separation between upper and lower
    // 2 layers of antennas is just under 4m.
    LAYER_VPOSITION[0]=3.5;
    LAYER_VPOSITION[1]=2.0;
    LAYER_VPOSITION[2]=-0.5;
    LAYER_VPOSITION[3]=-2.0;
    LAYER_VPOSITION[4]=-3.5;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1]=0.;
    LAYER_HPOSITION[2]=0.;
    LAYER_HPOSITION[3]=0.;
    LAYER_HPOSITION[4]=0.;
		
    LAYER_PHIPOSITION[0]=0.;
    LAYER_PHIPOSITION[1]=0.;
    LAYER_PHIPOSITION[2]=0.;
    LAYER_PHIPOSITION[3]=0.;
    LAYER_PHIPOSITION[4]=0.;
		
    // position of layers in z relative to vertical center of the payload
    // radius that antennas sit at on the payload
    for (int i=0;i<5;i++) {
      RRX[i]=0.5;
    }
		
  } //else if (Ross payload)
    // Smex payload (full ANITA flown 2006-2007)
  else if (fSettings->WHICH==Payload::Anita1Simple) {
    // fSettings->CYLINDRICALSYMMETRY=1;
		
    //these are physical layers
    // NRX_PHI[0]=8;
    // NRX_PHI[1]=8;
    // NRX_PHI[2]=16;
    // NRX_PHI[3]=8;
		
    PHITRIG[0]=16; // number of positions in phi in each *trigger* layer
    PHITRIG[1]=16;
    PHITRIG[2]=8;
		
    //these are physical layers again
    PHI_OFFSET[0]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
    PHI_OFFSET[1]=-2.*icemc::constants::PI/(double)NRX_PHI[0]/2.;
    PHI_OFFSET[2]=-2.*icemc::constants::PI/(double)NRX_PHI[0]/2.;
    PHI_OFFSET[3]=-2.*icemc::constants::PI/(double)NRX_PHI[0]/4.;
		
    //double INCLINE_NADIR=55; // this is set in the input file now EWG: so why not remove it?
		
    // sets their declination
    THETA_ZENITH[0]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[2]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[3]=icemc::constants::PI/2+INCLINE_NADIR*icemc::constants::RADDEG;
		
    // radius from center axis of the payload
    RRX[0] = 0.9210802;
    RRX[1] = 0.7553198;
    RRX[2] = 2.0645374;
    RRX[3]=1.40; // this is wrong, but I put this version of icemc here just to show one possible strategy to use V polarization and nadirs so it doesn't matter
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0;
    LAYER_VPOSITION[1] = -0.9670034;
    LAYER_VPOSITION[2] = -3.730752;
    LAYER_VPOSITION[3]=(-3.175-0.948-0.889); // this is wrong too, but nadirs are probably gone so doesn't matter.
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1] = 0.;
    LAYER_HPOSITION[2] = 0.;
    LAYER_HPOSITION[3]=0.; // this is wrong too, but nadirs are probably gone so doesn't matter.
		
		
  } //else if (SMEX payload - default)
  else if (fSettings->WHICH==Payload::Custom) {
		
    // cout << "Is this configuration cylindrically symmetric? Yes(1) or No(0)\n";
    // cin >> fSettings->CYLINDRICALSYMMETRY;
		
    // cout << "How many layers?\n";
    // cin >> fSettings->NLAYERS;
		
    for (int i=0;i<fSettings->NLAYERS;i++) {
			
      // cout << "How many antennas in the " << i << "th layer?\n";
      // cin >> NRX_PHI[i];
			
      std::cout << "What is the offset in phi for the " << i << "th layer?\n";
      std::cin >> PHI_OFFSET[i];
			
      std::cout << "What is the theta ascent for the " << i << "th layer (0 if pointed straight upwards, PI if pointed straight downwards)?\n";
      std::cin >> THETA_ZENITH[i];
			
      std::cout << "What is the vertical position of this layer relative to the vertical center of the payload?";
      std::cin >> LAYER_VPOSITION[i];
			
      std::cout << "What is the distance between of the vertical axis of the payload and the center of this layer?";
      std::cin >> LAYER_HPOSITION[i];
			
      std::cout << "What is the phi of this layer relative to the vertical center of the payload in the horizontal plane?";
      std::cin >> LAYER_PHIPOSITION[i];
			
			
			
      if (fSettings->CYLINDRICALSYMMETRY==0) {
	for (int j=0;j<NRX_PHI[i];j++) {
	  std::cout << "What is the phi of the " << j << "th antenna is this layer?\n";
	  std::cin >> PHI_EACHLAYER[i][j];
	} //for (read antenna phi)
      }//if (not cylindrically symmetric)
    } //for (antenna layers)

    
    // cout << "How many polarizations must pass a voltage threshold?\n";
    // cin >> fSettings->NFOLD;

    ///@todo restore maxthreshold!
    // std::cout << "How many times the expected noise level should the voltage threshold be?\n";
    // std::cin >> maxthreshold;
  } //else if (custom payload)

  else if (fSettings->WHICH==Payload::AnitaHill) {// anita hill

    //fSettings->NFOLD=1; // how many channels must pass the trigger
    //maxthreshold=2.3;

    // if (fSettings->NLAYERS!=2)
    //   cout << "Warning!!! Did not enter the right number of layers in the input file.  For Anita Hill, it's 2.";

    // fSettings->CYLINDRICALSYMMETRY=1;

    // NRX_PHI[0]=1; // this is how many antennas we have in phi on each "layer"
    // NRX_PHI[1]=1; // for anita hill, we are calling each station a different "layer"

    ///@todo reenable ANITA-HILL?
    // PHI_OFFSET[0]=(bn1->BN_LONGITUDE+0.)*icemc::constants::RADDEG; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
    // PHI_OFFSET[1]=(bn1->BN_LONGITUDE+0.)*icemc::constants::RADDEG;

    std::cout << "phi_offsets are " << PHI_OFFSET[0] << " " << PHI_OFFSET[1] << "\n";
		
    //double INCLINE_NADIR=55; SET IN INPUT NOW!!!!
		
    // sets their declination
    THETA_ZENITH[0]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
		
    // radius from center axis of the payload
    RRX[0] = 1.;
    RRX[1] = 1.;
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0.;
    LAYER_VPOSITION[1] = 0.;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1] = 100.; // in meters
		
    LAYER_PHIPOSITION[0]=0.;
    LAYER_PHIPOSITION[1] = 30.*icemc::constants::RADDEG;// in radians
  }
    
  else if(fSettings->WHICH==Payload::Anita1) { // Kurt's measurements for the first flight in elog 345
    //fSettings->NFOLD=8;
    // fSettings->CYLINDRICALSYMMETRY=0;
		
    // NRX_PHI[0]=8;
    // NRX_PHI[1]=8;
    // NRX_PHI[2]=16;
		
    PHITRIG[0]=16; // number of positions in phi in each trigger layer
    PHITRIG[1]=16;
    PHITRIG[2]=8;
		
		
    THETA_ZENITH[0]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[2]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
		
    PHI_OFFSET[0]=0.;
    PHI_OFFSET[1]=0.;
    PHI_OFFSET[2]=0.;
		
    ANTENNA_POSITION_START[0][0][0] = MINCH * TVector3(40.957,-39.29,125.38);
    ANTENNA_POSITION_START[0][0][0].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][1] = MINCH * TVector3(57.608,0.606,124.97);
    ANTENNA_POSITION_START[0][0][1].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][2] = MINCH * TVector3(41.049,40.605,124.896);
    ANTENNA_POSITION_START[0][0][2].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][3] = MINCH * TVector3(1.032,57.128,124.93);
    ANTENNA_POSITION_START[0][0][3].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][4] = MINCH * TVector3(-38.832,40.508,125.607);
    ANTENNA_POSITION_START[0][0][4].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][5] = MINCH * TVector3(-55.545,0.549,125.851);
    ANTENNA_POSITION_START[0][0][5].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][6] = MINCH * TVector3(-38.793,-39.423,126.105);
    ANTENNA_POSITION_START[0][0][6].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][7] = MINCH * TVector3(1.113,-55.918,125.731);
    ANTENNA_POSITION_START[0][0][7].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][0] = MINCH * TVector3(19.841,-45.739,87.738);
    ANTENNA_POSITION_START[0][1][0].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][1] = MINCH * TVector3(46.959,-18.601,87.364);
    ANTENNA_POSITION_START[0][1][1].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][2] = MINCH * TVector3(46.983,19.646,87.208);
    ANTENNA_POSITION_START[0][1][2].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][3] = MINCH * TVector3(19.823,46.633,87.194);
    ANTENNA_POSITION_START[0][1][3].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][4] = MINCH * TVector3(-18.429,46.496,87.486);
    ANTENNA_POSITION_START[0][1][4].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][5] = MINCH * TVector3(-45.439,19.34,88.0);
    ANTENNA_POSITION_START[0][1][5].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][6] = MINCH * TVector3(-45.446,-18.769,88.183);
    ANTENNA_POSITION_START[0][1][6].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][7] = MINCH * TVector3(-18.297,-45.857,88.066);
    ANTENNA_POSITION_START[0][1][7].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][0] = MINCH * TVector3(38.622,-94.184,-20.615);
    ANTENNA_POSITION_START[0][2][0].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][1] = MINCH * TVector3(71.662,-72.036,-20.966);
    ANTENNA_POSITION_START[0][2][1].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][2] = MINCH * TVector3(93.857,-39.130,-21.512);
    ANTENNA_POSITION_START[0][2][2].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][3] = MINCH * TVector3(101.664,-0.202,-21.942);
    ANTENNA_POSITION_START[0][2][3].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][4] = MINCH * TVector3(93.993,38.733,-22.351);
    ANTENNA_POSITION_START[0][2][4].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][5] = MINCH * TVector3(71.92,71.815,-22.386);
    ANTENNA_POSITION_START[0][2][5].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][6] = MINCH * TVector3(38.944,93.885,-22.274);
    ANTENNA_POSITION_START[0][2][6].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][7] = MINCH * TVector3(0.036,101.619,-21.813);
    ANTENNA_POSITION_START[0][2][7].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][8] = MINCH * TVector3(-38.991,93.809,-21.281);
    ANTENNA_POSITION_START[0][2][8].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][9] = MINCH * TVector3(-71.899,71.704,-20.777);
    ANTENNA_POSITION_START[0][2][9].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][10] = MINCH * TVector3(-94.121,38.754,-20.825);
    ANTENNA_POSITION_START[0][2][10].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][11] = MINCH * TVector3(-101.986,-0.133,-20.71);
    ANTENNA_POSITION_START[0][2][11].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][12] = MINCH * TVector3(-94.175,-39.024,-20.671);
    ANTENNA_POSITION_START[0][2][12].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][13] = MINCH * TVector3(-72.138,-72.132,-20.295);
    ANTENNA_POSITION_START[0][2][13].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][14] = MINCH * TVector3(-39.111,-94.25,-20.23);
    ANTENNA_POSITION_START[0][2][14].RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][15] = MINCH * TVector3(-0.163,-101.975,-20.229);
    ANTENNA_POSITION_START[0][2][15].RotateZ(-gps_offset);
    PHI_EACHLAYER[0][0] = -45.103 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][1] = -0.14 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][2] = 44.559 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][3] = 89.959 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][4] = 135.555 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][5] = 179.651 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][6] = -135.14 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[0][7] = -90.18 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][0] = -67.283 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][1] = -23.004 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][2] = 22.72 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][3] = 67.82 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][4] = 112.698 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][5] = 157.565 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][6] = -157.376 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[1][7] = -112.449 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][0] = -67.26 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][1] = -45.284 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][2] = -22.457 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][3] = 0.227 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][4] = 22.318 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][5] = 45.008 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][6] = 67.751 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][7] = 89.913 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][8] = 113.016 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][9] = 135.608 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][10] = 157.487 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][11] = 179.709 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][12] = -157.569 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][13] = -135.021 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][14] = -112.773 * icemc::constants::RADDEG - gps_offset;
    PHI_EACHLAYER[2][15] = -89.959 * icemc::constants::RADDEG - gps_offset;
    ANTENNA_DOWN[0][0] = 10.422 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][1] = 10.207 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][2] = 10.714 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][3] = 10.381 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][4] = 10.026 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][5] = 9.515 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][6] = 9.677 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][7] = 9.544 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][0] = 10.183 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][1] = 10.44 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][2] = 10.562 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][3] = 10.655 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][4] = 10.265 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][5] = 9.77 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][6] = 9.422 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][7] = 9.526 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][0] = 9.364 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][1] = 9.712 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][2] = 9.892 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][3] = 10.253 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][4] = 10.574 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][5] = 10.62 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][6] = 10.416 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][7] = 10.189 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][8] = 9.776 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][9] = 9.596 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][10] = 9.561 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][11] = 9.695 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][12] = 9.445 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][13] = 9.387 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][14] = 9.398 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][15] = 9.288 * icemc::constants::RADDEG;
    for(int iii = 0; iii < 3; iii++) // move from the square centers to the phase centers
      for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++)
	ANTENNA_POSITION_START[1][iii][jjj] = ANTENNA_POSITION_START[0][iii][jjj] = ANTENNA_POSITION_START[0][iii][jjj] - phase_center * TVector3(cos(PHI_EACHLAYER[iii][jjj])*sin(90.*icemc::constants::RADDEG+ANTENNA_DOWN[iii][jjj]), sin(PHI_EACHLAYER[iii][jjj])*sin(90.*icemc::constants::RADDEG+ANTENNA_DOWN[iii][jjj]), cos(90.*icemc::constants::RADDEG+ANTENNA_DOWN[iii][jjj]));
  }
  else if (fSettings->WHICH==Payload::EeVEX) {
		
    // Just one layer of antennas around balloon
    // EeVEX
		
    //   fSettings->NFOLD=1; //
    //maxthreshold=2.3;
		
    // fSettings->CYLINDRICALSYMMETRY=1;
		
    // NRX_PHI[0]=360;
    // NRX_PHI[1]=360;
    // NRX_PHI[2]=360;
    // NRX_PHI[3]=360;
    // NRX_PHI[4]=360;
		
    PHITRIG[0]=360;
    PHITRIG[1]=360;
    PHITRIG[2]=360;
    PHITRIG[3]=360;
    PHITRIG[4]=360;
		
		
    PHI_OFFSET[0]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[1]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[2]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[3]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[4]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
		
    // sets their declination
    THETA_ZENITH[0]=icemc::constants::PI/2+5.5*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+7.5*icemc::constants::RADDEG;
    THETA_ZENITH[2]=icemc::constants::PI/2+9.5*icemc::constants::RADDEG;
    THETA_ZENITH[3]=icemc::constants::PI/2+11.5*icemc::constants::RADDEG;
    THETA_ZENITH[4]=icemc::constants::PI/2+13.5*icemc::constants::RADDEG;
		
		
    // radius from center axis of the payload
    RRX[0] = 0.9210802;
    RRX[1] = 0.9210802;
    RRX[2] = 0.9210802;
    RRX[3] = 0.9210802;
    RRX[4] = 0.9210802;
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0;
    LAYER_VPOSITION[1]=10.;
    LAYER_VPOSITION[2]=20.;
    LAYER_VPOSITION[3]=30.;
    LAYER_VPOSITION[4]=40.;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1]=0.;
    LAYER_HPOSITION[2]=0.;
    LAYER_HPOSITION[3]=0.;
    LAYER_HPOSITION[4]=0.;
		
  } //else if (EeVEX)
    //anitaII or satellite
  else if (fSettings->WHICH==Payload::Anita2) {
    std::cout << "initializing and using anitaII payload geometry" << std::endl;
		
    // layer 0 is antennas 1-8 on the payload
    // layer 1 is antennas 9-15
    // layer 2 is antennas 16-32
		
    //fSettings->NFOLD=8;
    //maxthreshold=2.3;
		
    // fSettings->CYLINDRICALSYMMETRY=0;
		
    // NRX_PHI[0]=8;
    // NRX_PHI[1]=8;
    // NRX_PHI[2]=16;
    // NRX_PHI[3]=8;
		
    PHITRIG[0]=16; // number of positions in phi in each trigger layer
    PHITRIG[1]=16;
    PHITRIG[2]=8;
		
    PHI_OFFSET[0]=0.;
    PHI_OFFSET[1]=0;
    PHI_OFFSET[2]=0;
    PHI_OFFSET[3]=0;
		
    // sets their declination
    THETA_ZENITH[0]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[2]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[3]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
		
    ANTENNA_POSITION_START[0][0][0] = MINCH * TVector3(40.438,-36.958,147.227);
    ANTENNA_POSITION_START[0][0][1] = MINCH * TVector3(57.134,3.109,146.476);
    ANTENNA_POSITION_START[0][0][2] = MINCH * TVector3(40.549,43.106,145.871);
    ANTENNA_POSITION_START[0][0][3] = MINCH * TVector3(0.624,59.688,145.361);
    ANTENNA_POSITION_START[0][0][4] = MINCH * TVector3(-39.455,43.147,145.928);
    ANTENNA_POSITION_START[0][0][5] = MINCH * TVector3(-56.096,3.177,146.894);
    ANTENNA_POSITION_START[0][0][6] = MINCH * TVector3(-39.355,-36.753,147.757);
    ANTENNA_POSITION_START[0][0][7] = MINCH * TVector3(0.645,-53.539,147.876);
    ANTENNA_POSITION_START[0][1][0] = MINCH * TVector3(19.554,-43.890,109.531);
    ANTENNA_POSITION_START[0][1][1] = MINCH * TVector3(46.600,-16.625,108.889);
    ANTENNA_POSITION_START[0][1][2] = MINCH * TVector3(46.587,21.659,108.220);
    ANTENNA_POSITION_START[0][1][3] = MINCH * TVector3(19.476,48.539,107.671);
    ANTENNA_POSITION_START[0][1][4] = MINCH * TVector3(-18.798,48.502,107.852);
    ANTENNA_POSITION_START[0][1][5] = MINCH * TVector3(-45.899,21.424,108.516);
    ANTENNA_POSITION_START[0][1][6] = MINCH * TVector3(-45.895,-16.821,109.354);
    ANTENNA_POSITION_START[0][1][7] = MINCH * TVector3(-18.691,-43.864,109.843);
    ANTENNA_POSITION_START[0][2][0] = MINCH * TVector3(38.636,-93.988,2.636);
    ANTENNA_POSITION_START[0][2][1] = MINCH * TVector3(71.690,-72.108,1.953);
    ANTENNA_POSITION_START[0][2][2] = MINCH * TVector3(93.897,-39.211,0.498);
    ANTENNA_POSITION_START[0][2][3] = MINCH * TVector3(101.790,-0.212,-0.661);
    ANTENNA_POSITION_START[0][2][4] = MINCH * TVector3(94.047,38.773,-1.788);
    ANTENNA_POSITION_START[0][2][5] = MINCH * TVector3(72.080,71.816,-2.223);
    ANTENNA_POSITION_START[0][2][6] = MINCH * TVector3(39.065,93.999,-2.561);
    ANTENNA_POSITION_START[0][2][7] = MINCH * TVector3(0.121,101.815,-2.314);
    ANTENNA_POSITION_START[0][2][8] = MINCH * TVector3(-38.815,94.002,-2.034);
    ANTENNA_POSITION_START[0][2][9] = MINCH * TVector3(-71.809,71.912,-1.102);
    ANTENNA_POSITION_START[0][2][10] = MINCH * TVector3(-93.886,39.000,-0.673);
    ANTENNA_POSITION_START[0][2][11] = MINCH * TVector3(-101.885,0.048,0.102);
    ANTENNA_POSITION_START[0][2][12] = MINCH * TVector3(-94.017,-38.841,0.865);
    ANTENNA_POSITION_START[0][2][13] = MINCH * TVector3(-72.079,-71.902,1.864);
    ANTENNA_POSITION_START[0][2][14] = MINCH * TVector3(-39.152,-93.935,2.464);
    ANTENNA_POSITION_START[0][2][15] = MINCH * TVector3(-0.290,-101.771,2.991);
    ANTENNA_POSITION_START[0][3][0] = MINCH * TVector3(32.625,-82.045,-71.140);
    ANTENNA_POSITION_START[0][3][1] = MINCH * TVector3(79.071,-35.639,-72.809);
    ANTENNA_POSITION_START[0][3][2] = MINCH * TVector3(79.172,30.988,-74.893);
    ANTENNA_POSITION_START[0][3][3] = MINCH * TVector3(32.608,77.414,-75.342);
    ANTENNA_POSITION_START[0][3][4] = MINCH * TVector3(-33.398,78.088,-74.957);
    ANTENNA_POSITION_START[0][3][5] = MINCH * TVector3(-79.367,31.568,-73.922);
    ANTENNA_POSITION_START[0][3][6] = MINCH * TVector3(-78.900,-34.192,-72.645);
    ANTENNA_POSITION_START[0][3][7] = MINCH * TVector3(-33.046,-81.696,-70.907);
    PHI_EACHLAYER[0][0] = -45.012 * icemc::constants::RADDEG ;//ant 7
    PHI_EACHLAYER[0][1] = -0.588 * icemc::constants::RADDEG ;//ant 0
    PHI_EACHLAYER[0][2] = 45.694 * icemc::constants::RADDEG ;//ant 1
    PHI_EACHLAYER[0][3] = 90.310 * icemc::constants::RADDEG ;//ant 2
    PHI_EACHLAYER[0][4] = 135.161 * icemc::constants::RADDEG ;//ant3
    PHI_EACHLAYER[0][5] = 179.861 * icemc::constants::RADDEG ;//ant4
    PHI_EACHLAYER[0][6] = -134.930 * icemc::constants::RADDEG ;//ant5
    PHI_EACHLAYER[0][7] = -90.638 * icemc::constants::RADDEG ;//ant 6
    PHI_EACHLAYER[1][0] = -67.412 * icemc::constants::RADDEG ;//ant 15
    PHI_EACHLAYER[1][1] = -23.005 * icemc::constants::RADDEG ;//ant 8
    PHI_EACHLAYER[1][2] = 22.503 * icemc::constants::RADDEG ;//ant 9
    PHI_EACHLAYER[1][3] = 67.722 * icemc::constants::RADDEG ;//ant 10
    PHI_EACHLAYER[1][4] = 112.614 * icemc::constants::RADDEG ;//ant 11
    PHI_EACHLAYER[1][5] = 157.685 * icemc::constants::RADDEG ;//ant 12
    PHI_EACHLAYER[1][6] = -156.639 * icemc::constants::RADDEG ;//ant 13
    PHI_EACHLAYER[1][7] = -112.587 * icemc::constants::RADDEG ;//ant 14
    PHI_EACHLAYER[2][0] = -67.365 * icemc::constants::RADDEG ;//ant 29 
    PHI_EACHLAYER[2][1] = -45.135 * icemc::constants::RADDEG ;//ant 30
    PHI_EACHLAYER[2][2] = -23.002 * icemc::constants::RADDEG ;//ant 31
    PHI_EACHLAYER[2][3] = -1.013 * icemc::constants::RADDEG ;//ant 16
    PHI_EACHLAYER[2][4] = 21.934 * icemc::constants::RADDEG ;//ant 17
    PHI_EACHLAYER[2][5] = 44.467 * icemc::constants::RADDEG ;//ant 18
    PHI_EACHLAYER[2][6] = 67.288 * icemc::constants::RADDEG ;//ant 19
    PHI_EACHLAYER[2][7] = 89.971 * icemc::constants::RADDEG ;//ant 20
    PHI_EACHLAYER[2][8] = 112.390 * icemc::constants::RADDEG ;//ant 21
    PHI_EACHLAYER[2][9] = 134.988 * icemc::constants::RADDEG ;//ant 22
    PHI_EACHLAYER[2][10] = 157.387 * icemc::constants::RADDEG ;//ant 23
    PHI_EACHLAYER[2][11] = 179.843 * icemc::constants::RADDEG ;//ant 24
    PHI_EACHLAYER[2][12] = -157.444 * icemc::constants::RADDEG ;//ant 25
    PHI_EACHLAYER[2][13] = -134.877 * icemc::constants::RADDEG ;//ant 26
    PHI_EACHLAYER[2][14] = -112.406 * icemc::constants::RADDEG ;//ant 27
    PHI_EACHLAYER[2][15] = -90.081 * icemc::constants::RADDEG ;//ant 28
    PHI_EACHLAYER[3][0] = -67.997 * icemc::constants::RADDEG ;//ant 
    PHI_EACHLAYER[3][1] = -22.948 * icemc::constants::RADDEG ;//ant 
    PHI_EACHLAYER[3][2] = 22.382 * icemc::constants::RADDEG ;//ant
    PHI_EACHLAYER[3][3] = 67.583 * icemc::constants::RADDEG ;//ant
    PHI_EACHLAYER[3][4] = 112.844 * icemc::constants::RADDEG ;//ant
    PHI_EACHLAYER[3][5] = 157.761 * icemc::constants::RADDEG ;//ant 
    PHI_EACHLAYER[3][6] = -157.896 * icemc::constants::RADDEG ;//ant
    PHI_EACHLAYER[3][7] = -112.791 * icemc::constants::RADDEG ;//ant
    ANTENNA_DOWN[0][0] = 9.637 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][1] = 10.108 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][2] = 11.245 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][3] = 11.291 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][4] = 10.988 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][5] = 9.491 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][6] = 9.027 * icemc::constants::RADDEG;
    ANTENNA_DOWN[0][7] = 8.743 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][0] = 9.445 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][1] = 10.061 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][2] = 10.772 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][3] = 11.484 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][4] = 11.122 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][5] = 10.376 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][6] = 9.410 * icemc::constants::RADDEG;
    ANTENNA_DOWN[1][7] = 9.039 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][0] = 8.233 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][1] = 8.807 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][2] = 9.120 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][3] = 10.352 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][4] = 10.889 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][5] = 11.315 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][6] = 11.402 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][7] = 11.379 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][8] = 10.842 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][9] = 10.725 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][10] = 10.143 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][11] = 10.067 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][12] = 9.503 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][13] = 9.021 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][14] = 8.453 * icemc::constants::RADDEG;
    ANTENNA_DOWN[2][15] = 8.268 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][0] = 8.007 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][1] = 9.817 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][2] = 10.259 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][3] = 11.648 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][4] = 10.271 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][5] = 10.015 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][6] = 10.889 * icemc::constants::RADDEG;
    ANTENNA_DOWN[3][7] = 7.314 * icemc::constants::RADDEG;

    SIMON_DELTA_R[0][0] = -0.0384839;
    SIMON_DELTA_R[0][1] = 0.00634697;
    SIMON_DELTA_R[0][2] = -0.0861167;
    SIMON_DELTA_R[0][3] = 0.0461873;
    SIMON_DELTA_R[0][4] = 0.0153388;
    SIMON_DELTA_R[0][5] = -0.00927728;
    SIMON_DELTA_R[0][6] = 0.0239867;
    SIMON_DELTA_R[0][7] = 0.0125282;
    SIMON_DELTA_R[1][0] = -0.0111636;
    SIMON_DELTA_R[1][1] = -0.0959452;
    SIMON_DELTA_R[1][2] = -0.0330808;
    SIMON_DELTA_R[1][3] = -0.0475617;
    SIMON_DELTA_R[1][4] = 0.0196292;
    SIMON_DELTA_R[1][5] = -0.0190837;
    SIMON_DELTA_R[1][6] = -0.00922367;
    SIMON_DELTA_R[1][7] = -0.0294811;
    SIMON_DELTA_R[2][0] = 0.0140245;
    SIMON_DELTA_R[2][1] = -0.0621836;
    SIMON_DELTA_R[2][2] = -0.0379325;
    SIMON_DELTA_R[2][3] = -0.0108062;
    SIMON_DELTA_R[2][4] = -0.0601935;
    SIMON_DELTA_R[2][5] = -0.0968276;
    SIMON_DELTA_R[2][6] = -0.0348523;
    SIMON_DELTA_R[2][7] = 0.0121726;
    SIMON_DELTA_R[2][8] = 0.0405193;
    SIMON_DELTA_R[2][9] = 0.0239992;
    SIMON_DELTA_R[2][10] = -0.0405203;
    SIMON_DELTA_R[2][11] = -0.00401756;
    SIMON_DELTA_R[2][12] = -0.0362955;
    SIMON_DELTA_R[2][13] = -0.00587152;
    SIMON_DELTA_R[2][14] = -0.00611182;
    SIMON_DELTA_R[2][15] = -0.00321244;
    SIMON_DELTA_R[3][0] = -0.0437687;
    SIMON_DELTA_R[3][1] = -0.0643475;
    SIMON_DELTA_R[3][2] = -0.0804245;
    SIMON_DELTA_R[3][3] = -0.0112675;
    SIMON_DELTA_R[3][4] = 0.0337428;
    SIMON_DELTA_R[3][5] = -0.0525977;
    SIMON_DELTA_R[3][6] = -0.101587;
    SIMON_DELTA_R[3][7] = -0.0401037;

    SIMON_DELTA_PHI[0][0] = -0.0100608;
    SIMON_DELTA_PHI[0][1] = -0.00313443;
    SIMON_DELTA_PHI[0][2] = -0.015312;
    SIMON_DELTA_PHI[0][3] = 0.00206827;
    SIMON_DELTA_PHI[0][4] = -0.0227948;
    SIMON_DELTA_PHI[0][5] = 0.00750385;
    SIMON_DELTA_PHI[0][6] = 0.00388065;
    SIMON_DELTA_PHI[0][7] = -0.00131021;
    SIMON_DELTA_PHI[1][0] = -0.0299233;
    SIMON_DELTA_PHI[1][1] = -0.00165365;
    SIMON_DELTA_PHI[1][2] = -0.0107407;
    SIMON_DELTA_PHI[1][3] = 0.0145914;
    SIMON_DELTA_PHI[1][4] = -0.0150373;
    SIMON_DELTA_PHI[1][5] = -0.0121967;
    SIMON_DELTA_PHI[1][6] = -0.0038106;
    SIMON_DELTA_PHI[1][7] = 0.0106842;
    SIMON_DELTA_PHI[2][0] = -0.0087849;
    SIMON_DELTA_PHI[2][1] = 0.000682206;
    SIMON_DELTA_PHI[2][2] = -0.00516052;
    SIMON_DELTA_PHI[2][3] = -0.00770935;
    SIMON_DELTA_PHI[2][4] = -0.00862535;
    SIMON_DELTA_PHI[2][5] = -0.00920648;
    SIMON_DELTA_PHI[2][6] = 0.00037431;
    SIMON_DELTA_PHI[2][7] = 0.00310935;
    SIMON_DELTA_PHI[2][8] = -0.00546085;
    SIMON_DELTA_PHI[2][9] = -0.00901249;
    SIMON_DELTA_PHI[2][10] = -0.0145529;
    SIMON_DELTA_PHI[2][11] = -0.00666063;
    SIMON_DELTA_PHI[2][12] = -0.00372999;
    SIMON_DELTA_PHI[2][13] = 0.00197442;
    SIMON_DELTA_PHI[2][14] = -0.000789595;
    SIMON_DELTA_PHI[2][15] = 0.000188257;
    SIMON_DELTA_PHI[3][0] = -0.00289577;
    SIMON_DELTA_PHI[3][1] = -0.0203117;
    SIMON_DELTA_PHI[3][2] = -0.00503387;
    SIMON_DELTA_PHI[3][3] = -0.000220575;
    SIMON_DELTA_PHI[3][4] = -0.00416114;
    SIMON_DELTA_PHI[3][5] = -0.0223176;
    SIMON_DELTA_PHI[3][6] = 0.0058874;
    SIMON_DELTA_PHI[3][7] = 0.00899651;
		
    for(int iii = 0; iii < 4; iii++){ // move from the square centers to the phase centers
      for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++){
			 
	//ANTENNA_DOWN is measured from horiztonal. Put negatives in correct places. Verified with analysis code 
	ANTENNA_POSITION_START[0][iii][jjj] = ANTENNA_POSITION_START[0][iii][jjj] - phase_center_anita2_analysis * TVector3(cos(PHI_EACHLAYER[iii][jjj])*cos(-1*ANTENNA_DOWN[iii][jjj]), sin(PHI_EACHLAYER[iii][jjj])*cos(-1*ANTENNA_DOWN[iii][jjj]), sin(-1*ANTENNA_DOWN[iii][jjj]));
      }//jjj
    }//iii
    double r;
    double phi;
		
    double x;
    double y;
    double z;

	
    for(int iii = 0; iii < 4; iii++){ // move from the square centers to the phase centers
      for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++){
	x = ANTENNA_POSITION_START[0][iii][jjj][0];
	y = ANTENNA_POSITION_START[0][iii][jjj][1];
	z = ANTENNA_POSITION_START[0][iii][jjj][2];

	r = sqrt(pow(x,2)+pow(y,2));
	phi = atan2(y,x);
		   
	ANTENNA_POSITION_START[0][iii][jjj]= TVector3((r+SIMON_DELTA_R[iii][jjj])*cos(phi+SIMON_DELTA_PHI[iii][jjj]),(r+SIMON_DELTA_R[iii][jjj])*sin(phi+SIMON_DELTA_PHI[iii][jjj]),z);

	ANTENNA_POSITION_START[0][iii][jjj].RotateZ(-gps_offset_anita2);
	ANTENNA_POSITION_START[1][iii][jjj]=ANTENNA_POSITION_START[0][iii][jjj];
	PHI_EACHLAYER[iii][jjj]=atan2(ANTENNA_POSITION_START[0][iii][jjj][1],ANTENNA_POSITION_START[0][iii][jjj][0]);//set phi of each antennas to correct starting position

	//cout<<"Antenna pos is "<<ANTENNA_POSITION_START[0][iii][jjj]<<" PHI is "<<PHI_EACHLAYER[iii][jjj]<<"\n";
      }
    }

		
  } 
  else if (fSettings->WHICH==Payload::Anita3 || fSettings->WHICH==Payload::Anita4) { // ANITA-3 and ANITA-4
    std::cout << "initializing and using ANITA-III or IV payload geometry" << std::endl;
    // layer 0 is antennas 1-8 on the payload
    // layer 1 is antennas 9-15
    // layer 2 is antennas 16-32
    // layer 3 is antennas 32-48

    // fSettings->CYLINDRICALSYMMETRY=0;
      
    //these are physical layers
    // NRX_PHI[0]=8;
    // NRX_PHI[1]=8;
    // NRX_PHI[2]=16;
    // NRX_PHI[3]=16;
      
    PHITRIG[0]=16; // number of positions in phi in each *trigger* layer
    PHITRIG[1]=16;
    PHITRIG[2]=16;
      
    //these are physical layers again
    PHI_OFFSET[0]=0.; 
    PHI_OFFSET[1]=0.; // 2.*icemc::constants::PI/(double)NRX_PHI[0]/2.; // Linda: changed this offset to 0 as  it shouldn't be needed
    PHI_OFFSET[2]=0.;
    PHI_OFFSET[3]=0.;
      
    //double INCLINE_NADIR=55; // this is set in the input file now So should be removed
      
    // sets their declination
    THETA_ZENITH[0]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[2]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[3]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
      
    // Read photogrammetry positions
    std::string whichANITAroman="";
    if (fSettings->WHICH==Payload::Anita3) whichANITAroman+="III";
    else whichANITAroman+="IV";
    std::string photoFile;
#ifdef ANITA_UTIL_EXISTS
    const char* anitaUtilInstallEnv = getenv("ANITA_UTIL_INSTALL_DIR");
    if(!anitaUtilInstallEnv){
      icemc::report() << icemc::severity::error << "Unable to find ANITA_UTIL_INSTALL_DIR!" << std::endl;
    }
    photoFile += std::string(anitaUtilInstallEnv) +"/share/anitaCalib/anita"+whichANITAroman+"Photogrammetry.csv";
#else
    photoFile += (ICEMC_DATA_DIR+"/anita"+whichANITAroman+"Photogrammetry.csv");
#endif
    
    std::ifstream Anita3PhotoFile(photoFile.c_str());
    if (!Anita3PhotoFile){
      std::cerr << "Couldn't open photogrammetry!" << std::endl;
      return;
    }

    //First up are the antenna positions
    TString line;
    for(int i=0;i<2;i++) {
      line.ReadLine(Anita3PhotoFile);
      //std::cout << line.Data() << "\n";
    }

    //Array with photogrammetry values
    Double_t xAntPhoto[48]; //inch
    Double_t yAntPhoto[48]; //inch
    Double_t zAntPhoto[48]; //inch
    Double_t rAntPhoto[48]; //inch
    Double_t azCentrePhoto[48]; //deg
    Double_t apertureAzPhoto[48]; //deg
    Double_t apertureElPhoto[48]; //deg

    for(int ant=0;ant<48;ant++) {
      line.ReadLine(Anita3PhotoFile);
      //std::cout << "Seavey:\t" << line.Data() << "\n";
      TObjArray *tokens = line.Tokenize(",");
      for(int j=0;j<8;j++) {
	const TString subString = ((TObjString*)tokens->At(j))->GetString();
	//	TString *subString = (TString*) tokens->At(j);
	//	std::cout << j << "\t" << subString.Data() << "\n";
	switch(j) {
	case 0:
	  if (ant+1 != subString.Atoi()) {
	    std::cerr << "Antenna number mismatch\n";
	  }
	  break;
	case 1:	   
	  xAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 2:	   
	  yAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 3:	   
	  zAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 4:	   
	  rAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 5:	   
	  azCentrePhoto[ant]=subString.Atof(); //deg
	  break;
	case 6:	   
	  apertureAzPhoto[ant]=subString.Atof(); //deg
	  break;
	case 7:	   
	  apertureElPhoto[ant]=subString.Atof()*(-1); //deg // photogrammetry elevation defined as negative, here positive
	  break;
	default:	   
	  break;
	}
	  
      }
      tokens->Delete();
	
    }
    Anita3PhotoFile.close();

    // Fill photogrammetry position for top rings
    // for (int iant=0; iant<8;iant++){
    for(int ant=0; ant < 48; ant++){
      int ilayer = -1, ifold = -1;
      getLayerFoldFromTriggerRX(ant, ilayer, ifold);

      for(int pol = 0; pol < NPOL; pol++){
	ANTENNA_POSITION_START[pol][ilayer][ifold] = MINCH * TVector3(xAntPhoto[ant], yAntPhoto[ant], zAntPhoto[ant]);
	ANTENNA_POSITION_START[pol][ilayer][ifold].RotateZ(-gps_offset_anita3);  // top ring top antennas
      }
      PHI_EACHLAYER[ilayer][ifold] = azCentrePhoto[ant] * icemc::constants::RADDEG - gps_offset_anita3;
      ANTENNA_DOWN[ilayer][ifold] = apertureElPhoto[ant] * icemc::constants::RADDEG;

      // std::cout << ANTENNA_POSITION_START[0][ilayer][ifold].Phi()*TMath::RadToDeg() << std::endl;
    }

    // for (int iant=0; iant<8;iant++){
    //   ANTENNA_POSITION_START[0][0][iant] = MINCH * TVector3(xAntPhoto[iant*2], yAntPhoto[iant*2], zAntPhoto[iant*2]);
      
    //   ANTENNA_POSITION_START[0][1][iant] = MINCH * TVector3(xAntPhoto[iant*2+1], yAntPhoto[iant*2+1], zAntPhoto[iant*2+1]);
    //   ANTENNA_POSITION_START[0][0][iant].RotateZ(-gps_offset_anita3);  // top ring top antennas
    //   ANTENNA_POSITION_START[0][1][iant].RotateZ(-gps_offset_anita3);  // top ring bottom antennas

    //   PHI_EACHLAYER[0][iant] = azCentrePhoto[iant*2] * icemc::constants::RADDEG - gps_offset_anita3;
    //   PHI_EACHLAYER[1][iant] = azCentrePhoto[iant*2+1] * icemc::constants::RADDEG - gps_offset_anita3;
    //   ANTENNA_DOWN[0][iant] = apertureElPhoto[iant*2] * icemc::constants::RADDEG; 
    //   ANTENNA_DOWN[1][iant] = apertureElPhoto[iant*2+1] * icemc::constants::RADDEG; 
 
    // }

    // Fill photogrammetry position for middle and bottom rings
    // for (int iant=0; iant<16;iant++){
    //   ANTENNA_POSITION_START[0][2][iant] = MINCH * TVector3(xAntPhoto[iant+16], yAntPhoto[iant+16], zAntPhoto[iant+16]);
    //   ANTENNA_POSITION_START[0][2][iant].RotateZ(-gps_offset_anita3);	    // middle ring antennas
    //   ANTENNA_POSITION_START[0][3][iant] = MINCH * TVector3(xAntPhoto[iant+32], yAntPhoto[iant+32], zAntPhoto[iant+32]);
    //   ANTENNA_POSITION_START[0][3][iant].RotateZ(-gps_offset_anita3);  // bottom ring antennas

    //   PHI_EACHLAYER[2][iant] = azCentrePhoto[iant+16] * icemc::constants::RADDEG - gps_offset_anita3;
    //   ANTENNA_DOWN[2][iant] = apertureElPhoto[iant+16] * icemc::constants::RADDEG; 

    //   PHI_EACHLAYER[3][iant] = azCentrePhoto[iant+32] * icemc::constants::RADDEG - gps_offset_anita3;
    //   ANTENNA_DOWN[3][iant] = apertureElPhoto[iant+32] * icemc::constants::RADDEG; 

    // }
      
    // HERE HPOL IS 0 AND VPOL IS 1
    std::string whichANITAcard="";
    if (fSettings->WHICH==Payload::Anita3) {
      whichANITAcard+="3";
    }
    else{
      whichANITAcard+="4";
    }
    std::string phaseCenterName;
#ifdef ANITA_UTIL_EXISTS
    phaseCenterName += ( (std::string)getenv("ANITA_UTIL_INSTALL_DIR") +"/share/anitaCalib/phaseCenterPositionsRelativeToPhotogrammetryAnita"+whichANITAcard+".dat");
#else
    phaseCenterName += (ICEMC_DATA_DIR+"/phaseCenterPositionsRelativeToPhotogrammetryAnita"+whichANITAcard+".dat");
#endif
    
    std::ifstream PhaseCenterFile(phaseCenterName.c_str());
    Int_t antNum, tpol, pol;
    Double_t deltaR,deltaPhi,deltaZ;
    char firstLine[180];
    Double_t deltaRPhaseCentre[2][4][16]; //Relative to photogrammetry + ring offset
    Double_t deltaZPhaseCentre[2][4][16]; //Relative to photogrammetry + ring offset
    Double_t deltaPhiPhaseCentre[2][4][16]; //Relative to photogrammetry + ring offset
      
    PhaseCenterFile.getline(firstLine,179);
    // HERE HPOL IS 0 AND VPOL IS 1 that's why we invert pol here
    while(PhaseCenterFile >> antNum >> tpol >> deltaR >> deltaPhi >> deltaZ) {
      int ilayer, ifold;
      getLayerFoldFromTriggerRX(antNum, ilayer, ifold);      
      
      // int ilayer = (antNum<16)*((antNum%2==0)*0 + (antNum%2==1)*1)+ (antNum>15)*(antNum<32)*2+(antNum>31)*3;
      // int ifold = (ilayer<2)*((antNum-ilayer)/2)+(ilayer>1)*(antNum%16);

      if (tpol==1) pol=0;
      else if (tpol==0) pol=1;
	
      deltaRPhaseCentre[pol][ilayer][ifold] = deltaR;
      deltaPhiPhaseCentre[pol][ilayer][ifold] = deltaPhi*TMath::DegToRad();
      deltaZPhaseCentre[pol][ilayer][ifold] = deltaZ;
    } 
    PhaseCenterFile.close();

    std::ifstream relativePhaseCenterToAmpaDelaysFile((ICEMC_DATA_DIR+"/relativePhaseCenterToAmpaDelaysAnita"+whichANITAcard+".dat").c_str());

#ifdef ANITA_UTIL_EXISTS
    AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
    AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance(stoi(whichANITAcard));
    (void) fGeomTool;
    int tempAnt, intTempPol;
    AnitaPol::AnitaPol_t tempPol;
    for(Int_t surf=0; surf<NUM_SURF; surf++){
      for(Int_t chan=0; chan<RFCHAN_PER_SURF; chan++){
	AnitaGeomTool::getAntPolFromSurfChan(surf,chan, tempAnt, tempPol);
    	if (tempAnt!=-1){
	  // in EventReaderRoot 0: HPOL, 1: VPOL, in icemc it's the opposite
	  intTempPol = (tempPol==0) ? 1 : 0;
	  extraCableDelays[intTempPol][tempAnt] = cal->relativePhaseCenterToAmpaDelays[surf][chan]*1e-9;
	}
      }
    }

#else
    for(int ipol=0; ipol<2; ipol++){
      for (int iant=0; iant<48; iant++){
	extraCableDelays[ipol][iant] = 0;
      }
    }
    
#endif

    
    double x, y, z, r, phi;
    for (int ipol = 0; ipol < 2; ipol++){
      for(int ilayer = 0; ilayer < 4; ilayer++){ 
	for(int ifold = 0; ifold < NRX_PHI[ilayer]; ifold++){

	  // First attempt to define phase-centers is done by coming 20 cm
	  // inwards from the antenna face end (that was only defined for VPOL)
	  // that's why the initial ipol=0
	  ANTENNA_POSITION_START[ipol][ilayer][ifold] = ANTENNA_POSITION_START[0][ilayer][ifold] - phase_center_anita3 * TVector3(cos(PHI_EACHLAYER[ilayer][ifold])*cos(ANTENNA_DOWN[ilayer][ifold]), sin(PHI_EACHLAYER[ilayer][ifold])*cos(ANTENNA_DOWN[ilayer][ifold]), sin(ANTENNA_DOWN[ilayer][ifold]));
	  x = ANTENNA_POSITION_START[ipol][ilayer][ifold].X();
	  y = ANTENNA_POSITION_START[ipol][ilayer][ifold].Y();
	  
	  r = sqrt(pow(x,2)+pow(y,2)) + deltaRPhaseCentre[ipol][ilayer][ifold];
	  phi = atan2(y,x) + deltaPhiPhaseCentre[ipol][ilayer][ifold];
	  z = ANTENNA_POSITION_START[ipol][ilayer][ifold].Z() + deltaZPhaseCentre[ipol][ilayer][ifold];	  

	  // if(phi<0) phi+=TMath::TwoPi();
	  // if(phi>TMath::TwoPi()) phi-=TMath::TwoPi();
	  
	  ANTENNA_POSITION_START[ipol][ilayer][ifold]= TVector3(r*cos(phi),r*sin(phi),z);

	  PHI_EACHLAYER[ilayer][ifold]=phi;
	  
	}
      }
    }
  }
    
    
  else if (fSettings->WHICH==Payload::Satellite) { // satellite

    // layer 0 is antennas 1-8 on the payload
    // layer 1 is antennas 9-15
    // fSettings->CYLINDRICALSYMMETRY=1;
		
    // //these are physical layers
    // NRX_PHI[0]=8;
    // NRX_PHI[1]=8;
		
    PHITRIG[0]=8; // number of positions in phi in each *trigger* layer
    PHITRIG[1]=8;
		
    //these are physical layers again
    PHI_OFFSET[0]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
    PHI_OFFSET[1]=-2.*icemc::constants::PI/(double)NRX_PHI[0]/2.;
		
    // sets their declination
    THETA_ZENITH[0]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
    THETA_ZENITH[1]=icemc::constants::PI/2+INCLINE_TOPTHREE*icemc::constants::RADDEG;
		
    // radius from center axis of the payload
    RRX[0] = 0.9210802;
    RRX[1] = 0.7553198;   
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0;
    LAYER_VPOSITION[1] = -7.5;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1] = 0.;
		
		
  } //else if (satellite)

  std::cout << "nantennas is " << fSettings->NANTENNAS << "\n";
    
  // gets noise (vrms) for each bandwidth slice and antenna layer according to antenna theta
    

}//GetPayload
