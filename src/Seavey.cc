#include "Seavey.h"
#include "EnvironmentVariable.h"
#include "Report.h"
#include "Constants.h"
#include "AnitaSimSettings.h"

#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TFile.h" // @todo for debugging

#include <algorithm>

#if defined(ANITA_UTIL_EXISTS) and defined(VECTORIZE)
#include "vectormath_trig.h"
#endif



std::ostream& operator<<(std::ostream& os, anitaSim::Seavey::Pol pol){
  switch(pol){
  case anitaSim::Seavey::Pol::H:
    os << "H";
    break;
  case anitaSim::Seavey::Pol::V:
    os << "V";
    break;
  }
  return os;
}






 // number of points within bandwidth that gain is measured.
constexpr int numGainPoints = 131;
std::array<double, numGainPoints> freqHz_measured; ///< Frequency values for arrays with numGainPoints, goes from 200 MHz to 1.5 GHz with step size of 10 MHz (units in Hz)

/**
 * Ped's famous gains
 */
std::array<double, numGainPoints> gain_dBi_hh_measured; ///< gain for HPol component of signal to HPol feed (co-pol)
std::array<double, numGainPoints> gain_dBi_vv_measured; ///< gain for VPol component of signal to VPol feed (co-pol)
std::array<double, numGainPoints> gain_dBi_hv_measured; ///< gain for HPol component of signal to VPol feed (cross-pol)
std::array<double, numGainPoints> gain_dBi_vh_measured; ///< gain for VPol component of signal to HPol feed (cross-pol)


/**
 * The antenna heights (m), calculated from Ped's gains, converts an E-field (V/m) to a voltage (V)
 */
std::array<double, numGainPoints> heightHH_m; ///< converted from gain_dBi_hh_measured
std::array<double, numGainPoints> heightVV_m; ///< converted from gain_dBi_vv_measured
std::array<double, numGainPoints> heightHV_m; ///< converted from gain_dBi_hv_measured
std::array<double, numGainPoints> heightVH_m; ///< converted from gain_dBi_vh_measured

// If one day icemc becomes even more modular, you come up with some geometry where some are in ice and others not or something
// refractiveIndexUsedForHeightArrays will need to change to be a per-antenna member variable (rather than a global)
// maybe it should be anyway because this is a wee bit ugly
double refractiveIndexUsedForHeightArrays = 0; ///< Required to find height from gains

constexpr int numAnglePoints = 7;
constexpr std::array<double, numAnglePoints> referenceAnglesDeg {0, 5, 10, 20, 30, 45, 90}; // the off axis measurements are have these step sizes (Degrees)
const std::array<double, numAnglePoints> referenceAnglesRad {referenceAnglesDeg[0]*icemc::constants::RADDEG,
							     referenceAnglesDeg[1]*icemc::constants::RADDEG,
							     referenceAnglesDeg[2]*icemc::constants::RADDEG,
							     referenceAnglesDeg[3]*icemc::constants::RADDEG,
							     referenceAnglesDeg[4]*icemc::constants::RADDEG,
							     referenceAnglesDeg[5]*icemc::constants::RADDEG,
							     referenceAnglesDeg[6]*icemc::constants::RADDEG};

std::array<std::array<double, numGainPoints>, numAnglePoints> gain_v_angle_az;
std::array<std::array<double, numGainPoints>, numAnglePoints> gain_h_angle_az;
std::array<std::array<double, numGainPoints>, numAnglePoints> gain_v_angle_el;
std::array<std::array<double, numGainPoints>, numAnglePoints> gain_h_angle_el;

bool doneLoadGains = false; ///< So we only read the Seavey gains and calculate derived quantities once.







/** 
 * Utility function for loadGains
 * 
 * @param fileName the file to open
 * @param failHard exit on failure if true
 * 
 * @return the (hopefully) opened ifstream
 */
std::ifstream openCarefully(const char* fileName, bool failHard = true){
  std::ifstream dataFile;
  dataFile.open(fileName);
  if(dataFile.fail()){
    if(failHard){
      icemc::report() << icemc::severity::error << "Quiting because I can't open " << fileName << "\n";
      exit(1);
    }
    else{
      icemc::report() << icemc::severity::warning << "I can't open "  << fileName << "\n";
    }
  }
  return dataFile;
}





/** 
 * Read the gains into the arrays.
 * @todo This is NOT THREAD SAFE, but icemc won't be threaded for the forseeable future
 */
void loadGains(){
  if(!doneLoadGains){

    const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
    const std::string ICEMC_DATA_DIR=ICEMC_SRC_DIR+"/data/";

    std::vector<std::string> boresightFileNames {"hh_0", "vv_0",  "hv_0", "vh_0"};
    bool firstFile = true; // fill freqHz_measured first time, check against it on all others
    for(const auto& fileName : boresightFileNames){
      std::string fullPath = ICEMC_DATA_DIR + fileName;
      std::ifstream boresight_file = openCarefully(fullPath.c_str());

      auto& boresightGain = (fileName == "vv_0" ? gain_dBi_hh_measured :
			     fileName == "hh_0" ? gain_dBi_vv_measured :
			     fileName == "hv_0" ? gain_dBi_hv_measured :
			                          gain_dBi_vh_measured);
      
      for(int i = 0; i < numGainPoints; ++i) {
	double f;
	boresight_file >> f >> boresightGain.at(i);
	if(firstFile){
	  freqHz_measured.at(i) = f;
	}
	else if (f != freqHz_measured.at(i)){
	  icemc::report() << icemc::severity::warning << "frequency = " << f << ", freqHz_measured[i] = " << freqHz_measured.at(i) << "\n";
	}
      }
      firstFile = false;
    }

    std::vector<std::string> angleFileNames {"vv_az", "hh_az", "vv_el", "hh_el"};

    for(const auto& fileName : angleFileNames){

      auto& gainVsAngle = (fileName == "vv_az" ? gain_v_angle_az :
			   fileName == "hh_az" ? gain_h_angle_az :
			   fileName == "vv_el" ? gain_v_angle_el :
			                         gain_h_angle_el );

      std::string filePath = ICEMC_DATA_DIR + fileName;
      std::ifstream angle_file = openCarefully(filePath.c_str());

      for(auto& g0 : gainVsAngle.at(0)){
	g0 = 1; // 0th bin is on boresight, so no off-axis effects... so multiplier is 1.
      }

      for(int j = 1; j < numAnglePoints; j++){
	int k = 0;
	double f = 0;
	for(auto& g : gainVsAngle.at(j)){	  
	  angle_file >> f >> g;
	  if(f != freqHz_measured.at(k)){
	    icemc::report() << icemc::severity::warning << "Check off-axis frequencies for " << filePath << std::endl;
	  }
	  k++;
	}
      }
      
      angle_file.close();

    }

    doneLoadGains = true;
  }
}


/** 
 * Convert a gain in dBi to antenna height, requires knowing the refractive index of the receiving material
 * 
 * @param gain_dB 
 * @param freq 
 * @param nmedium_receiver 
 * 
 * @return 
 */
inline double getAntennaHeightFromGain_dB(double gain_dB, double freqHz, double nmedium_receiver) {
  // from gain=4*pi*A_eff/lambda^2
  // and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
  // gain_dB is in dB
  using namespace icemc::constants;
  const double c = icemc::constants::CLIGHT;
  const double term1 = (gain_dB*c*c)/(4*icemc::constants::PI*freqHz*freqHz);
  const double Z_air = nmedium_receiver*icemc::constants::Z0;
  const double term2 = icemc::constants::Zr/Z_air;
  return 2*sqrt(term1*term2);
}




/** 
 * We need the refractive index of the antenna medium to calculate the heights.
 * 
 * @param refractiveIndex is the refractive index of the medium
 * @todo This also isn't thread safe
 * 
 */
void fillHeightArrays(double refractiveIndex){
  
  loadGains();
  
  // then we recalculate everything
  if(refractiveIndex != refractiveIndexUsedForHeightArrays){
    
    for(const std::string& polString : {"vv", "hh", "hv",  "vh"}){

      auto& boresightGain = (polString == "vv" ? gain_dBi_vv_measured :
    			     polString == "hh" ? gain_dBi_hh_measured :
    			     polString == "hv" ? gain_dBi_hv_measured :
      			                         gain_dBi_vh_measured);

      auto& antennaHeight = (polString == "hh" ? heightHH_m :
    			     polString == "vv" ? heightVV_m :
    			     polString == "hv" ? heightHV_m :
        			                 heightVH_m);
      int i=0;
      for(auto g : boresightGain){
	antennaHeight.at(i) = getAntennaHeightFromGain_dB(g, freqHz_measured.at(i), refractiveIndex);
	i++;
      }
    }
    refractiveIndexUsedForHeightArrays = refractiveIndex;
  }
}













TCanvas* anitaSim::Seavey::plotInterpolatedGains(double freq, const int nBins){

  auto c = new TCanvas();

  double minAngle = *std::min_element(referenceAnglesDeg.begin(), referenceAnglesDeg.end());
  double maxAngle = *std::max_element(referenceAnglesDeg.begin(), referenceAnglesDeg.end());
 
  TString hName = TString::Format("hInterpGains_%d_%d", nBins, nBins);
  auto hH = new TH2D(hName, "HPol Gain;#deltaAzimuth(Degree);#deltaElevation(degree);Gain factor (no units)", nBins, minAngle, maxAngle,  nBins, minAngle, maxAngle);
  // auto hV = new TH2D("hV", "VPol Gain;#deltaAzimuth(Degree);#deltaElevation(degree);Gain factor (no units)", nBins, minAngle, maxAngle,  nBins, minAngle, maxAngle);

  for(int by=1; by <= hH->GetNbinsY(); by++){
    double dEl = hH->GetYaxis()->GetBinCenter(by);
    for(int bx=1; bx <= hH->GetNbinsX(); bx++){
      double dAz = hH->GetXaxis()->GetBinCenter(bx);
      double w = dAz + dEl;
      hH->Fill(dEl, dAz, w);
    }
  }  
  
  hH->Draw("colz");
  c->SetLogz(1);
  
  return c;
}



TCanvas* anitaSim::Seavey::plotGains() {
  loadGains();

  auto c = new TCanvas();

  c->Divide(1,  2);
  auto c1 = c->cd(1);

  TMultiGraph* grBoresight = new TMultiGraph();

  TGraph* gr_h =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_hh_measured.data());  
  TGraph* gr_v =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_vv_measured.data());
  TGraph* gr_hv = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_hv_measured.data());
  TGraph* gr_vh = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_vh_measured.data());

  gr_h->SetLineColor(kRed);
  gr_v->SetLineColor(kBlue);
  gr_hv->SetLineColor(kGreen);
  gr_vh->SetLineColor(kYellow);
  
  grBoresight->Add(gr_h);
  grBoresight->Add(gr_v);
  grBoresight->Add(gr_hv);
  grBoresight->Add(gr_vh);

  grBoresight->SetTitle("Seavey Boresight Gains Model in icemc;Frequency (Hz);Gain (dBi)");
  
  auto l1 = new TLegend();
  l1->AddEntry(gr_h, "HPol Gain", "l");
  l1->AddEntry(gr_v, "VPol Gain", "l");
  l1->AddEntry(gr_hv, "H->V cross-pol Gain", "l");
  l1->AddEntry(gr_vh, "V->H cross-pol Gain", "l");  

  grBoresight->SetBit(kCanDelete);
  grBoresight->Draw("al");

  l1->SetBit(kCanDelete);
  l1->Draw();

  c1->SetLogy(1);

  auto c2 = c->cd(2);
  c2->Divide(2);
  auto c2_1 = c2->cd(1);
  
  TMultiGraph* grAz = new TMultiGraph();
  TMultiGraph* grEl = new TMultiGraph();
  auto l2_1 = new TLegend();
  l2_1->SetNColumns(2);
  auto l2_2 = new TLegend();
  l2_2->SetNColumns(2);
  for(int j=0; j < numAnglePoints; j++){
    TGraph* grH = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_h_angle_az.at(j).data());
    TGraph* grV = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_v_angle_az.at(j).data());

    TGraph* grH2 = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_h_angle_el.at(j).data());
    TGraph* grV2 = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_v_angle_el.at(j).data());
    
    grH->SetLineColor(kRed);
    grV->SetLineColor(kBlue);
    grH->SetLineStyle(j+1);
    grV->SetLineStyle(j+1);
    grH2->SetLineColor(kRed);
    grV2->SetLineColor(kBlue);
    grH2->SetLineStyle(j+1);
    grV2->SetLineStyle(j+1);
    
    grAz->Add(grV);
    grAz->Add(grH);
    grEl->Add(grV2);
    grEl->Add(grH2);

    TString legTextV = TString::Format("VPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    TString legTextH = TString::Format("HPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    l2_1->AddEntry(grV, legTextV, "l");    
    l2_1->AddEntry(grH, legTextH, "l");

    TString legTextV2 = TString::Format("VPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    TString legTextH2 = TString::Format("HPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    l2_2->AddEntry(grV2, legTextV2, "l");    
    l2_2->AddEntry(grH2, legTextH2, "l");
  }
  
  grAz->SetTitle("Off boresight gains - Azimuth;Frequency (Hz); Relative gain (unitless)");
  grAz->Draw("al");
  grAz->SetBit(kCanDelete);
  l2_1->SetBit(kCanDelete);
  l2_1->Draw();
  c2_1->SetLogy(1);
  auto c2_2 = c2->cd(2);
  grEl->SetTitle("Off boresight gains - Elevation;Frequency (Hz); Relative gain (unitless)");
  grEl->Draw("al");
  grEl->SetBit(kCanDelete);
  l2_2->SetBit(kCanDelete);
  l2_2->Draw();
  c2_2->SetLogy(1);
  
  return c;
}


TCanvas* anitaSim::Seavey::plotHeights(double refractiveIndex) {

  fillHeightArrays(refractiveIndex);

  auto c = new TCanvas();

  TMultiGraph* grBoresight = new TMultiGraph();

  TGraph* gr_h =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightHH_m.data());  
  TGraph* gr_v =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightVV_m.data());
  TGraph* gr_hv = new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightHV_m.data());
  TGraph* gr_vh = new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightVH_m.data());

  gr_h->SetLineColor(kRed);
  gr_v->SetLineColor(kBlue);
  gr_hv->SetLineColor(kGreen);
  gr_vh->SetLineColor(kYellow);
  
  grBoresight->Add(gr_h);
  grBoresight->Add(gr_v);
  grBoresight->Add(gr_hv);
  grBoresight->Add(gr_vh);

  grBoresight->SetTitle(TString::Format("Seavey Antenna Heights in icemc with n=%4.2lf;Frequency (Hz);Height (m)", refractiveIndex));
  
  auto l1 = new TLegend();
  l1->AddEntry(gr_h, "HPol Height", "l");
  l1->AddEntry(gr_v, "VPol Height", "l");
  l1->AddEntry(gr_hv, "H->V cross-pol Height", "l");
  l1->AddEntry(gr_vh, "V->H cross-pol Height", "l");  

  grBoresight->SetBit(kCanDelete);
  grBoresight->Draw("al");

  l1->SetBit(kCanDelete);
  l1->Draw();

  c->SetLogy(1);
  return c;
}


/** 
 * Linearly interpolate in one of the arrays of length numPoints
 * 
 * @param arr some type with .at() and .size() members whose values correspond to the freqHz_measured array (try std::vector or std::array)
 * @param freqHz the frequency at which to interpolate for
 */

template <typename T>
double linearInterp(const T& arr, double freqHz) {

  const double f0 = freqHz_measured.at(0);
  const double df = freqHz_measured.at(1) - f0;
  
  const int i1 = floor((freqHz - f0)/df);
  const int i2 = i1 + 1;

  const double f1 = f0 + i1*df;

  const double a1 = i1 < 0 || i1 >= arr.size() ? 0 : arr.at(i1);
  const double a2 = i2 < 0 || i2 >= arr.size() ? 0 : arr.at(i2);

  const double m = (a2 - a1)/df;
  const double a_interp = a1 + m*(freqHz - f1);

  return a_interp;
}


double anitaSim::Seavey::getHeight(Pol pol, double freqHz) const { 
  fillHeightArrays(fRefractiveIndex);

  switch(pol){
  case Pol::V:
    return linearInterp(heightVV_m, freqHz);
  case Pol::H:
    return linearInterp(heightHH_m, freqHz);
  default:
    icemc::report() << icemc::severity::warning << "Requested height for unknown pol, return 0" << std::endl;
    return 0;
  }
}


double anitaSim::Seavey::getHeight(XPol xPol, double freqHz) const {
  fillHeightArrays(fRefractiveIndex);

  switch(xPol){
  case XPol::VtoH:
    return linearInterp(heightVH_m, freqHz);
  case XPol::HtoV:
    return linearInterp(heightHV_m, freqHz);
  default:
    icemc::report() << icemc::severity::warning << "Requested height for unknown cross-pol, returning 0" << std::endl;
    return 0;
  }
}


/** 
 * Gets the index of the last reference angle less than angleRad
 * 
 * Whether you want Az/El doesn't matter since both are the same.
 * To be used for interpolating off-axis response.
 * If this returns j, then interpolate between j and j+1
 *  
 * @param angleRad is the off-axis angle in radians
 * 
 * @return the index of the last reference angle
 */
int getLowerAngleBin(double angleRad){

  int j1 = -1;
  angleRad = fabs(angleRad);
  for(int j=0; j < referenceAnglesRad.size(); j++){
    if(angleRad >= referenceAnglesRad.at(j)){
      j1 = j;
    }
    else {
      break;
    }
  }
  return j1;
}






bool anitaSim::Seavey::freqAllowedByPassBands(double freqHz) const {

  if(fPassBandsHz.size()==0){
    return true;
  }

  for(auto p : fPassBandsHz){
    // if(freqHz >= p.first && freqHz < p.second){
    if(freqHz >= p.first && freqHz <= p.second){ ///@todo make a choice here
      return true;
    }
  }
  return false;
}





anitaSim::Seavey::Seavey(const TVector3& positionV, const TVector3& positionH,
		      const TVector3& ePlane, const TVector3& hPlane, const TVector3& normal,
		      const Settings* settings, double refractiveIndexOfMedium) : // remove include if you move this
  fPositionV(positionV),
  fPositionH(positionH),
  fEPlane(ePlane),
  fHPlane(hPlane),
  fNormal(normal),
  fRefractiveIndex(refractiveIndexOfMedium)
{
  if(settings){
    fPassBandsHz.emplace_back(std::make_pair(settings->FREQ_LOW_SEAVEYS, settings->FREQ_HIGH_SEAVEYS));
    //std::cout << "The pass band (Hz) is " << fPassBandsHz.back().first << "\t" << fPassBandsHz.back().second << std::endl;
  }
}






double anitaSim::Seavey::getOffAxisResponse(Pol pol, AngleDir dir, double freqHz, double angleRad) const {
  loadGains();

  const int j1 = getLowerAngleBin(angleRad); // calls fabs()
  const int j2 = j1 + 1;

  if(j1 < 0 || j2 >= referenceAnglesRad.size()){
    return 0;
  }
  
  double g1 = 0;
  double g2 = 0;
  switch(pol){
  case Pol::V:
    g1 = dir == AngleDir::Azimuth ? linearInterp(gain_v_angle_az.at(j1), freqHz) : linearInterp(gain_v_angle_az.at(j1), freqHz);
    g2 = dir == AngleDir::Azimuth ? linearInterp(gain_v_angle_az.at(j2), freqHz) : linearInterp(gain_v_angle_az.at(j2), freqHz);
    break;
  case Pol::H:
    g1 = dir == AngleDir::Azimuth ? linearInterp(gain_h_angle_az.at(j1), freqHz) : linearInterp(gain_h_angle_az.at(j1), freqHz);
    g2 = dir == AngleDir::Azimuth ? linearInterp(gain_h_angle_az.at(j2), freqHz) : linearInterp(gain_h_angle_az.at(j2), freqHz);
    break;
  default:
    icemc::report() << icemc::severity::warning << "Requested off-axis gain for for unknown pol, " << (int)pol <<", returning 0\n";
    break;    
  }

  const double a1 = referenceAnglesRad.at(j1); // since looked up from j1, fabs() has been called...
  const double a2 = referenceAnglesRad.at(j2); // since looked up from j2, fabs() has been called...

  const double da = a2 - a1;
  const double m = (g2 - g1)/da;

  //... so don't forget to put a fabs() in here!
  const double g = m*(fabs(angleRad) - a1) + g1;

  // if(fDebug){
  //   std::cout << a2 << "\t" << a1 << "\t" << g2 << "\t" << g1 << "\t" << m << "\t" << fabs(angleRad) << std::endl;
  // }
  
  return g;
}





void anitaSim::Seavey::addSignal(const icemc::PropagatingSignal& s) {

  // if (fDebug == true)
  //   std::cout << "Debugging!" << std::endl;

  
  ///@todo Put the factor of 0.5 for "voltage dividing" elsewhere, like in an actual splitter class downstream of the Seavey...

  double e_component_kvector = 0;
  double h_component_kvector = 0;
  double n_component_kvector = 0;
  anitaSim::Seavey::GetEcompHcompkvector(fEPlane.global,  fHPlane.global,  fNormal.global, s.poynting,
					 e_component_kvector, h_component_kvector, n_component_kvector);

  double e_component = 0;
  double h_component = 0;
  double n_component = 0;
  anitaSim::Seavey::GetEcompHcompEvector(fEPlane.global, fHPlane.global, s.polarization, e_component, h_component, n_component);

  double hitangle_e = 0;
  double hitangle_h = 0;
  anitaSim::Seavey::GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector,
				 hitangle_e, hitangle_h);

  const double df_Hz = s.waveform.getDeltaF();
  const double dt = s.waveform.getDeltaT();

  // Make copies for the VPol and HPol feeds
  icemc::FTPair thisHPol = s.waveform;
  icemc::FTPair thisVPol = s.waveform;

  if(fDebug){
    static int ant = -1;
    ant++;
    const char* opt = ant == 0 ? "recreate" : "update";
    TFile* f = TFile::Open("fSeaveysDebug.root", opt);
    f->cd();

    TGraph grV = thisVPol.getTimeDomain();
    grV.SetName(TString::Format("%d_grV_before", ant));
    grV.Write();

    TGraph grH = thisHPol.getTimeDomain();
    grH.SetName(TString::Format("%d_grH_before", ant));
    grH.Write();
    
    f->Write();
    f->Close();
  }

  std::vector<std::complex<double> >& vPolFreqs = thisVPol.changeFreqDomain();
  std::vector<std::complex<double> >& hPolFreqs = thisHPol.changeFreqDomain();  

  /**
   * In anita.cc, they loop over 0,1 for get_gain_angle(hitangle_e)
   * then they loop over 2, 3 for get_gain_angle(hitangle_h)		 
   * the old index over looping variables goes into the gain_angle arrays
   * 0 is vv_az, 1 is hh_az, 2 is hh_el, 3 is vv_el.
   */

  /**
   * @todo So this is what I think anita.cc was doing
   * I'm not sure it's right, but if it's wrong, it's the wrong off-axis
   * response of the cross-pol, so probably a small effect
   */

  bool temp = fDebug;
  fDebug = false;
  double freqHz = 0;
  for(auto& c : vPolFreqs){

    if(freqAllowedByPassBands(freqHz)){

      const double heightVV = getHeight(Pol::V, freqHz); // VPol component of signal to VPol feed
      const double heightHV = getHeight(XPol::HtoV, freqHz); // HPol component of signal to VPol feed via cross-pol antenna response

      const double offAxisResponseV  = getOffAxisResponse(Pol::V, AngleDir::Azimuth, freqHz, hitangle_e);
      // const double offAxisResponseHV = getOffAxisResponse(Pol::H, AngleDir::Azimuth, freqHz, hitangle_e);
      const double offAxisResponseHV = getOffAxisResponse(Pol::H, AngleDir::Elevation, freqHz, hitangle_h);

      
      if(fDebug){
        std::cout << "Seavey V: " << freqHz/1e6 << "MHz" 
		  << " heightVV=" << heightVV
		  << " oarV=" << offAxisResponseV
		  << " e_comp=" << e_component
		  << " hitangle_e=" << hitangle_e
		  << " totalGainFactor=" << 0.5*sqrt(heightVV*heightVV*e_component*e_component*offAxisResponseV)
		  << " dt=" << dt << "\n";
	
        // std::cout << "Seavey V: " << freqHz/1e9 << ", heights=("
	// 	<< heightVV << ", " << heightHV << "), oars=("
	// 	<< offAxisResponseV << ", " << offAxisResponseHV << "), e_comp = "
	// 	<< e_component << ", h_comp" << h_component << ", hitangle_e = "
	// 	<< hitangle_e << "\n";
      }
 
      // 0.5 is for voltage dividing apparently, it doesn't happen in the Seavey... but it does happen downstream... maybe
      //const double totalGainFactorV = 0.5*sqrt(  heightVV*heightVV*e_component*e_component*offAxisResponseV
      //					       + heightHV*heightHV*h_component*h_component*offAxisResponseHV );
      double totalGainFactorV = 0.5*sqrt(heightVV*heightVV*e_component*e_component*offAxisResponseV);
      c *= totalGainFactorV/(sqrt(2)*dt*1.E6); // factor from icemc::ChanTrigger right before AntennaGain is applied
    }
    else{
      c = 0;
    }
    
    freqHz += df_Hz;
  }

  //fDebug = temp;
  //fDebug=false;
  freqHz = 0; // freqHz is incremented in the loop, so reset
  for(auto& c : hPolFreqs){

    if(freqAllowedByPassBands(freqHz)){

      // get everything going into the HPol feed... via direct and cross-pol.
      const double heightHH = getHeight(Pol::H, freqHz);
      const double heightVH = getHeight(XPol::VtoH, freqHz);

      // then you need to take acconut of how far off boresight you are... i.e. the off-axis reponse of the antennas.
      const double offAxisResponseH  = getOffAxisResponse(Pol::H, AngleDir::Elevation, freqHz, hitangle_h);
      const double offAxisResponseVH = getOffAxisResponse(Pol::V, AngleDir::Elevation, freqHz, hitangle_h);

      // 0.5 is for voltage dividing apparently, it doesn't happen in the Seavey... but it does happen downstream... maybe
      // double totalGainFactorH = 0.5*sqrt(  heightHH*heightHH*h_component*h_component*offAxisResponseH
      // 					 + heightVH*heightVH*e_component*e_component*offAxisResponseVH);

      double totalGainFactorH = 0.5*sqrt(heightHH*heightHH*h_component*h_component*offAxisResponseH);
					   
      
      if(fDebug){
           std::cout << "Seavey H: " << freqHz/1e6 << "MHz" 
		     << " heightHH=" << heightHH
		     << " oarH=" << offAxisResponseH
		     << " e_comp=" << e_component
		     << " hitangle_e=" << hitangle_e
		     << " totalGainFactor=" << 0.5*sqrt(heightHH*heightHH*h_component*h_component*offAxisResponseH)
		     << " dt=" << dt << "\n";

      }
      c *= totalGainFactorH/(sqrt(2)*dt*1.E6);
    }
    else{
      c = 0;
    }
    
    freqHz += df_Hz;
  }

  /**
   * @todo In order to make SCREEN stuff work, make this additive rather than just the most recent
   * This will require doing some addition of the waveform. And FTPair doesn't do a simple +=.
   */
  fHPol = thisHPol;
  fVPol = thisVPol;

  // if(startEnergy > 0){
  //   TGraph gr = fVPol.makePowerSpectrumGraph();
  //   double integral = 0;
  //   for(int i=0; i < gr.GetN(); i++){
  //     integral += gr.GetY()[i];
  //   }
  //   if(integral > 0){
  //     std::cout << startEnergy/integral << "\t" << startEnergy << "\t" << integral << std::endl;
  //   }
  // }

  // {
  //   TGraph gr = fVPol.makePowerSpectrumGraph();
  //   double integral = 0;
  //   for(int i=0; i < gr.GetN(); i++){
  //     integral += gr.GetY()[i];
  //   }
  //   std::cout << integral << std::endl;
  // }  

  //fDebug = temp;
  if(fDebug){
    static int ant = -1;
    ant++;
    const char* opt = "update";
    TFile* f = TFile::Open("fSeaveysDebug.root", opt);
    f->cd();

    TGraph grV = fVPol.getTimeDomain();
    grV.SetName(TString::Format("%d_grV_after", ant));
    grV.Write();

    TGraph grH = fHPol.getTimeDomain();
    grH.SetName(TString::Format("%d_grH_after", ant));
    grH.Write();
    
    f->Write();
    f->Close();
  }

}





const icemc::FTPair& anitaSim::Seavey::get(Pol pol) const {
  
  if(pol==Pol::V){
    return fVPol;
  }
  else if(pol==Pol::H){
    return fHPol;
  }
  else{
    icemc::report() << icemc::severity::warning << "Pol for unknown pol requested, returning V.\n";
    return fVPol;
  }  
}









void anitaSim::Seavey::GetEcompHcompkvector(const TVector3& n_eplane, const TVector3& n_hplane, const TVector3& n_normal,
					 const TVector3 n_exit2bn,
					 double& e_component_kvector, double& h_component_kvector, double& n_component_kvector) {

  // find component along e-plane for the purpose of finding hit angles, that is, in direction of k vector, direction of radio wave)
  e_component_kvector = -(n_exit2bn.Dot(n_eplane));
  // find component along h-plane for the purpose of finding hit angles, that is, in direction of k vector, direction of radio wave)
  h_component_kvector = -(n_exit2bn.Dot(n_hplane));
  // find the component normal
  n_component_kvector = -(n_exit2bn.Dot(n_normal));

} // end GetEcompHcompkvector



void anitaSim::Seavey::GetEcompHcompEvector(const TVector3& n_eplane, const TVector3& n_hplane, const TVector3& n_pol,
					 double& e_component, double& h_component, double& n_component) {

  // find component along e-plane in direction of polarization, that is in direction of the E field   
  e_component = n_pol.Dot(n_eplane);
  //    std::cout << "n_component : " << n_exit2bn << " " << n_normal << " " << n_component << std::endl;
    
  // find component along h-plane in direction of polarization, that is in direction of the E field 
  h_component = n_pol.Dot(n_hplane);


  ///@todo maybe restore this at some point?
  // if (settings1->REMOVEPOLARIZATION) {
  //   //Trying to remove effects of polarization at antenna. Stephen
  //   e_component = n_pol.Dot(n_pol);
  //   h_component = 0.001;
  //   n_component = 0.001;
  // } //if
  
} // end GetEcompHcompEvector


void anitaSim::Seavey::GetHitAngles(double e_component_kvector, double h_component_kvector, double n_component_kvector, double& hitangle_e, double& hitangle_h) {
#if defined(ANITA_UTIL_EXISTS) and defined(VECTORIZE)
  Vec2d y(e_component_kvector, h_component_kvector); 
  Vec2d x(n_component_kvector, n_component_kvector); 
  Vec2d answer = atan2(y,x); 
  hitangle_h = answer[0]; 
  hitangle_e = answer[1]; 

#else
  hitangle_e=atan2(h_component_kvector,n_component_kvector);
  hitangle_h=atan2(e_component_kvector,n_component_kvector);
#endif

}


void anitaSim::Seavey::updatePosition(const Geoid::Position& position, double heading, double pitch, double roll){

  std::vector<VectorPair*> forRotation {&fPositionV, &fPositionH, &fEPlane,  &fHPlane,  &fNormal};

  for(auto& vp : forRotation){

    const double BalloonTheta = position.Theta();
    double BalloonPhi = position.Phi();
  
    if(BalloonPhi > icemc::constants::PI){
      BalloonPhi = BalloonPhi - icemc::constants::TWOPI;
    }

    vp->global = vp->payload;

    const TVector3 zaxis(0.,0.,-1.);
    TVector3 xaxis(1.,0.,0.);  // roll axis
    TVector3 yaxis(0.,-1.,0.); // pitch axis for positive rotation to the clockwise of roll

    // do heading...
    vp->global.Rotate(heading*icemc::constants::RADDEG,zaxis);
    xaxis.Rotate(heading*icemc::constants::RADDEG,zaxis);
    yaxis.Rotate(heading*icemc::constants::RADDEG,zaxis);

    // do pitch...
    vp->global.Rotate(pitch*icemc::constants::RADDEG,yaxis);
    xaxis.Rotate(pitch*icemc::constants::RADDEG,yaxis);

    // do roll...
    vp->global.Rotate(roll*icemc::constants::RADDEG,xaxis);//roll and pitch

    ////now place balloon at proper lat and lon
    // BalloonPhi =latitude*icemc::constants::RADDEG;
    vp->global.RotateY(BalloonTheta);
    vp->global.RotateZ(BalloonPhi);    
  }
  
  // finally translate to payload position
  std::vector<VectorPair*> forTranslation {&fPositionV, &fPositionH};  
  for(auto& vp : forTranslation){
    vp->global += position;
  }  
}



