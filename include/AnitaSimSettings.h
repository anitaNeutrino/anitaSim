#ifndef ANITA_SIM_SETTINGS_H
#define ANITA_SIM_SETTINGS_H

#include "Settings.h"
#include "FlightDataManager.h"

namespace anitaSim {
  class Anita;

  enum class Payload : int
    {
     AnitaLite    = 0,
     Ross         = 1,
     Anita1Simple = 2,
     Custom       = 3,
     AnitaHill    = 4,
     SLAC         = 5, //?
     Anita1       = 6,
     EeVEX        = 7,
     Anita2       = 8,
     Anita3       = 9,
     Anita4       = 10,
     Satellite    = 11
  };
  

  class Settings : public icemc::Settings {
  public:
    Settings();
    virtual ~Settings(){};
    void Initialize();

    virtual void ReadInputs(const char* fileName , std::ofstream &foutput) override;
    void ApplyInputs(Anita*) const;
    
    FlightPath WHICHPATH;    
    int NLAYERS;
    int NANTENNAS;
    int NRX_PHI[Anita::NLAYERS_MAX];
    Payload WHICH; // which payload to use 0=Anita-lite,1=Ross,2=Smex,3=make your own
    int ANITAVERSION;
    int CYLINDRICALSYMMETRY; // is it cylindrically symmetric =1 if which=1,2, =0 if which=0
    // if which=3 then 0 or 1
    int TRIGGERSCHEME;  // frequency domain voltage, frequency domain energy, time domain diode integration
    double INCLINE_TOPTHREE;
    double INCLINE_NADIR;
    int GAINS;
    int BANDING;
    int NBANDS;
    int trigRequirements[4];//  0th element - L1 - how many channels per antenna should pass
    // 1st element- L2 - how many antennas on a layer
    // 2nd element - L3 - how many L2 triggers should be coincident
    int REQUIRE_CENTRE; // require centre antenna in clump to be one of those hit
    double INCLUDE_NADIRONLY; // cant angle of nadir (bottom) layer of antennas
    double SIGMA_THETA; // resolution on the polar angle of the signal
    double FREQ_LOW;       ///< lowest frequency
    double FREQ_HIGH;
    
    int trigEffScanPhi;                      // central phi sector of trigger efficiency scan

    int BN_LATITUDE;
    int BN_LONGITUDE;
    int BN_ALTITUDE;
    int RANDOMIZE_BN_ORIENTATION;
    // int CENTER;                                                                ///< whether or not to center one phi sector of the payload on the incoming signal (for making signal efficiency curves)
    double MAXHORIZON;

    int LCPRCP; // 1 for circular polarization trigger, 0 for V and H
    int JUSTVPOL; // 0 for both polarizations, 1 for just V polarization
    // doesn't allow for both LCPRCP=1 and JUSTVPOL=1
    //int FIFTHBAND; // 1 to include 0.2-1.2 GHz as a frequency band if JUSTVPOL==1
    //int NFOLD=3;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite
    int NFOLD;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite


    //int CHMASKING=1; // whether or not to include channel masking
    //int PHIMASKING=1; // whether or not to include phi masking
    int CHMASKING; // whether or not to include channel masking
    int PHIMASKING; // whether or not to include phi masking

    //int DISCONES=1; // whether or not to use discones
    int DISCONES; // whether or not to use discones


    //double NDISCONES_PASS=3; // number of discones needed to pass
    double NDISCONES_PASS; // number of discones needed to pass

    int BORESIGHTS; // whether to loop over boresights
    int SLAC; // whether or not we are simulating the slac run
    double SLACSLOPE; // slope of the ice
    double SLACICELENGTH;  // length of the block of ice
    double SLAC_HORIZDIST; // horizontal distance from interaction to center of payload at slac beam test
    double SLAC_DEPTH; // vertical depth of interaction at slac beam test
    double SLAC_HORIZ_DEPTH; // horizontal depth of interaction at slac



    double FREQ_LOW_SEAVEYS; // min frequency for seaveys
    double FREQ_HIGH_SEAVEYS; // max frequency for seaveys
    double BW_SEAVEYS;
    //int FORSECKEL=1; // Make array of strength of signal across frequencies for different viewing angles.

    int antennaclump; //number of antenna in clump (L2)
    // End of the once-global varibles.
    double COHERENT_THRESHOLD;
    int APPLYIMPULSERESPONSEDIGITIZER;       // apply impulse response in the digitizer path
    int APPLYIMPULSERESPONSETRIGGER;         // apply impulse response in the trigger path
    int USETIMEDEPENDENTTHRESHOLDS;          // use time-dependent thresholds
    int USEDEADTIME;                         // use dead time from flight
    int NOISEFROMFLIGHTTRIGGER;              // use thermal noise from flight in trigger path
    int NOISEFROMFLIGHTDIGITIZER;            // use thermal noise from flight in digitizer path
    int TRIGGEREFFSCAN;                      // do a trigger efficiency scan
    int TRIGGEREFFSCAPULSE;                  // Apply pulse at AMPA (0) or at SURF (1)

    int TUFFSON;                             // Are the TUFFs on for the whole flight?

    int ADDCW;                               // Add CW
  
  private:
    /** 
     * @brief Some of the logic from Anita moved to here to make Settings const correct
     */
    void setNrxPhiAndNantennasFromWhich();


    std::vector<double> efficiencyScanOffAxisAttenuations;
    std::vector<double> efficiencyScanPhiSectorDelay;
    std::vector<double> efficiencyScanRingDelay;
    std::vector<int> efficiencyScanRingsUsed;
    std::vector<int> efficiencyScanApplyRingDelay;
    std::vector<int> whichTUFFsON;
    std::vector<double> tempThresholds;  
    std::vector<double> bandLowEdgesMHz;
    std::vector<double> bandHighEdgesMHz;
    std::vector<int> requiredBands;
    std::vector<int> allowedBands;
    std::vector<double> notchFilterLimitsMHz;
    std::vector<int> channelRequirePol;
    std::vector<int> channelAllowedPol;
    
    
  };
};




/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param which is the Payload class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const anitaSim::Payload& which);


#endif// ANITA_SIM_SETTINGS_H
