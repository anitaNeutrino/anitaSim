#ifndef ANITA_SIM_GLOBAL_TRIGGER_H
#define ANITA_SIM_GLOBAL_TRIGGER_H

////////////////////////////////////////////////////////////////////////////////////////////////
//class GlobalTrigger:
////////////////////////////////////////////////////////////////////////////////////////////////
#include "TVector3.h"
#include "TriggerState.h"

namespace anitaSim {

  class Settings;
  class Anita;  

  //!  Global Trigger
  class GlobalTrigger {
    
  private:
    const Settings* fSettings = nullptr; 
    
  public:
    GlobalTrigger(const Settings *settings1, const Anita* anita1);
    // void GetArrivalTimes(int inu,Anita* anita1, const TVector3 &rf_direction);
  
    // these are not really used now that we bin in frequency, but we keep them anyway.
  
    double TRIGTIMESTEP;
    int nstepback;  // when a trigger is fired, how many bins you step back to start a coincidence window
    static const int NTRIGPHISECTORS = 16; /// moved from anita.h to where it's actually used...
    // this is used for Anita 3
    //  const double L1_COINCIDENCE[3]={8.E-9,8.E-9,8.E-9}; // L1 coincidence window, in seconds  
    std::array<double, 3> L1_COINCIDENCE_ANITA3; // L1 coincidence window, in seconds  
    // in this scenario B->M is the same as M->B for example
    //  const double L3_COINCIDENCE=22.5E-9; // L3 is now among neighboring phi sectors  
    double LASTTIMETOTESTL1_ANITA3; // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.
    double L3_COINCIDENCE; // L3 is now among neighboring phi sectors.  For Anita 3 and Anita 4 where there is an l3 coincidence  

    // used for an Anita 4 scenario where B->M doesn't have the same window as M->B necessarily
    double L1_COINCIDENCE_MOREGENERAL[3][2]; // L1 coincidence window, in seconds  
    // in this scenario B->M is *not* the same as M->B for example
    double LASTTIMETOTESTL1_ANITA4;
    double LASTTIMETOTESTL1_ANITA4LR_SCA;



    // In LR scenario A, a bottom antenna initiates the trigger, then require coincidences with hits in other layers
    std::array<double, 2> L1_COINCIDENCE_LR_SCA; // L1 coincidence window, in seconds  
    // which layers we're considering LCP, RCP polarizations instead of V,H in scenario A
    double WHICHLAYERSLCPRCP[Anita::NTRIGGERLAYERS_MAX];


    // This if a coincidence between LCP and RCP for a single antenna for Anita 4
    double L1_COINCIDENCE_ANITA4LR_SCB; // L1 coincidence window, in seconds  
    double LASTTIMETOTESTL1_ANITA4LR_SCB;

    std::array<double, 3> L2_COINCIDENCE_ANITA4LR_SCB; // L2 coincidence window, in seconds  
    double LASTTIMETOTESTL2_ANITA4LR_SCB;

    double L3_COINCIDENCE_ANITA4LR_SCB; // L3 coincidence window, in seconds  

    double DELAYS[3]; // delay top, middle, bottom

    std::vector<int> flag_e_L1[Anita::NPHI_MAX];
    std::vector<int> flag_h_L1[Anita::NPHI_MAX];
  
    UShort_t phiTrigMask[Anita::NPOL]; // which phi sector is masked for Anita 2 or Anita-3 (V-POL)
    //UShort_t phiTrigMaskH; // which phi sector is masked for Anita 3 (H-POL)
    UShort_t l1TrigMask[Anita::NPOL]; // which channel is masked for Anita-3 (V-POL)
    //UShort_t l1TrigMaskH; // which channel is masked for Anita 3 (H-POL)

    UShort_t thresholds_eachant[2][48]; // thresholds as read from the surf file: first index is pol, second is antenna number (only working for Anita3)

    double volts[2][Anita::NLAYERS_MAX][Anita::NPHI_MAX];                        // voltage (1st index=antenna,2nd index=pol., lcp=0. rcp=1)
    double volts_em[2][Anita::NLAYERS_MAX][Anita::NPHI_MAX];                        // component of voltage from em shower (1st index=antenna,2nd index=pol., lcp=0. rcp=1)
    double volts_original[2][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //added djg
    // e and h polarizations for each antenna summed across bandwidth
    
    
    int nchannels_perrx_triggered[48]; //Records number of first level triggers on each antenna for a single neutrino
    
    int nchannels_perband_triggered[48][8];  //Records individual channel triggers for each antenna.  (Either a 0 or a 1.) Channels 0-3 are lcp, 4-7 are rcp.
    
    
    
    
    
    // which channels on the payload pass.  
    // 1st component:  layers on the payload
    // 2nd component:  counting around the payload in phi
    // 3rd component:  which polarization
    // 4th component:  slices of bandwidth
    //  int channels_passing[4][16][2][5]; // keeps track of which channels pass
    int channels_passing[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2][Anita::NBANDS_MAX]; // keeps track of which channels pass
    int channels_passing_justNoise[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2][Anita::NBANDS_MAX]; // keeps track of which channels pass from noise triggers
    // make this an array of vectors instead so that we can have an arbitrary number of bands for satellite
    std::vector<int> vchannels_passing[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2];
    std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhits;
    std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhitsall;
    std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhitsnoise; 


    std::array<int, Anita::NTRIG> triggerbits;
    // for the nadir studies
    
    // this is L2 and L3 triggers
    // void PassesTrigger(const Settings *settings1,Anita* anita1,int discones_passing,int mode,
    // 		       int *l3trig,
    // 		       int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],
    // 		       int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],
    // 		       int antennaclump,
    // 		       int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],
    // 		       int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],
    // 		       int inu, int *thispasses);

    void PassesTrigger(Anita* anita1, int mode, TriggerState& triggerState, bool noiseOnly=false);
    
    
    void PassesTrigger(Anita* anita1,int mode,
		       TriggerState& triggerState,
		       double this_threshold,
		       bool noiseOnly=false);

    // void PassesTrigger(const Settings *settings1,Anita* anita1,int discones_passing,int mode,
    // 		       int *l3trig,
    // 		       int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],
    // 		       int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],
    // 		       int antennaclump,
    // 		       int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],
    // 		       int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],
    // 		       int inu,double this_threshold,int *thispasses);
    
  
    // void PassesTriggerBasic(const Settings *settings1,Anita* anita1,int discones_passing,int mode,int *l3trig,
    // 			    int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],
    // 			    int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],
    // 			    int antennaclump,
    // 			    int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],
    // 			    int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],
    // 			    int *thispasses, int inu);

    void PassesTriggerBasic(Anita* anita1,int mode, TriggerState& triggerState, bool noiseOnly=false);
    
    void PassesTriggerCoherentSum(const Anita* anita1, int inu, int *thispasses);

    void PassesTriggerSummedPower(const Anita* anita1);

    void PassesTriggerScheme5(Anita* anita1, double this_threshold, int *thispasses);
  
    void L3Trigger(const Anita* anita1,TriggerState& triggerState);
    // void L3Trigger(const Settings *settings1,Anita* anita1,
    // 		   int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],
    // 		   int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],
    // 		   int discones_passing,int *l3trig,int *thispasses);
    
    
    int GetPhiSector(int i,int j); // given trigger layer and index, get phi sector.
    // for the upper two layers, the phi sector is just j
    // for the nadir layer, the phi sector is 2*j+1
    void GetAnitaLayerPhiSector(int i,int j,int &whichlayer,int &whichphisector);
    void FillInNadir(const Anita* anita1,int *ant);
    void FillInNadir(const Anita* anita1,int ant);
    
    // The following functions are related to the coherent-sum trigger
    std::vector < std::vector < std::vector <double> > > volts_rx_rfcm_trigger; // This is used for easy access to the waveforms of specific phi sectors and layers, and it combines physical layers into their trigger layer.
    // Accessed by volts_rx_rfcm_trigger[phi_sector][layer][timestep]
    int first_phi_sector_hit; // This is used by the coherent waveform sum trigger scheme to make processing more efficient
    double three_bit_round(double input, bool round_zero_up = true, bool allow_zero = false);
    void convert_wfm_to_3_bit(const std::vector <double>& wfm, double rms, std::vector <double>& output);
    void delay_align_antenna_waveforms(const std::vector< std::vector < std::vector <double> > >& waveforms, const std::vector < std::vector <unsigned int> >& delays, std::vector < std::vector <double> >& output);
    void sum_aligned_waveforms(const std::vector < std::vector <double> >& waveforms, std::vector <double>& output);
    double summed_power_window(const std::vector <double>& waveform, unsigned int start_index, unsigned int length);
    // End of functions relating to coherent-sum trigger

    int findahit(const std::vector<int>& myvector,int first,int last);
    int findanl3(int *l3,int NPHISECTORS);

    int L1Anita3_OnePhiSector(int IZERO,std::vector<int> &vl0_realtime_bottom, std::vector<int> &vl0_realtime_middle, std::vector<int> &vl0_realtime_top,
			      std::vector<int> &vl1_realtime_bottom, std::vector<int> &vl1_realtime_middle, std::vector<int> &vl1_realtime_top);

    void L1Anita3_AllPhiSectors(const Anita* anita1,std::array<std::array<std::vector<int>,16>,2> &l1trig);  


    int L1Anita4_OnePhiSector(int IZERO,std::vector<int> &vl0_realtime_bottom, std::vector<int> &vl0_realtime_middle, std::vector<int> &vl0_realtime_top,
			      std::vector<int> &vl1_realtime_bottom, std::vector<int> &vl1_realtime_middle, std::vector<int> &vl1_realtime_top);
    void L1Anita4_AllPhiSectors(const  Anita* anita1,std::array<std::array<std::vector<int>,16>,2> &l1trig);

    void L2Anita3and4(const Anita* anita1,std::array<std::array<std::vector<int>,16>,2> l1trig,
		      std::array<std::array<std::vector<int>,16>,2> &l2trig);

    void L3Anita3and4(const Anita* anita1,std::array<std::array<std::vector<int>,16>,2> vl2trig,
		      int vl3trig[2][16],int *thispasses);


    int PartofL1Anita4LR_ScA_TwoPhiSectors(int ilayerreverse,int IZERO,int ipolar,
					   std::vector<int> &v1l0_realtime_left, std::vector<int> &v2l0_realtime_left, 
					   std::vector<int> &v1l0_realtime_right, std::vector<int> &v2l0_realtime_right, 
					   std::vector<int> &vl1_realtime);
  
    int L1Anita4LR_ScA_TwoPhiSectors(int IZERO,int ipolar,
				     std::vector<int> &vl0_realtime_bottomleft, std::vector<int> &v2l0_realtime_bottomleft, 
				     std::vector<int> &vl0_realtime_bottomright, std::vector<int> &v2l0_realtime_bottomright, 
				     std::vector<int> &vl0_realtime_middleleft, std::vector<int> &v2l0_realtime_middleleft,
				     std::vector<int> &vl0_realtime_middleright, std::vector<int> &v2l0_realtime_middleright,
				     std::vector<int> &vl0_realtime_topleft, std::vector<int> &v2l0_realtime_topleft,
				     std::vector<int> &vl0_realtime_topright, std::vector<int> &v2l0_realtime_topright,
				     std::vector<int> &vl1_realtime_bottom, 
				     std::vector<int> &vl1_realtime_middle, 
				     std::vector<int> &vl1_realtime_top);
  
    void L1Anita4LR_ScA_AllPhiSectors(const Anita* anita1,std::array<std::array<std::vector<int>,16>,2> &vl1trig);
    void L3Anita4LR_ScA(const Anita* anita1,std::array<std::array<std::vector<int>,16>,2> vl2trig,
			int *thispasses);




    void L1Anita4LR_ScB_AllAntennas_OneBin(int IZERO, const Anita* anita1, std::array< std::array< std::vector<int>,16>,3> &vl1trig_anita4lr_scb,int &npassesl1);
    // L1 trigger is at the antenna level again.  Just require coincidence between LCP and RCP
    void L1Anita4LR_ScB_OneBin(int IZERO,const std::vector<int>& vleft,const std::vector<int>& vright,
			       std::vector<int> &vl1trig);

    void L2Anita4LR_ScB_AllPhiSectors_OneBin(int IZERO,
					     const Anita* anita1,
					     std::array< std::array< std::vector<int>,16>,3> vl1trig_anita4lr_scb,
					     std::array<std::array<std::vector<int>,3>,16> &vl2_realtime_anita4_scb,
					     int &npassesl2,
					     int &npassesl2_type0);

  
    // keep track of whether you get a coincidence between 1, 2 or 3 antennas in a phi sector with the right windows.

    void L2Anita4LR_ScB_OnePhiSector_OneBin(int IZERO,std::vector<int> vl1_bottom, 
					    std::vector<int> vl1_middle,
					    std::vector<int> vl1_top,
					    std::array<std::vector<int>,3> &vl2_realtime_anita4_scb,int &npassesl2,int &npassesl2_type0);

    int L3or30Anita4LR_ScB_TwoPhiSectors_OneBin(int IZERO, 
						std::array<std::vector<int>,3> vl2_realtime_anita4_scb, // 3 neighbors, whether 1, 2 or 3 pass
						std::array<std::vector<int>,3> vl2_realtime_anita4_scb_other, // 3 neighbors, whether 1, 2 or 3 pass
						int npass1,int npass2);
  
    void L3Anita4LR_ScB_OneBin(int IZERO, const Anita* anita1,std::array<std::array<std::vector<int>,3>,16> vl2_realtime_anita4_scb,
			       std::array<std::vector<int>,16> &vl3trig_type0, std::array<std::vector<int>,16> &vl3trig_type1,
			       int &thispasses_l3type0,int &thispasses_l3type1);

  
    void delayL0(std::vector<int> &vl0,double delay);
    void delay_AllAntennas(const Anita* anita1);








  private:
    unsigned cwst_event_number;
    unsigned cwst_center_phi_sector;
    double cwst_rms_noise;
    double cwst_actual_rms;
    double cwst_threshold;
    unsigned cwst_window_start;
    unsigned cwst_window_end;
    double cwst_deg_theta;
    double cwst_deg_phi;
    double cwst_actual_deg_theta;
    double cwst_actual_deg_phi;
    TVector3 cwst_rf_direction;
    TVector3 cwst_0th_sector_position;
    double cwst_timesteps[Anita::HALFNFOUR];
    RX cwst_RXs[48];
    RX cwst_aligned_wfms[9];
    std::vector <double> cwst_summed_wfm;
    std::vector <double> cwst_power_of_summed_wfm;
    double cwst_power;
    


  };
}

#endif //ANITA_SIM_GLOBAL_TRIGGER_H
