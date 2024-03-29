////////////////////////////////////////////////////////////////////////////////////////////////
//class Anita:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ANITA_SIM_ANITA_H
#define ANITA_SIM_ANITA_H

#include "rx.h"
#include <array>
#include "TVector3.h"
#include "PayloadGeometry.h"
#include "RicianRNG.h"

#ifdef ANITA_UTIL_EXISTS
#include "FFTtools.h"
#endif

#include "TF1.h"

class TGraph;
class TFile;
class TTree;
class TH1F;
class TRandom3;

namespace icemc{
  class Antarctica;
}

namespace anitaSim {
  
  class RX;
  class FlightDataManager;
  class Settings;
  
  //! Contains everything about positions within payload and signals it sees for each event, in both the trigger and signal paths.
  class Anita : public PayloadGeometry {

  private:
    TGraph *gshort[4];
    void setTrigRequirement(int WHICH);
  public:
    ///@todo find a way to remove the balloon from this constructor!
    Anita(const Settings* settings, const char* outputdir, const FlightDataManager* bn1); // constructor
    virtual ~Anita();

    void Initialize(const Settings *settings1, std::ofstream &foutput, TString outputdir); ///< initialize a bunch of stuff
    void initializeFixedPowerThresholds(std::ofstream &foutput);
    void readVariableThresholds(const Settings *settings1);
    void readAmplification();
    void getDiodeDataAndAttenuation(const Settings *settings1, TString outputdir);
    void getPulserData();

    // takes arrays that span NFREQ and turn them into arrays that span HALFNFOUR
    void MakeArrayforFFT(double *vsignalarray_e,double *vsignal_e_forfft, double phasedelay, bool useconstantdelay, bool debug=false) const;
    void GetArrayFromFFT(double *tmp_fftvhz, double *vhz_rx) const;


    int getLabAttn(int NPOINTS_LAB, double *freqlab, double *labattn);
    void labAttn(double *vhz) const;
    void SetNoise(const Settings *settings1, FlightDataManager *bn1, const icemc::Antarctica *antarctica);
    int tuffIndex; // keith edits

    static const int NBANDS_MAX=100;                            ///< max number of bands
    static const int NFREQ=128;
    static const int NTRIG=5;
    static const int NTRIGGERLAYERS_MAX=3;

    double THERMALNOISE_FACTOR;                                 ///< factor to multiply thermal noise for error analysis

    TFile *fnoise;
    TTree *tdiode;

    static constexpr int NFOUR=1024; // Number of fourier point
    static constexpr int HALFNFOUR=NFOUR/2; // Half of the number of fourier points

    // these are used for the satellite thing
    int NBANDS;                                                                                        ///< number of frequency sub-bands (not counting full band)

    // these variables are for filling the tsignals tree
    // double signal_vpol_inanita[5][HALFNFOUR];                                                          ///< this is the signal waveform in the vertical polarization, before converting to LCP, RCP where applicable
    //double noise_vpol_inanita[5][HALFNFOUR];                                                         ///< this is the noise waveform in the vertical polarization, before converting to LCP, RCP where applicable
    // double total_vpol_inanita[5][HALFNFOUR];                                                           ///< this is the sum of the signal and noise in the vertical polarization, before converting to LCP, RCP where applicable

    // double timedomainsignal_rfcm[HALFNFOUR];
    // double timedomainsignal_lab[HALFNFOUR];

    TTree *turfratechain;
    TTree *surfchain;
    TFile *fturf;
    TFile *fsurf;

    UShort_t phiTrigMask;
    UShort_t phiTrigMaskH;
    UShort_t l1TrigMask;
    UShort_t l1TrigMaskH;
    Double_t deadTime;                                                                                 ///< fractional deadTime
    unsigned int realTime_turfrate;                                                                    ///< realtime from the turf rate file
    unsigned int realTime_tr_min;                                                                      ///< min realtime from the turf rate file
    unsigned int realTime_tr_max;                                                                      ///< max realtime from the turf rate file
    unsigned int realTime_surf;                                                                        ///< realtime from the surf file
    unsigned int realTime_surf_min;                                                                    ///< min realtime from the surf file
    unsigned int realTime_surf_max;                                                                    ///< max realtime from the surf file
    UShort_t thresholds[NPOL][48];                                                                        ///< thresholds as read from the surf file: first index is pol, second is antenna number (only working for Anita3)
    UShort_t scalers[NPOL][48];                                                                           ///< scalers as read from the surf file: first index is pol, second is antenna number (only working for Anita3)
    Double_t fakeThresholds[NPOL][48];                                                                    ///< Fake thresholds (coming from converting fake scalers to thresholds)
    Double_t fakeThresholds2[NPOL][48];                                                                   ///< Fake thresholds 2 (coming from converting flight scalers to thresholds)
    Double_t fakeScalers[NPOL][48];                                                                       ///< Fake scalers (coming from converting threhsolds during flight to scalers using threshold scan)

    int iturf;// for indexing
    int isurf;
    int iturfevent;

    static const int npointThresh = 1640;
    Float_t threshScanThresh[NPOL][48][npointThresh];                                                     ///< adc thresholds from threshold scan
    Float_t threshScanScaler[NPOL][48][npointThresh];                                                     ///< scalers from threshold scan
    Float_t minadcthresh[NPOL][48];
    Float_t maxadcthresh[NPOL][48];

    void setphiTrigMaskAnita3(UInt_t realTime_flightdata);
    void setphiTrigMask(UInt_t realTime_flightdata);
    void setTimeDependentThresholds(UInt_t realTime_flightdata);

    double total_diodeinput_1_inanita[5][HALFNFOUR];                                                   ///< this is the waveform that is input to the tunnel diode in the first (LCP or vertical) polarization
    double total_diodeinput_2_inanita[5][HALFNFOUR];                                                   ///< this is the waveform that is input to the tunnel diode in the second (RCP or horizontal) polarization

    double total_diodeinput_1_allantennas[48][HALFNFOUR];                                              ///< this is across all antennas, just the full band
    double total_diodeinput_2_allantennas[48][HALFNFOUR];                                              ///< needs comment

    // these are just for the antenna that receives the signal first
    double timedomain_output_inanita[NPOL][5][HALFNFOUR];                                                 ///< this is just for writing out to the following tree

    double time_trig[HALFNFOUR];
    double weight_inanita; // weight of the event
    int arrayofhits_inanita[3][16][2][HALFNFOUR];
    // same as arrayofhits_inanita but it's time reversed
    int arrayofhits_forgaryanderic[3][16][2][HALFNFOUR];
    int l1trig_anita3and4_inanita[NPOL][16][HALFNFOUR];
    int l1trig_anita4lr_inanita[3][16][HALFNFOUR];
    int l1trig_anita4lr_forgaryanderic[3][16][HALFNFOUR];
    int l2trig_anita4lr_inanita[16][3][HALFNFOUR];
    int l2trig_anita4lr_forgaryanderic[16][HALFNFOUR];                                                           ///< when it passes 2/3
    int l3type0trig_anita4lr_inanita[16][HALFNFOUR];
    int l3trig_anita4lr_inanita[16][HALFNFOUR];
    int l3type0trig_anita4lr_forgaryanderic[16][HALFNFOUR];
    int l3type1trig_anita4lr_forgaryanderic[16][HALFNFOUR];
    double timedomain_output_corrected_forplotting[2][6][HALFNFOUR];                                            ///< this is just for writing out to the following tree
    double timedomain_output_allantennas[2][48][HALFNFOUR];                                                     ///< this is across all antennas, just the full band
    int flag_e_inanita[5][HALFNFOUR];
    int flag_h_inanita[5][HALFNFOUR];
    double dangle_inanita,emfrac_inanita,hadfrac_inanita;
    double ston[5];                                                                                             ///< signal to noise;
    std::array<int, 5> iminbin;                                                                                             ///< this is the minimum bin to start
    std::array<int, 5> imaxbin;
    int maxbin_fortotal[5];                                                                                     ///< when it sums the noise and signal together it shortens the waveform
    int channels_passing[2][5];                                                                                 ///< channels passing.  This is reset for every antenna for every event
    int channels_passing_justNoise[2][5];
    int l1_passing; // l1 passing
    int l1_passing_allantennas[48]; // l1 passing

    double vmmhz_banding[NFREQ];                                                                                ///< V/m/MHz after banding
    double vmmhz_banding_rfcm[NFREQ];                                                                           ///< V/m/MHz after banding and rfcms

    // Note: The following 4 RMS noise variables are for all antennas of all events.
    // In fact, they don't represent RMS until after all events are finished!
    double rms_rfcm[2];                                                                                         ///< rms noise just after rfcm's
    double rms_lab[2];                                                                                          ///< rms noise at lab chip


    TH1F *hsignals[5];                                                                                          ///< s/n (max diode output/mean diode output) for vertical polarization in each band

    std::array<double, NFOUR/4> f_pulser;
    std::array<double, NFOUR/4> f_phases;
    std::array<double, NFOUR/4> f_noise;
    std::array<double, NFOUR/4> v_pulser;
    std::array<double, NFOUR/4> v_phases;
    std::array<double, NFOUR/4> v_noise;



    double cumulat_prob[9];
    double cumulat_prob_plus1[9];


    // for filling tsignals tree
    double timedomainnoise_rfcm_banding[2][5][HALFNFOUR];
    double timedomainnoise_rfcm[2][HALFNFOUR];
    double timedomainnoise_lab[2][HALFNFOUR];

    double phases[5][HALFNFOUR];

    // for filling tglob
    // for each polarization
    int passglobtrig[NPOL];
    double integral_vmmhz_foranita;


    int nnoiseevents;                                                                                          ///< total number of noise events we're choosing from
    int noiseeventcounter;                                                                                     ///< counts which event we're on so we go in order

    double FREQ_LOW;                                                                                           ///< lowest frequency
    double FREQ_HIGH;                                                                                          ///< highest frequency

    double NOTCH_MIN;                                                                                          ///< low edge of notch filter.  This is set in the input file
    double NOTCH_MAX;                                                                                          // upper edge of notch filter

    int BANDING;// set in the input file
    // whether or not you set your own banding (1)
    // or use anita-1 banding

    double freq[NFREQ];  // frequency for each bin
    double freq_forfft[NFOUR]; // frequencies for taking fft of signal
    double freq_forplotting[NFOUR/4]; // just one entry for frequency, unlike the above.

    double freqdomain_rfcm_banding[5][HALFNFOUR/2]; // average noise in frequency domain
    double freqdomain_rfcm[HALFNFOUR/2]; // average noise in frequency domain
    double avgfreqdomain_lab[HALFNFOUR/2]; // average noise in frequency domain
    
    double phases_rfcm_banding[2][5][HALFNFOUR/2];
    double phases_rfcm[2][HALFNFOUR/2];
    double phases_lab[2][HALFNFOUR/2];



    void getDiodeModel();
    void setDiodeRMS(const Settings *settings1, TString outputdir);
    int GetRxTriggerNumbering(int ilayer, int ifold) const;                                                           ///< get antenna number based on which layer and position it is
  
    TF1 fdiode;
    double maxt_diode;
    int idelaybeforepeak[5];
    int iwindow[5];
    double diode_real[5][NFOUR]; // This is the time domain of the diode response. (actually NFOUR/2 array is used.)
    double fdiode_real[5][NFOUR]; // This is the fft of the diode response. (use all NFOUR array. This is for doubling array size for zero padding)

    void myconvlv(double *timedomain_forconvl,const int NFOUR,double *fdiode,double &maxdiodeconvl,double &onediodeconvl,double *power_noise,double *diodeconv);

    static int AntennaNumbertoSurfNumber(int ilayer,int ifold); // find surf where this antenna is triggered
    static int GetAntennaNumber(int ilayer,int ifold);
    static int WhichBand(int ibw,int ipol); // which band, 1-8, in order as they are on the surf

    void Banding(int j, const double *freq_noise,double *powerperfreq,int NPOINTS_NOISE) const;
    void RFCMs(int ilayer,int ifold,double *vmmhz) const;
    void normalize_for_nsamples(double *spectrum, double nsamples, double nsamp);
    void convert_power_spectrum_to_voltage_spectrum_for_fft(int length,double *spectrum, const double domain[], const double phase[]);
    void GetNoiseWaveforms(); // make time domain noise waveform based on avgnoise being the v^2
    void GetPhases();
    
    // each of the above graphs has 601 bins in it
    static const int NPOINTS_BANDS=601;

    double freq_bands[5][NPOINTS_BANDS]; // a frequency array for each of the four bands
    double attn_bands[5][NPOINTS_BANDS]; // attn array for each of the four bands in dB
    double bandsattn[5][NPOINTS_BANDS]; // as a fraction
    //double correl[4][NPOINTS_BANDS]; // correlation between each band and the fullband
    double correl_banding[5][NPOINTS_BANDS]; // correlation between each band and the fullband
    double correl_lab[NPOINTS_BANDS]; // correlation between each band and the fullband
    //double correl[5][NPOINTS_BANDS]; // correlation between each band and the fullband


    static const int NPOINTS_AMPL=58;// bins in amplification
    double freq_ampl[NANTENNAS_MAX][NPOINTS_AMPL]; // frequencies for each amp bin
    double ampl[NANTENNAS_MAX][NPOINTS_AMPL]; // amplification
    double ampl_notdb[NANTENNAS_MAX][NPOINTS_AMPL];// amplification again, but as a fraction
    double noisetemp[NANTENNAS_MAX][NPOINTS_AMPL]; // noise temp each amp bin


    static const int NPOINTS_NOISE=2000;

    //double bwslice_thresholds[4]={2.319,2.308,2.300,2.290}; // this allows you to set different thresholds for each band
    double bwslice_vnoise[NLAYERS_MAX][5]; // expected noise voltage for antenna layer and
    //each slice in bandwidth

    double probability[5];
    double bwslice_enoise[5]; // average integrated power in each band
    double bwslice_fwhmnoise[5]; // 1/2 of fwhm of hnoise
    double bwslice_rmsdiode[5]; // average rms diode output across noise waveforms in each band
    double bwslice_meandiode[5]; // mean diode output across all samples in a sample of noise waveforms generated for each band
    double bwslice_vrms[5]; // rms noise voltage for this bandwidth slice
    double bwslice_dioderms_fullband_allchan[2][48][6]; // diode rms for noise read from flight
    double bwslice_diodemean_fullband_allchan[2][48][6]; // diode rms for noise read from flight
    double freq_noise[5][NPOINTS_NOISE]; // frequency array that goes with vnoise array

    double powerthreshold[5];
    double powerthreshold_nadir[5];
    int NCH_PASS; // for ANITA 3 trigger - requires some number of channels pass

    double l1window; // time window where we require coincidences at L1

    double INTEGRATIONTIME; // integration time of the tunnel diode
    static const int nsamp=100; // number of samples that were used to measure the noise data
    double TIMESTEP; // time step between samples for digitization


    double DEADTIME;

    double TRIG_TIMESTEP; // this is the l1 trigger window for the anita 3 trigger.
    unsigned N_STEPS_PHI;
    unsigned N_STEPS_THETA;

    static const unsigned N_SUMMED_PHI_SECTORS = 4;
    static const unsigned N_SUMMED_LAYERS = 3;

    int GAINS;
    
    double energythreshold; // relative to expected energy from noise
    double MIN_PHI_HYPOTHESIS;
    double MAX_PHI_HYPOTHESIS;
    double MIN_THETA_HYPOTHESIS;
    double MAX_THETA_HYPOTHESIS;

    int USEPHASES;

    int NTRIGGERLAYERS; // number of layers considered by the trigger.  may be different from nlayers, the number of physical layers on the payload.
    // In Anita 1 and Anita 2, the number of physical layers were 3 while the number of trigger layers were 2.
    int REQUIRE_CENTRE; // require centre antenna in clump to be one of those hit

    double diffraction[2][89][NFREQ];
    void SetDiffraction();
    double GetDiffraction(int ilayer, double zenith_angle, int ifreq);



    static const int NPOINTS_LAB=272; // from note 137
    double freqlab[NPOINTS_LAB]; // frequency for each lab attn. bin
    double labattn[NPOINTS_LAB]; // lab attenuation
    double VNOISE[NLAYERS_MAX]; // noise calculated for each antenna layer depending on cant angle- this is only used right now for the chance in hell cuts
    int trigRequirements[NLAYERS_MAX];//  0th element - L1 - how many channels per antenna should pass
    std::array<double,5> bwslice_thresholds;  // thresholds for each band -- this is just an initialization- this is set in the input file
    std::array<double,5> bwslice_allowed; // these bands are allowed to contribute to the trigger sum -- this is set in the input file
    int bwslice_required[5]; // these bands are required to be among the channels that pass -- this is set in the input file
    int pol_allowed[NPOL];// which polarisations are allowed to have channels that fire (V,H)
    int pol_required[NPOL];// which polarisations are required to have channels that fire (V,H)



    double bwslice_center[5]; // center frequencies
    double bwslice_width[5]; // 3 dB bandwidths, without overlap
    double bwslice_min[5]; //minimum of each bandwidth slice
    double bwslice_max[5]; //minimum of each bandwidth slice
    double bwmin; // minimum width of any allowed bandwidth slice

    static const unsigned int NUM_COHERENT_ANTENNAS = 9;
    unsigned hypothesis_offsets[16][200][200][4][3]; // Time bin offsets for each hypothesis - [center_phi_sector_index][phi_angle_index][theta_angle_index][phi_sector_index][layer_index]
    std::vector< std::vector< std::vector <double> > > hypothesis_angles; // Time bin offsets for each hypothesis - [center_phi_sector_index][phi_angle_index][theta_angle_index][phi_sector_index][layer_index]
    //unsigned antenna_indices[16][9];	// These are the indices of the antennas used for a given hypothesis' center phi center index - [center_phi_sector-index][which of the nine]
    std::vector< std::vector <int> > vdifferent_offsets;
    std::vector< std::vector <double> > vdifferent_angles;

    std::array<double, NPHI_MAX> VNOISE_ANITALITE; // noise for each antenna, for the anita-lite trigger configuration.
    double LIVETIME;

    double SIGMA_THETA; // resolution on the polar angle of the signal

#ifdef ANITA_UTIL_EXISTS
    RFSignal *fSignalChainResponseDigitizerTuffs[NPOL][3][16][6]; // 0:VPOL, 1:HPOL ---- 0:TOP, 1:MIDDLE, 2:BOTTOM------- 0:configA, 1:configB, 2:configC, 3:configG, 4:configO, 5:configP
    RFSignal *fSignalChainResponseTriggerTuffs[NPOL][3][16][6];  // same as for DigitizerTuffs
    void readImpulseResponseDigitizer(const Settings *settings1);
    void readImpulseResponseTrigger(const Settings *settings1);
    void readTuffResponseDigitizer(const Settings *settings1);
    void readTuffResponseTrigger(const Settings *settings1);
    void readTriggerEfficiencyScanPulser(const Settings *settings1);
    void readNoiseFromFlight(const Settings *settings1);
    void getQuickTrigNoiseFromFlight(double justNoise[HALFNFOUR], int ipol, int iant, int ituff);
    TGraph *RayleighFits[NPOL][48];
    Int_t numFreqs;
    Double_t *freqs;
    TGraph *gPulseAtAmpa;
    RFSignal *fSignalChainResponseDigitizer[NPOL][3][16]; // 0:VPOL, 1:HPOL ---- 0:TOP, 1:MIDDLE, 2:BOTTOM
    RFSignal *fSignalChainResponseTrigger[NPOL][3][16]; // 0:VPOL, 1:HPOL ---- 0:TOP, 1:MIDDLE, 2:BOTTOM
#endif
    void calculateDelaysForEfficiencyScan();

    Double_t fTimes[HALFNFOUR];
    Double_t fSignalChainResponseA3DigitizerFreqDomain[NPOL][3][16][400];
    Double_t fSignalChainResponseDigitizerFreqDomain[NPOL][3][16][6][400];
    Double_t fSignalChainResponseTriggerFreqDomain[NPOL][3][16][6][400];
    Double_t fRatioTriggerToA3DigitizerFreqDomain[NPOL][3][16][6][400];
    Double_t fRatioDigitizerToA3DigitizerFreqDomain[NPOL][3][16][6][400];
    Double_t deltaT;

    // Trigger efficiency scan parameters
    Int_t trigEffScanPhi;                      // central phi sector of trigger efficiency scan
    std::array<Double_t, 5> trigEffScanAtt;                // attenuations to apply to central and adjecent antennas
    std::array<Double_t, 5> trigEffScanPhiDelay;           // delays between phi sectors
    std::array<Double_t, 3> trigEffScanRingDelay;          // delays between rings
    std::array<Int_t   , 5> trigEffScanApplyRingDelay;     // to which phi sectors apply ring delays 
    std::array<Int_t   , 3> trigEffScanRingsUsed;          // to which rings apply scan
    std::array<Double_t, HALFNFOUR> trigEffScanPulseAtAmpa;
    std::array<Double_t, NFOUR> trigEffScanPulseAtAmpaUpsampled;
    std::array<Double_t, NFREQ> trigEffScanAmplitudeAtAmpa;

    std::array<std::array<Double_t, HALFNFOUR>, 250> trigEffScanPulseAtSurf;
    std::array<int, 3> TUFFstatus;
    int ntuffs;

    icemc::RicianRNG fRician;

  };
}

#endif //ANITA_SIM_ANITA_H
