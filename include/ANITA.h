#ifndef ANITA_SIM_ANITA_FULL_H
#define ANITA_SIM_ANITA_FULL_H

#include "Detector.h"
#include "PropagatingSignal.h"

#include "anita.h"
#include "FlightDataManager.h"
#include "Seavey.h"
#include "AnitaSimOutput.h"

namespace anitaSim {

  /**
   * @class ANITA
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class ANITA : public icemc::Detector, public FlightDataManager, public Anita {
  public:
    ANITA(const Settings* settings);
    virtual ~ANITA();

    virtual double getStartTime() const override {return fFirstRealTime;} // currently in balloon
    virtual double getEndTime() const override {return fLastRealTime;}  // currently in balloon

    virtual int getNumRX() const override {return fNumRX;}
    virtual TVector3 getPositionRX(int i) const override;

    virtual const Geoid::Position& getPosition(double time = TMath::QuietNaN()) override;
    virtual bool applyTrigger() override {return  applyTrigger(-1);}
    virtual bool applyTrigger(int inu);
    virtual void write(const icemc::Event& event) override;
    
    virtual void getDesiredNDt(int& n, double& dt) const override {
      n = 1024;
      dt = 1e-9*1./2.6;
    }

    virtual bool chanceInHell(const icemc::PropagatingSignal& signal) override;

    virtual void addSignalToRX(const icemc::PropagatingSignal& signal, int rx) override {
      addSignalToRX(signal, rx, -1);
    }
    virtual void addSignalToRX(const icemc::PropagatingSignal& signal, int rx, int inu); // just for debugging

    double GetAverageVoltageFromAntennasHit(const Settings *settings1,  int *nchannels_perrx_triggered,  const double *voltagearray,  double& volts_rx_sum) const;


  private:

    const Settings* fSettings;
    int fNumRX;

    std::vector<Seavey> fSeaveys; ///< The set of Seavey antennas on the payload

    ///@todo maybe move or remove these things...? Although they are used in the ROOTified ANITA style data trees...
    VoltsRX fVoltsRX;
    int fL3trig[Anita::NPOL] = {0};  // 16 bit number which says which phi sectors pass L3 V-POL
    // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
    int fL2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX] = {{0}};
    //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
    int fL1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX] = {{0}};

    static const int nAnt = 48; ///@todo I hate this.
    std::vector<double> justNoise_trig[Anita::NPOL][nAnt];
    std::vector<double> justSignal_trig[Anita::NPOL][nAnt];
    std::vector<double> justNoise_dig[Anita::NPOL][nAnt];
    std::vector<double> justSignal_dig[Anita::NPOL][nAnt];
    
    friend class AnitaSimOutput; ///@todo Can I do this and respect privacy with getters?
    AnitaSimOutput fAnitaOutput; ///< Handles converting the MC output into the same format as real ANITA data

    double fLastPositionTime = TMath::QuietNaN();


    void initSeaveys(const Settings *settings1, const Anita *anita1);

    // /** 
    //  * What's the ilayer/ifold of given RX?
    //  * 
    //  * @param rx index of the fSeaveys
    //  * @param ilayer layer of ANITA 
    //  * @param ifold index in phi, maybe...
    //  */
    // void getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const; ///@todo temp public for test

    
  };
}

#endif //ANITASIM_ANITA_FULL_H
