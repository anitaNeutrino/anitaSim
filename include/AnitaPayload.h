#ifndef ANITA_SIM_ANITA_FULL_H
#define ANITA_SIM_ANITA_FULL_H

#include "Detector.h"
#include "PropagatingSignal.h"

#include "anita.h"
#include "FlightDataManager.h"
#include "Seavey.h"
#include "AnitaSimOutput.h"
#include "TriggerState.h"

namespace anitaSim {

  /**
   * @class AnitaPayload
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class AnitaPayload : public icemc::Detector, public FlightDataManager, public Anita {
  public:
    AnitaPayload(const Settings* settings);
    virtual ~AnitaPayload();

    virtual double                 getStartTime()         const override {return fFirstRealTime;}
    virtual double                 getEndTime()           const override {return fLastRealTime;}
    virtual int                    getNumRX()             const override {return fNumRX;}
    virtual TVector3               getPositionRX(int i)   const override;

    virtual const Geoid::Position& getPosition(double time = TMath::QuietNaN()) override;

    virtual bool applyTrigger()         override {return  applyTrigger(-1);}
    virtual bool applyTrigger(int inu);
    virtual void write(const icemc::Event& event) override;
    
    virtual void getDesiredNDt(int& n, double& dt) const override {
      n = 1024;
      dt = 1e-9*1./2.6;
    }

    virtual bool chanceInHell(const icemc::PropagatingSignal& signal) override;
    virtual double signalThreshold() const override;
    
    virtual void addSignalToRX(const icemc::PropagatingSignal& signal, int rx) override {
      addSignalToRX(signal, rx, -1);
    }
    virtual void addSignalToRX(const icemc::PropagatingSignal& signal, int rx, int inu); // just for debugging

  private:

    const Settings* fSettings;
    int fNumRX;

    std::vector<Seavey> fSeaveys; ///< The set of Seavey antennas on the payload

    VoltsRX fVoltsRX;
    TriggerState fTriggerState;
    std::vector<std::array<std::array<double, 5>, Anita::NPOL> > fThresholdsAnt;
    // double fThresholdsAnt[nAnt][Anita::NPOL][5] = {{{0}}};
    
    friend class AnitaSimOutput; ///@todo Can I do this and respect privacy with getters?
    AnitaSimOutput fAnitaOutput; ///< Handles converting the MC output into the same format as real ANITA data

    double fLastPositionTime = TMath::QuietNaN();

    void initSeaveys();

  };
}

#endif //ANITASIM_ANITA_FULL_H
