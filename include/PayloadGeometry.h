#ifndef ANITA_SIM_PAYLOAD_GEOMETRY_H
#define ANITA_SIM_PAYLOAD_GEOMETRY_H

#include "TVector3.h"


namespace anitaSim {

  class Settings;

  /**
   * @class PayloadGeometry
   * @brief Represent the payload geometry
   */
  class PayloadGeometry {
  private:

  public:
    PayloadGeometry(const Settings* settings);

    /** 
     * What's the ilayer/ifold of given trigger RX?
     * 
     * @param rx index of the fSeaveys
     * @param ilayer layer of ANITA 
     * @param ifold index in phi, maybe...
     */
    void getLayerFoldFromTriggerRX(int rx, int& ilayer, int& ifold) const; ///@todo temp public for test

    int GetRxTriggerNumbering(int ilayer, int ifold) const;// get antenna number based on which layer and position it is
    void GetPayload();
    
    static const int NPOL=2;                                    ///< number of polarizations
    static const int NANTENNAS_MAX=2000;
    static const int NLAYERS_MAX=5;                             ///< max number of layers (in smex design, it's 4)
    static const int NPHI_MAX=400;                              ///< max number of antennas around in phi (in smex, 16)

    int number_all_antennas;                                    ///< this keeps count of the number of antennas for use with timing calculations, etc.

    TVector3 ANTENNA_POSITION_START[NPOL][NLAYERS_MAX][NPHI_MAX]; ///< antenna positions from Kurt's measurements
    double ANTENNA_DOWN[NLAYERS_MAX][NPHI_MAX];                 ///< down angles of antennas from Kurt's measurements
    double SIMON_DELTA_R[NLAYERS_MAX][NPHI_MAX];                ///< measurements by Simon used in analysis ANITA-2
    double SIMON_DELTA_PHI[NLAYERS_MAX][NPHI_MAX];              ///< measurements by Simon used in analysis ANITA-2

    // TVector3 antenna_positions[NPOL][NLAYERS_MAX * NPHI_MAX];     ///< these are the antenna positions in space in a coordinate system where x=north and y=west and the origin is at the center of the payload

    int NRX_PHI[NLAYERS_MAX] = {0};                                   ///< number of antennas around in each layer. (radians)
    double PHI_EACHLAYER[NLAYERS_MAX][NPHI_MAX] = {{0}};                ///< phi of the center of each antenna on each layer
  
    //before correcting for offset for the layer.
    //only used if it is cylindrically symmetric (radians)
    double PHI_OFFSET[NLAYERS_MAX] = {0};                             ///< antenna offset in phi for each layer (radians)
    double THETA_ZENITH[NLAYERS_MAX] = {0};                           ///< how the antenna is tilted in theta (in radians with 0=up)
    // 0=horizontal,+90=down

    double LAYER_VPOSITION[NLAYERS_MAX] = {0};                 ///< position of layers in z relative to vertical center of the payload

    // anita proposal "says that the separation between upper and lower
    // 2 layers of antennas is just under 4m.
    // for anita hill, consider the positions of the "layers" of the "payload" (the stations) to be sitting on the horizontal grid defined by polar coordinates
    
    double LAYER_HPOSITION[NLAYERS_MAX];                 ///< distance in horizontal plane between center axis of the "payload" and each "layer".
    double LAYER_PHIPOSITION[NLAYERS_MAX];               ///< phi corresponding to the position of each "layer" on the "payload"
    double RRX[NLAYERS_MAX] = {0};                             ///< radius that the antenna sits from the axis of the payload (feedpoint)
    Double_t deltaTPhaseCentre[NPOL][NLAYERS_MAX][NPHI_MAX];    ///< Relative to photogrammetry + ring offset

    double INCLINE_TOPTHREE; // cant angle of top three layers of antennas
    double INCLINE_NADIR; // cant angle of nadir (bottom) layer of antennas
    int PHITRIG[NLAYERS_MAX]; // number of positions in phi for each trigger layer
    double extraCableDelays[NPOL][48];

  private:
    const Settings* fSettings;

    
  };
  

}

#endif
