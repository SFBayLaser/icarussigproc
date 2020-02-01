#ifndef WAVEFORMPARAMSALG_H
#define WAVEFORMPARAMSALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       WaveformParamsAlg
// Module Type: producer
// File:        WaveformParamsAlg.h
//
//              The intent of this module is to provide methods for
//              characterizing an input RawDigit waveform
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//                         computing the means and rms
// RMSRejectionCutHi     - vector of maximum allowed rms values to keep channel
// RMSRejectionCutLow    - vector of lowest allowed rms values to keep channel
// RMSSelectionCut       - vector of rms values below which to not correct
// MaxPedestalDiff       - Baseline difference to pedestal to flag
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
// Based on work done by Brian Kirby, Mike Mooney and Jyoti Joshi
//
////////////////////////////////////////////////////////////////////////
#include "lardataobj/RawData/RawDigit.h"

namespace icarussigproc
{
class WaveformParamsAlg
{
public:

    // Copnstructors, destructor.
    WaveformParamsAlg();
    ~WaveformParamsAlg();

    // Provide definitions of the raw waveforms for internal use
      using VectorFloat = std::vector<float>;
      using ArrayFloat  = std::vector<VectorFloat>;
    
    // Basic waveform mean and rms
    void getMeanAndRms(VectorFloat& rawWaveform,
                       float&       aveVal,
                       float&       rmsVal,
                       int&         numBins) const;
    
    // Basic waveform mean and rms plus trunated rms
    void getMeanAndTruncRms(VectorFloat& rawWaveform,
                            float&       aveVal,
                            float&       rmsVal,
                            float&       rmsTrunc,
                            int&         numBins) const;

    // Truncated rms calculation
    void getTruncatedRMS(VectorFloat& rawWaveform,
                         float&       pedestal,
                         float&       truncRms) const;
    
private:

};

} // end of namespace caldata

#endif
