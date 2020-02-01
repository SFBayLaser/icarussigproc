
#include "WaveformParamsAlg.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace icarussigproc
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
WaveformParamsAlg::WaveformParamsAlg()
{
    return;
}
    
//----------------------------------------------------------------------------
/// Destructor.
WaveformParamsAlg::~WaveformParamsAlg()
{}

void WaveformParamsAlg::getTruncatedRMS(VectorFloat& rawWaveform,
                                        float&       pedestal,
                                        float&       truncRms) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(rawWaveform.size());
    
    // Fill the vector
    std::transform(rawWaveform.begin(),rawWaveform.end(),adcLessPedVec.begin(),std::bind(std::minus<short>(),std::placeholders::_1,pedestal));
    
    // sort in ascending order so we can truncate the sume
    std::sort(adcLessPedVec.begin(), adcLessPedVec.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    int minNumBins = 0.6 * rawWaveform.size();
    
    // Get the truncated sum
    truncRms = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.begin() + minNumBins, adcLessPedVec.begin(), 0.);
    truncRms = std::sqrt(std::max(0.,truncRms / double(minNumBins)));
    
    return;
}

void WaveformParamsAlg::getMeanAndRms(VectorFloat& rawWaveform,
                                      float&       aveVal,
                                      float&       rmsVal,
                                      int&         numBins) const
{
    // The strategy for finding the average for a given wire will be to
    // find the most populated bin and the average using the neighboring bins
    // To do this we'll use a map with key the bin number and data the count in that bin
    std::pair<VectorFloat::const_iterator,VectorFloat::const_iterator> minMaxValItr = std::minmax_element(rawWaveform.begin(),rawWaveform.end());
    
    int minVal = std::floor(*minMaxValItr.first);
    int maxVal = std::ceil(*minMaxValItr.second);
    int range  = maxVal - minVal + 1;
    
    std::vector<int> frequencyVec(range, 0);
    int              mpCount(0);
    int              mpVal(0);
    
    for(const auto& val : rawWaveform)
    {
        int intVal = std::round(val) - minVal;
        
        frequencyVec[intVal]++;
        
        if (frequencyVec.at(intVal) > mpCount)
        {
            mpCount = frequencyVec[intVal];
            mpVal   = intVal;
        }
    }
    
    // take a weighted average of two neighbor bins
    int meanCnt  = 0;
    int meanSum  = 0;
    int binRange = std::min(16, int(range/2 + 1));

    for(int idx = mpVal-binRange; idx <= mpVal+binRange; idx++)
    {
        if (idx < 0 || idx >= range) continue;
        
        meanSum += (idx + minVal) * frequencyVec[idx];
        meanCnt += frequencyVec[idx];
    }
    
    aveVal = float(meanSum) / float(meanCnt);
    
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> adcLessPedVec(rawWaveform.size());
    
    std::transform(rawWaveform.begin(),rawWaveform.end(),adcLessPedVec.begin(),std::bind(std::minus<float>(),std::placeholders::_1,aveVal));
    
    // recalculate the rms for truncation
    rmsVal  = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.);
    rmsVal  = std::sqrt(std::max(float(0.),rmsVal / float(adcLessPedVec.size())));
    numBins = meanCnt;
    
    return;
}

void WaveformParamsAlg::getMeanAndTruncRms(VectorFloat& rawWaveform,
                                           float&       aveVal,
                                           float&       rmsVal,
                                           float&       rmsTrunc,
                                           int&         numBins) const
{
    // The strategy for finding the average for a given wire will be to
    // find the most populated bin and the average using the neighboring bins
    // To do this we'll use a map with key the bin number and data the count in that bin
    std::pair<VectorFloat::const_iterator,VectorFloat::const_iterator> minMaxValItr = std::minmax_element(rawWaveform.begin(),rawWaveform.end());
    
    int minVal = std::floor(*minMaxValItr.first);
    int maxVal = std::ceil(*minMaxValItr.second);
    int range  = maxVal - minVal + 1;
    
    std::vector<int> frequencyVec(range, 0);
    int              mpCount(0);
    int              mpVal(0);
    
    for(const auto& val : rawWaveform)
    {
        int intVal = std::round(val) - minVal;
        
        frequencyVec[intVal]++;
        
        if (frequencyVec.at(intVal) > mpCount)
        {
            mpCount = frequencyVec[intVal];
            mpVal   = intVal;
        }
    }
    
    // take a weighted average of two neighbor bins
    int meanCnt  = 0;
    int meanSum  = 0;
    int binRange = std::min(16, int(range/2 + 1));
    
    for(int idx = mpVal-binRange; idx <= mpVal+binRange; idx++)
    {
        if (idx < 0 || idx >= range) continue;
        
        meanSum += (idx + minVal) * frequencyVec[idx];
        meanCnt += frequencyVec[idx];
    }
    
    aveVal = float(meanSum) / float(meanCnt);
    
    // Subtract the pedestal 
    std::transform(rawWaveform.begin(),rawWaveform.end(),rawWaveform.begin(),std::bind(std::minus<float>(),std::placeholders::_1,aveVal));
    
    // do rms calculation - the old fashioned way and over all adc values
    VectorFloat adcLessPedVec = rawWaveform;
    
    // recalculate the rms for truncation
    rmsVal  = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.);
    rmsVal  = std::sqrt(std::max(float(0.),rmsVal / float(adcLessPedVec.size())));
    
    // Drop the "large" rms values and recompute
    std::vector<float>::iterator newEndItr = std::remove_if(adcLessPedVec.begin(),adcLessPedVec.end(),[rmsVal](const auto& val){return std::abs(val) > 2.5*rmsVal;});
    
    rmsTrunc = std::inner_product(adcLessPedVec.begin(), newEndItr, adcLessPedVec.begin(), 0.);
    numBins  = std::distance(adcLessPedVec.begin(),newEndItr);
    rmsTrunc = std::sqrt(std::max(float(0.),rmsTrunc / float(numBins)));
    
    return;
}
    
}
