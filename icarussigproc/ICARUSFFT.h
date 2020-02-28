///////////////////////////////////////////////////////////////////////
///
/// \file   ICARUSFFT.h
///
/// \brief  This is meant to provide an interface to fftw
///
/// \author Tracy Usher (usher@SLAC.Stanford.edu)
///
////////////////////////////////////////////////////////////////////////

#ifndef ICARUSFFT_H
#define ICARUSFFT_H

#include <vector>
#include <complex>
#include <algorithm>
#include <stdexcept>

#include "fftw3.h"

namespace icarussigproc
{

template <class T> class ICARUSFFT
{
public:

    using TimeVec      = std::vector<T>;
    using FrequencyVec = std::vector<std::complex<T>>;

    ICARUSFFT(int = 4096);
    ~ICARUSFFT();

    // Basic FFT interface for forward, inverse FFT
    void forwardFFT(TimeVec&, FrequencyVec&) const;
    void inverseFFT(FrequencyVec&, TimeVec&) const;

    // Convolution and Deconvolution
    void convolute(  TimeVec&, const FrequencyVec&, int) const;
    void deconvolute(TimeVec&, const FrequencyVec&, int) const;

    // Recover power spectral density
    void getFFTPower(const TimeVec&, TimeVec&)   const;
    
private:
    // The above are essentially the same operation with different 
    // kernels... so do the common part here
    void convolute(TimeVec&, const FrequencyVec&) const;

    // The plans for running the FFT 
    fftw_plan    fForwardPlan;
    fftw_plan    fInversePlan;

    // Local arrays 
    TimeVec      fTimeVec;
    FrequencyVec fFrequencyVec;
};

template <class T> inline ICARUSFFT<T>::ICARUSFFT(int numTimeSamples) 
{
    // First size our two local vectors
    fTimeVec.resize(numTimeSamples, 0.);
    fFrequencyVec.resize(numTimeSamples, std::complex<T>(0.,0.)); // One extra bin to reflect

    // Now get the plans
    fForwardPlan = fftw_plan_dft_r2c_1d(numTimeSamples, fTimeVec.data(), reinterpret_cast<fftw_complex*>(fFrequencyVec.data()), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    fInversePlan = fftw_plan_dft_c2r_1d(numTimeSamples, reinterpret_cast<fftw_complex*>(fFrequencyVec.data()), fTimeVec.data(), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

    return;
}

template <class T> inline ICARUSFFT<T>::~ICARUSFFT()
{
    fftw_destroy_plan(fForwardPlan);
    fftw_destroy_plan(fInversePlan);

    fTimeVec.clear();
    fFrequencyVec.clear();

    return;
}

template <class T> inline void ICARUSFFT<T>::forwardFFT(TimeVec& timeVec, FrequencyVec& frequencyVec) const
{
    if (timeVec.size() != fTimeVec.size()) 
        throw std::runtime_error("ICARUSFFT: Input time vector size does not match expected");

    frequencyVec.resize(timeVec.size());

    fftw_execute_dft_r2c(fForwardPlan, timeVec.data(), reinterpret_cast<fftw_complex*>(frequencyVec.data()));

    // Reflect the output frequency vector
    size_t vecSize    = timeVec.size();
    size_t nyquistBin = vecSize/2 + 1;

    for(size_t idx = nyquistBin; idx < timeVec.size(); idx++)
        frequencyVec[idx] = std::conj(frequencyVec[vecSize - idx]);

    return;
}

template <class T> inline void ICARUSFFT<T>::inverseFFT(FrequencyVec& frequencyVec, TimeVec& timeVec) const
{
    if (frequencyVec.size() < fFrequencyVec.size()) 
        throw std::runtime_error("ICARUSFFT: Input frequency vector size does not match expected");

    timeVec.resize(frequencyVec.size());

    fftw_execute_dft_c2r(fInversePlan, reinterpret_cast<fftw_complex*>(frequencyVec.data()), timeVec.data());

    // Now normalize
    T normFactor = 1. / T(timeVec.size());

    std::transform(timeVec.begin(),timeVec.end(),timeVec.begin(),std::bind1st(std::multiplies<T>(),normFactor));

    return;
}

template <class T> inline void ICARUSFFT<T>::convolute(TimeVec& timeVec, const FrequencyVec& kernel) const
{
    // Check for consistency
    if (timeVec.size() != fTimeVec.size()) throw std::runtime_error("ICARUSFFT: Input time vector size does not match expected");

    // Get a local frequency vector (we are, after all, a const function)
    FrequencyVec frequencyVec(timeVec.size());

    // Ok, now get teh transform of the input vector
    forwardFFT(timeVec, frequencyVec);

    // Convolute
    std::transform(frequencyVec.begin(),frequencyVec.end(),kernel.begin(),frequencyVec.begin(),std::multiplies<std::complex<T>>());

    // Now transoform back to the time domain
    inverseFFT(frequencyVec, timeVec);

    return;
}

template <class T> inline void ICARUSFFT<T>::convolute(TimeVec& timeVec, const FrequencyVec& kernel, int timeOffset) const
{
    // First do the actual convolution
    convolute(timeVec, kernel);

    // Now rotate the output vector
    std::vector<T> temp;

    if (timeOffset <=0)
    {
        temp.assign(timeVec.begin(),timeVec.begin()-timeOffset);
        timeVec.erase(timeVec.begin(),timeVec.begin()-timeOffset);
        timeVec.insert(timeVec.end(),temp.begin(),temp.end());
    }
    else
    {
        temp.assign(timeVec.end()-timeOffset,timeVec.end());
        timeVec.erase(timeVec.end()-timeOffset,timeVec.end());
        timeVec.insert(timeVec.begin(),temp.begin(),temp.end());
    }

    return;
}

template <class T> inline void ICARUSFFT<T>::deconvolute(TimeVec& timeVec, const FrequencyVec& kernel, int timeOffset) const
{
    // First do the operation
    convolute(timeVec, kernel);

    // Now rotate for the time offset
    std::vector<T> temp;

    if (timeOffset <=0)
    {
        temp.assign(timeVec.end()+timeOffset,timeVec.end());
        timeVec.erase(timeVec.end()+timeOffset,timeVec.end());
        timeVec.insert(timeVec.begin(),temp.begin(),temp.end());
    }
    else
    {
        temp.assign(timeVec.begin(),timeVec.begin()+timeOffset);
        timeVec.erase(timeVec.begin(),timeVec.begin()+timeOffset);
        timeVec.insert(timeVec.end(),temp.begin(),temp.end());
    }

    return;
}

template <typename T> inline void ICARUSFFT<T>::getFFTPower(const TimeVec& inputVec, TimeVec& outputPowerVec) const
{
    // Get a temporary FFT vector
    FrequencyVec fftOutputVec(inputVec.size());

    // Local copy of input vector... sadly
    TimeVec locInputVec = inputVec;

    // Now we compute the convolution kernel which is a straigtforward operation
    forwardFFT(locInputVec, fftOutputVec);
    
    size_t halfFFTDataSize(inputVec.size()/2 + 1);

    // Initialize the output power vec (resize it mostly)
    outputPowerVec.resize(halfFFTDataSize, 0.);
    
    std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTDataSize, outputPowerVec.begin(), [](const auto& val){return std::abs(val);});
    
    return;
}

}

#endif
