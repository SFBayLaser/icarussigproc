/**
 * \file Deconvolve.h
 *
 * \ingroup icarussigproc
 * 
 * \brief Class def header for a class Deconvolve
 *
 * @author koh0207
 */

/** \addtogroup icarussigproc

    @{*/
#ifndef __SIGPROC_TOOLS_DECONVOLVE_H__
#define __SIGPROC_TOOLS_DECONVOLVE_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "MiscUtils.h"

namespace icarussigproc {

  /**
     \class Deconvolve
     User defined class Deconvolve ... these comments are used to generate
     doxygen documentation!
  */
  class Deconvolve{
    
    public:
      
      /// Default constructor
      Deconvolve(){}

      void filterLee(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void filterLee(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void filterLee(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int
      );


      void MMWF(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void MMWF(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void MMWF(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int
      );


      void MMWFStar(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const unsigned int,
        const unsigned int
      );

      void MMWFStar(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const unsigned int,
        const unsigned int
      );

      void MMWFStar(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const unsigned int,
        const unsigned int
      );


      void filterLeeEnhanced(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );

      void filterLeeEnhanced(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );

      void filterLeeEnhanced(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );
      
      /// Default destructor
      ~Deconvolve(){}

    private:

      template <typename T>
      void filterLee(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=7,
        const unsigned int sy=7);

      template <typename T>
      void MMWF(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=7,
        const unsigned int sy=7);

      template <typename T>
      void MMWFStar(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const unsigned int sx=7,
        const unsigned int sy=7);
    
      template <typename T>
      void filterLeeEnhanced(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=7,
        const unsigned int sy=7,
        const float a=1,
        const float epsilon=2.5);

      template <typename T>
      void inverseFilter(
        std::vector<std::vector<T>>& deconvolvedWaveform,
        const std::vector<std::vector<T>>& waveLessCoherent,
        const std::vector<std::vector<T>>& responseFunction
      );
  };
}

#endif
/** @} */ // end of doxygen group 
