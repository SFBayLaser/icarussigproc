#ifndef ICARUSSIGPROCDEFS_H
#define ICARUSSIGPROCDEFS_H
////////////////////////////////////////////////////////////////////////
//
// File:        ICARUSSigProcDefs.h
//
//              This file provides data structure definitions for filtering
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>

namespace icarussigproc
{
      using VectorShort  = std::vector<short>;
      using VectorInt    = std::vector<int>;
      using VectorFloat  = std::vector<float>;
      using VectorDouble = std::vector<double>;
      using VectorBool   = std::vector<bool>;
      using ArrayShort   = std::vector<VectorShort>;
      using ArrayInt     = std::vector<VectorInt>;
      using ArrayFloat   = std::vector<VectorFloat>;
      using ArrayDouble  = std::vector<VectorDouble>;
      using ArrayBool    = std::vector<VectorBool>;

} // end sigproc namespace
#endif
