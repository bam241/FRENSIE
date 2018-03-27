//---------------------------------------------------------------------------//
//!
//! \file   Utility.__init__.i
//! \author Alex Robinson
//! \brief  The utility module swig interface file
//!
//---------------------------------------------------------------------------//

%define %utility_docstring
"
PyFrensie.Utility is the python interface to the FRENSIE utilities package.
"
%enddef

%module(package   = "PyFrensie.Utility",
        autodoc   = "1",
        docstring = %utility_docstring) __init__

%pythonbegin
%{
  # With python 2.7.14 and swig version 3.0.10 the default generated import
  # code cannot find the ___init__.so file. Unfortunately the 'moduleimport'
  # option of the %module macro does not seem to work with this version
  # of swig either. To get the import code working we have to manually add
  # the current directory to the system path temporarily.

  import os.path
  import sys
  current_dir,file_name = os.path.split(__file__)
  sys.path.insert(0, current_dir)
%}

// Set the special python variables
%pythoncode
%{
# Remove the local current directory from the sys path (added to help
# import code - see comment above)
sys.path.pop(0)

__all__ = ['Distribution',
           'Prng',
           'Interpolation'
           ]
%}

%{
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "Utility_RandomNumberGenerator.hpp"
%}

// Add the shortcut for initializing the random number generator
%feature("autodoc")
initFrensiePrng
"
This method can be used to initialize the Utility.Prng.RandomNumberGenerator
instead of calling 'Utility.Prng.RandomNumberGenerator.createStreams()'.
"

%inline %{
//! Initialize the random number generator
void initFrensiePrng()
{
  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}
%}

//---------------------------------------------------------------------------//
// end Utility.__init__.i
//---------------------------------------------------------------------------//
