//---------------------------------------------------------------------------//
//!
//! \file   Utility_ArrayView.hpp
//! \author Alex Robinson
//! \brief  The array view class declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_ARRAY_VIEW_HPP
#define UTILITY_ARRAY_VIEW_HPP

// Std Lib Includes
#include <vector>
#include <array>

// FRENSIE Includes
#include "Utility_View.hpp"

namespace Utility{

/*! The array view class
 * 
 * This class was inspired by the Teuchos::ArrayView class found in
 * the Trilinos Software package 
 * (https://trilinos.org/docs/dev/packages/teuchos/doc/html/index.html).
 */
template<typename T>
class ArrayView : public View<T*>
{

public:

  //! Default constructor
  ArrayView();

  //! Iterator constructor
  ArrayView( T* start, T* end );

  //! Range constructor
  ArrayView( T* array_start, const typename ArrayView<T>::size_type array_size );

  //! Vector constructor
  ArrayView( std::vector<T>& vector );

  //! Const vector constructor
  template<typename U>
  ArrayView( const std::vector<U>& vector );

  //! Array constructor
  template<size_t N>
  ArrayView( std::array<T,N>& array );

  //! Const array constructor
  template<typename U,size_t N>
  ArrayView( const std::array<U,N>& array );

  //! Copy constructor
  ArrayView( ArrayView<T>& other_view );

  //! Const array view copy constructor
  template<typename U>
  ArrayView( const ArrayView<U>& other_view );

  //! Assignment operator
  ArrayView<T>& operator=( ArrayView<T>& other_view );

  //! Const view assignment operator
  template<typename U>
  ArrayView<T>& operator=( const ArrayView<U>& other_view );

  //! Destructor
  ~ArrayView()
  { /* ... */ }

  //! Return a sub-array view
  ArrayView<T> operator()( const typename ArrayView<T>::size_type offset,
                           const typename ArrayView<T>::size_type size ) const;

  //! Return a const array view
  ArrayView<const typename std::remove_const<T>::type> toConst() const;

  //! Implicitly convert to a const array view
  operator ArrayView<const typename std::remove_const<T>::type>() const;
};

//! Create an array view of a std::vector
template<typename T>
inline ArrayView<T> arrayView( std::vector<T>& vector )
{
  return ArrayView<T>( vector );
}

//! Create a const array view of a std::vector
template<typename T>
inline ArrayView<const T> arrayView( const std::vector<T>& vector )
{
  return ArrayView<const T>( vector );
}

//! Create a const array view of a std::vector
template<typename T>
inline ArrayView<const T> arrayViewOfConst( const std::vector<T>& vector )
{
  return ArrayView<const T>( vector );
}

//! Create an array view of a std::array
template<typename T, size_t N>
inline ArrayView<T> arrayView( std::array<T,N>& array )
{
  return ArrayView<T>( array );
}

//! Create a const array view of a std::array
template<typename T, size_t N>
inline ArrayView<const T> arrayView( const std::array<T,N>& array )
{
  return ArrayView<const T>( array );
}

//! Create a const array view of a std::array
template<typename T, size_t N>
inline ArrayView<const T> arrayViewOfConst( const std::array<T,N>& array )
{
  return ArrayView<const T>( array );
}

//! Const cast array view
template<typename T2, typename T1>
inline ArrayView<T2> av_const_cast( const ArrayView<T1>& array_view )
{
  return ArrayView<T2>( const_cast<T2*>( array_view.begin() ),
                        array_view.size() );
}

//! Reinterpret cast array view
template<typename T2, typename T1>
inline ArrayView<T2> av_reinterpret_cast( const ArrayView<T1>& array_view )
{
  return ArrayView<T2>( reinterpret_cast<T2*>( const_cast<ArrayView<T1>&>(array_view).begin() ),
                        reinterpret_cast<T2*>( const_cast<ArrayView<T1>&>(array_view).end() ) );
}
  
} // end Utility namespace

//---------------------------------------------------------------------------//
// Template includes
//---------------------------------------------------------------------------//

#include "Utility_ArrayView_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_ARRAY_VIEW_HPP

//---------------------------------------------------------------------------//
// end Utility_ArrayView.hpp
//---------------------------------------------------------------------------//
