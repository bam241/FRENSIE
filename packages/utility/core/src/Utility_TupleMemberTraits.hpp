//---------------------------------------------------------------------------//
//!
//! \file   Utility_TupleMemberTraits.hpp
//! \author Alex Robinson
//! \brief  Tuple member traits class specializations
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_TUPLE_MEMBER_TRAITS_HPP
#define UTILITY_TUPLE_MEMBER_TRAITS_HPP

// FRENSIE Includes
#include "Utility_TupleMemberTraitsDecl.hpp"
#include "Utility_Tuple.hpp"

namespace Utility{

namespace Traits{

/*! The partial specialization of the TupleMemberTraits for FIRST
 * \ingroup tuple_member_traits
 */
template<typename Tuple>
struct TupleMemberTraits<Tuple,FIRST>
{
  typedef typename Tuple::firstType tupleMemberType;
  static inline tupleMemberType get( const Tuple &tuple )
  { return tuple.first; }
  static inline void set( Tuple &tuple, const tupleMemberType &value )
  { tuple.first = value; }
};

/*! The partial specialization of the TupleMemberTraits for SECOND
 * \ingroup tuple_member_traits
 */
template<typename Tuple>
struct TupleMemberTraits<Tuple,SECOND>
{
  typedef typename Tuple::secondType tupleMemberType;
  static inline tupleMemberType get( const Tuple &tuple )
  { return tuple.second; }
  static inline void set( Tuple &tuple, const tupleMemberType &value )
  { tuple.second = value; }
};

/*! The partial specialization of the TupleMemberTraits for THIRD
 * \ingroup tuple_member_traits
 */
template<typename Tuple>
struct TupleMemberTraits<Tuple,THIRD>
{
  typedef typename Tuple::thirdType tupleMemberType;
  static inline tupleMemberType get( const Tuple &tuple )
  { return tuple.third; }
  static inline void set( Tuple &tuple, const tupleMemberType &value )
  { tuple.third = value; }
};

/*! The partial specialization of the TupleMemberTraits for FOURTH
 * \ingroup tuple_member_traits
 */
template<typename Tuple>
struct TupleMemberTraits<Tuple,FOURTH>
{
  typedef typename Tuple::fourthType tupleMemberType;
  static inline tupleMemberType get( const Tuple &tuple )
  { return tuple.fourth; }
  static inline void set( Tuple &tuple, const tupleMemberType &value )
  { tuple.fourth = value; }
};

/*! The specialization of the TupleMemberTraits for double and FIRST
 * \ingroup tuple_member_traits
 */
template<>
struct TupleMemberTraits<double,FIRST>
{
  typedef double tupleMemberType;
  static inline double get( const double value )
  { return value; }
  static inline void set( double &value, const double &new_value )
  { value = new_value; }
};

/*! The specialization of the TupleMemberTraits for const double and FIRST
 * \ingroup tuple_member_traits
 */
template<>
struct TupleMemberTraits<const double,FIRST>
{
  typedef const double tupleMemberType;
  static inline const double get( const double value )
  { return value; }
  static inline void set( const double &value, const double &new_value )
  { /* ... */ }
};

/*! The specialization of the TupleMemberTraits for unsigned and FIRST
 * \ingroup tuple_member_traits
 */
template<>
struct TupleMemberTraits<unsigned,FIRST>
{
  typedef unsigned tupleMemberType;
  static inline unsigned get( const unsigned value )
  { return value; }
  static inline void set( unsigned &value, const unsigned &new_value )
  { value = new_value; }
};

/*! The specialization of the TupleMemberTraits for const unsigned and FIRST
 * \ingroup tuple_member_traits
 */
template<>
struct TupleMemberTraits<const unsigned,FIRST>
{
  typedef unsigned tupleMemberType;
  static inline unsigned get( const unsigned value )
  { return value; }
  static inline void set( const unsigned &value, const unsigned &new_value )
  { /* ... */ }
};

} // end Traits namespace

} // end Utility namespace

#endif // end UTILITY_TUPLE_MEMBER_TRAITS_HPP

//---------------------------------------------------------------------------//
// end Utility_TupleMemberTraits.hpp
//---------------------------------------------------------------------------//

