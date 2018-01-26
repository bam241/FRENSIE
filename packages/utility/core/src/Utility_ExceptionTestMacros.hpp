//---------------------------------------------------------------------------//
//!
//! \file   Utility_ExceptionTestMacros.hpp
//! \author Alex Robinson
//! \brief  Macros that test if an exception has occurred and throw if so
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_EXCEPTION_TEST_MACROS_HPP
#define UTILITY_EXCEPTION_TEST_MACROS_HPP

// Std Lib Includes
#include <sstream>
#include <string>

/*! Exception test macro used to throw an exception when a required condition
 * fails.
 *
 * This macro is based off of the Teuchos_TestForException macro. This macro
 * should be used anywhere that the failure of a specified conditions
 * warrants the throwing of an exception.
 * \ingroup exception_macros
 */
#define TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)	\
{									\
 const bool throw_exception = (throw_exception_test);			\
 if( throw_exception ){							\
   std::ostringstream detailed_msg;					\
   detailed_msg << "\n" << __FILE__ << ":" << __LINE__ << ":\n"       \
       << "Throw test that evaluated to true: "#throw_exception_test	\
       << "\n" << msg;						\
   const std::string &detailed_msg_str = detailed_msg.str();		\
   throw Exception(detailed_msg_str);					\
 }									\
}

/*! Throw an exception always
 *
 * This macros should be used in conditional execution blocks that should never
 * be reached (e.g. default case statement).
 * \ingroup exception_macros
 */
#define THROW_EXCEPTION( Exception, msg ) \
{					  \
 std::ostringstream detailed_msg;	  \
 detailed_msg << "\n" << __FILE__ << ":" << __LINE__ << ":\n"  \
              << msg;                                            \
 throw Exception(detailed_msg.str());                            \
}

#endif // end UTILITY_EXCEPTION_TEST_MACROS_HPP

//---------------------------------------------------------------------------//
// end Utility_ExceptionTestMacros.hpp
//---------------------------------------------------------------------------//

