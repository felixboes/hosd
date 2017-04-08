/// @addtogroup struct
/// @{

/////////////////////////////////////////////////////////////////////////////
///
/// @file inthash.h
///
/// This file contains the definition of two hashing methods for integers.
/// Such methods are necessary when using a hashedset container.
///
/// @author Pawel Pilarczyk
///
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 1997-2013 by Pawel Pilarczyk.
//
// This file is part of the Homology Library.  This library is free software;
// you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this software; see the file "license.txt".  If not, write to the
// Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA.

// Started in January 2002. Last revision: October 24, 2013.


#ifndef _CHOMP_STRUCT_INTHASH_H_
#define _CHOMP_STRUCT_INTHASH_H_


//namespace chomp {
//namespace homology {

// --------------------------------------------------
// ------------- hash keys for integers -------------
// --------------------------------------------------

/// The first hash key for an unsigned int number.
inline int_t hashkey1 (const unsigned long &number)
{
	return static_cast<int_t> (number);
} /* hashkey1 */

/// The second hash key for an unsigned int number.
inline int_t hashkey2 (const unsigned long &number)
{
	return static_cast<int_t> (number ^ 0xFA5A75A7ul) << 8;
} /* hashkey2 */

/// This macro is used to define a hash key for any type that can be
/// cast onto an unsigned inteter type, in particular, in the library
/// it is used to define the hash key functions for int, long, and char,
/// both signed and unsigned.
#define DEFHASHKEYS(type) \
inline int_t hashkey1 (const type &number) \
{ return hashkey1 (static_cast<const unsigned long> (number)); } \
inline int_t hashkey2 (const type &number) \
{ return hashkey2 (static_cast<const unsigned long> (number)); }

DEFHASHKEYS(unsigned int)
DEFHASHKEYS(signed int)
//DEFHASHKEYS(unsigned long)
DEFHASHKEYS(signed long)
DEFHASHKEYS(unsigned short)
DEFHASHKEYS(signed short)
DEFHASHKEYS(unsigned char)
DEFHASHKEYS(signed char)

#undef DEFHASHKEYS

/// A template function that extracts the first hash key from an object
/// which has the method "hashkey1".
/// Provided for backwards compatibility with some data types.
//template <class T>
//inline int_t hashkey1 (const T &x)
//{
//	return x. hashkey1 ();
//}

/// A template function that extracts the second hash key from an object
/// which has the method "hashkey2".
/// Provided for backwards compatibility with some data types.
//template <class T>
//inline int_t hashkey2 (const T &x)
//{
//	return x. hashkey2 ();
//}

//} // namespace homology
//} // namespace chomp

#endif // _CHOMP_STRUCT_INTHASH_H_

/// @}

