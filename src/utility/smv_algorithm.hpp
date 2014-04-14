// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//                       Agricultural Research Service
//                  United States Department of Agriculture
//
//   This software/database is a "United States Government Work" under the
//   terms of the United States Copyright Act.  It was written as part of
//   the author's official duties as a United States Government employee
//   and thus cannot be copyrighted.  This software/database is freely
//   available to the public for use. The Department of Agriculture (USDA)
//   and the U.S. Government have not placed any restriction on its use or
//   reproduction.
//
//   Although all reasonable efforts have been taken to ensure the accuracy
//   and reliability of the software and data, the USDA and the U.S.
//   Government do not and cannot warrant the performance or results that
//   may be obtained by using this software or data. The USDA and the U.S.
//   Government disclaim all warranties, express or implied, including
//   warranties of performance, merchantability or fitness for any
//   particular purpose.
//
//   Please cite the author in any work or product based on this material.
//
//   =========================================================================
#ifndef _SMV_ALGORITHM_H_
#define _SMV_ALGORITHM_H_ 1


namespace afidd
{
namespace smv
{


// This helps make a correct less than operator for places and transitions.
// I've gotten this logic wrong too many times before out of laziness. Hence.
template<typename T>
bool LazyLess(const T& a, const T& b)
{
  return (a<b);
}


template<typename T, typename...Ts>
bool LazyLess(const T& a, const T& b, const Ts&... args)
{
  if (a<b)
  {
    return true;
  }
  else if (a==b)
  {
    return LazyLess(args...);
  }
  else
  {
    return false;
  }
}

}
}

#endif
