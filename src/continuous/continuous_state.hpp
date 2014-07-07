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
#ifndef _CONTINUOUS_STATE_H_
#define _CONTINUOUS_STATE_H_ 1

#include <map>
#include "stochnet.hpp"
#include "gspn.hpp"


namespace afidd
{
namespace smv
{
template<typename Mark, typename TransitionKey,
    typename UserState=detail::NoExtraState>
class GSPNState
{
public:
  typedef Mark Marking;
  Marking marking;
  double absolute_time;
  double CurrentTime() const { return absolute_time; }
  void SetTime(double timeprime) { absolute_time=timeprime; }
  UserState user;
  TransitionKey last_transition;
  GSPNState() : absolute_time{0} {}
};

}
}



#endif // _CONTINUOUS_STATE_H_
