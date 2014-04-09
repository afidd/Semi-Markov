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
#ifndef _COW_TOKEN_H_
#define _COW_TOKEN_H_ 1

class Cow
{
public:
  size_t id;
  double birthday;
  int sex;
  int parity;
  double disease_time;

  Cow(size_t id, double time, int sex)
  : id{id}, birthday{time}, sex{sex}, parity{0}, disease_time{0}
  {}

  Cow() {}
};

namespace afidd
{
namespace smv
{
template<>
struct color_type<Cow>
{
  typedef size_t type;
};


template<>
struct unique_color<Cow>
{
  static const bool value=true;
};


size_t color(const Cow& cow)
{
  return cow.id;
}

} // smv
} // end namespace afidd


#endif
