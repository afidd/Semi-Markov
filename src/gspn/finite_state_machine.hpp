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
#ifndef _FINITE_STATE_MACHINE_H_
#define _FINITE_STATE_MACHINE_H_ 1

#include <utility>
#include "logging.hpp"


namespace afidd
{
namespace smv
{
template<class Dynamics,class OutputFunction>
class FiniteStateMachine
{
 public:
  using InputToken=typename Dynamics::InputToken;
  using State=typename Dynamics::State;
  using OutputToken=typename OutputFunction::OutputToken;

  //! The dynamics and output function are fixed, and const, for the lifetime.
  FiniteStateMachine(Dynamics& dynamics,
      OutputFunction& output_function)
  : _dynamics(dynamics), _output_function(output_function) {
    BOOST_LOG_TRIVIAL(debug) << "FiniteStateMachine dynamics "
      << &dynamics << std::endl
      << "FiniteStateMachine _dynamics " << &_dynamics;
  }

  ~FiniteStateMachine() {}

  //! Copy constructor uses same dynamics and output function.
  FiniteStateMachine(const FiniteStateMachine& obj)
  : _dynamics(obj._dynamics), _output_function(obj._output_function),
    _state(obj._state)
  {}


  //! Assignment operator just copies the state.
  FiniteStateMachine& operator=(const FiniteStateMachine& other) {
    _state=other._state;
    return *this;
  }


  //! Move constructor only effective if State is MoveConstructible.
  FiniteStateMachine(FiniteStateMachine&& other)
  : _dynamics(other._dynamics), _output_function(other._output_function) {
    _state=std::move(other._state);
  }


  //! Move assignment
  FiniteStateMachine& operator=(FiniteStateMachine&& other) {
    _state=std::move(other._state);
    return *this;
  }


  //! A machine is considered initialized when it receives its initial state.
  void Initialize(State& initial_state) {
    BOOST_LOG_TRIVIAL(debug) << "fsm initialize() begin";
    _state=std::move(initial_state);
    _dynamics.Initialize(_state);
    _initialized=true;
    BOOST_LOG_TRIVIAL(debug) << "fsm initialize() end";
  }


  /*! The output function measures current state before dynamics computes
   *  a new state.
   */
  OutputToken operator()(const InputToken& x) {
    if (!_initialized)
    {
      BOOST_LOG_TRIVIAL(info)
          << "Please call initalize() before calling the FiniteStateMachine";
      throw std::logic_error(
          "Please call initalize() before calling the FiniteStateMachine");
    }
    const auto& y=_output_function(_state);
    _dynamics(_state, x);
    return(y);
  }
 private:
  //! The dynamics determine the next step. They are const.
  Dynamics& _dynamics;
  //! The output function reads the state and token to determine what to return.
  OutputFunction& _output_function;
  //! This is the only mutable part.
  State _state;
  bool _initialized=false;
};
}
}


#endif // _FINITE_STATE_MACHINE_H_
