#ifndef _SMV_ALGORITHM_H_
#define _SMV_ALGORITHM_H_ 1


namespace afidd
{
namespace smv
{


// This helps make a correct less than operator for places and transitions.
// I've gotten this logic wrong too many times before out of laziness. Hence.
template<typename T>
bool lazy_less(const T& a, const T& b)
{
  return (a<b);
}


template<typename T, typename...Ts>
bool lazy_less(const T& a, const T& b, const Ts&... args)
{
  if (a<b)
  {
    return true;
  }
  else if (a==b)
  {
    return lazy_less(args...);
  }
  else
  {
    return false;
  }
}

}
}

#endif
