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
