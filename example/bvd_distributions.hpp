#ifndef _BVD_DISTRIBUTIONS_HPP_
#define _BVD_DISTRIBUTIONS_HPP_ 1

#include "gspn_random.hpp"

namespace smv=afidd::smv;

template<typename RNG>
class HeiferCalvingDistribution : public smv::TransitionDistribution<RNG> {
  double te_;
  smv::TriangularDistribution<RNG> triangular0_;
  smv::TriangularDistribution<RNG> triangular1_;
  smv::TriangularDistribution<RNG> triangular2_;
 public:
  HeiferCalvingDistribution(double te) : te_(te) {
    double Z1min=17;
    double Z1mid=21.5;
    double Z1max=28;
    double Z2min=28;
    double Z2mid=44;
    double Z2max=57;
    double Z3min=57;
    double Z3mid=65;
    double Z3max=150;
    triangular0_=smv::TriangularDistribution<RNG>(Z1min, Z1mid, Z1max, 0);
    triangular1_=smv::TriangularDistribution<RNG>(Z2min, Z2mid, Z2max, 0);
    triangular2_=smv::TriangularDistribution<RNG>(Z3min, Z3mid, Z3max, 0);
  }

  virtual double Sample(double current_time, RNG& rng) const {
    double X_param=30.0;
    double G_param=282;

    double total=G_param+smv::uniform(rng)*X_param;
    double which_E=smv::uniform(rng);
    if (which_E <= 0.7) {
      total+=0;
    } else if (which_E <= 0.85) {
      total+=triangular0_.Sample(0, rng);
    } else if (which_E <= 0.95) {
      total+=triangular1_.Sample(0, rng);
    } else {
      total+=triangular2_.Sample(0, rng);
    }
    return current_time+total;
  }


  virtual bool BoundedHazard() const { return false; }

  
  virtual double HazardIntegral(double t0, double t1) const { return 0.0; }
  virtual double ImplicitHazardIntegral(double xa, double t0) const {
    return 0.0;
  }

  virtual double EnablingTime() const { return te_; }


  bool CheckSamples(const std::vector<double>& samples, double dt) {
    return false;
  };
};



template<typename RNG>
class DairyCalvingDistribution : public smv::TransitionDistribution<RNG> {
  double te_;
  smv::TriangularDistribution<RNG> triangular0_;
  smv::TriangularDistribution<RNG> triangular1_;
  smv::TriangularDistribution<RNG> triangular2_;
 public:
  DairyCalvingDistribution() {
    double Z1min=17;
    double Z1mid=21.5;
    double Z1max=28;
    double Z2min=28;
    double Z2mid=44;
    double Z2max=57;
    double Z3min=57;
    double Z3mid=65;
    double Z3max=150;
    triangular0_=smv::TriangularDistribution<RNG>(Z1min, Z1mid, Z1max, 0);
    triangular1_=smv::TriangularDistribution<RNG>(Z2min, Z2mid, Z2max, 0);
    triangular2_=smv::TriangularDistribution<RNG>(Z3min, Z3mid, Z3max, 0);
  }

  virtual double Sample(double current_time, RNG& rng) const {
    double X_param=50.0;
    double G_param=282;

    double total=G_param+smv::uniform(rng)*X_param;
    double which_E=smv::uniform(rng);
    if (which_E <= 0.7) {
      total+=0;
    } else if (which_E <= 0.85) {
      total+=triangular0_.Sample(0, rng);
    } else if (which_E <= 0.95) {
      total+=triangular1_.Sample(0, rng);
    } else {
      total+=triangular2_.Sample(0, rng);
    }
    return current_time+total;
  }

  virtual bool BoundedHazard() const { return false; }

  
  virtual double HazardIntegral(double t0, double t1) const { return 0.0; }
  virtual double ImplicitHazardIntegral(double xa, double t0) const {
    return 0.0;
  }

  virtual double EnablingTime() const { return te_; }


  bool CheckSamples(const std::vector<double>& samples, double dt) {
    return false;
  };
};


#endif
