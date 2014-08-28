#ifndef _BVD_DISTRIBUTIONS_HPP_
#define _BVD_DISTRIBUTIONS_HPP_ 1


template<typename RNG>
class HeiferCalvingDistribution {
  double te_;
  TriangularDistribution triangular0_;
  TriangularDistribution triangular1_;
  TriangularDistribution triangular2_;
 public:
  HeiferCalvingDistribution() {
    Z1min=17;
    Z1mid=21.5;
    Z1max=28;
    Z2min=28;
    Z2min=44;
    Z2max=57;
    Z3min=57;
    Z3mid=65;
    Z3max=150;
    triangular0_=TriangularDistribution(Z1min, Z1mid, Z1max, 0);
    triangular1_=TriangularDistribution(Z2min, Z2mid, Z2max, 0);
    triangular2_=TriangularDistribution(Z3min, Z3mid, Z3max, 0);
  }

  virtual double Sample(double current_time, RNG& rng) const {
    double X_param=30.0;
    double G_param=282;

    double total=G_param+uniform(rng)*X_param;
    double which_E=uniform(rng);
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

  virtual double EnablingTime() const { return te_; }
};



template<typename RNG>
class DairyCalvingDistribution {
  double te_;
  TriangularDistribution triangular0_;
  TriangularDistribution triangular1_;
  TriangularDistribution triangular2_;
 public:
  DairyCalvingDistribution() {
    Z1min=17;
    Z1mid=21.5;
    Z1max=28;
    Z2min=28;
    Z2min=44;
    Z2max=57;
    Z3min=57;
    Z3mid=65;
    Z3max=150;
    triangular0_=TriangularDistribution(Z1min, Z1mid, Z1max, 0);
    triangular1_=TriangularDistribution(Z2min, Z2mid, Z2max, 0);
    triangular2_=TriangularDistribution(Z3min, Z3mid, Z3max, 0);
  }

  virtual double Sample(double current_time, RNG& rng) const {
    double X_param=50.0;
    double G_param=282;

    double total=G_param+uniform(rng)*X_param;
    double which_E=uniform(rng);
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

  virtual double EnablingTime() const { return te_; }
};


#endif
