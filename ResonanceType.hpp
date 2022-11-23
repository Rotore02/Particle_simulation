#include "ParticleType.hpp"

#ifndef RESONANCETYPE
#define RESONANCETYPE

class ResonanceType : public ParticleType {
 private:
  double const fWidth;

 protected:
 public:
  ResonanceType(char* Name, double Mass, double Charge, double Width);
  double GetWidth() const override;
  void Print() const override;
};

#endif
