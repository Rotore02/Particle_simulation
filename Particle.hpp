#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "vector"

#ifndef PARTICLE
#define PARTICLE

struct Impulse {
  double fPx;
  double fPy;
  double fPz;
  double Norm();
};

class Particle {
 private:
  Impulse fImp;
  int fIndex;
  static int const fMaxNumParticleType{10};
  static std::vector<ParticleType*> fParticleType;
  static int fNParticleType;
  static int FindParticle(const char name);
  void Boost(double bx, double by, double bz);

 protected:
 public:
  Particle(char name, Impulse Imp = {0, 0, 0});
  Impulse GetImp() const;
  int GetIndex();
  void SetIndex(char Name_);
  void SetImp(Impulse Imp);
  static void AddParticleType(char* Name, double Mass, double Charge);
  static void AddParticleType(char* Name, double Mass, double Charge,
                              double Width);
  double GetCharge() const;
  double GetMass() const;
  static void PrintParticleType();
  void PrintParticle();
  double GetEnergy() const;
  double MassInv(Particle& p) const;
  int Decay2body(Particle& dau1, Particle& dau2) const;
};

#endif
