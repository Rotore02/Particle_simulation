#include "string"

#ifndef PARTICLETYPE
#define PARTICLETYPE

class ParticleType {
 private:
  double const fMass;
  double const fCharge;
  char* const fName;

 protected:
 public:
  ParticleType(char* Name, double Mass, double Charge);
  char* GetName() const;
  double GetCharge() const;
  double GetMass() const;
  virtual double GetWidth() const;
  virtual void Print() const;
};

#endif
