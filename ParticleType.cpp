#include "ParticleType.hpp"

#include "iostream"

ParticleType::ParticleType(char* Name, double Mass, double Charge)
    : fName{Name}, fMass{Mass}, fCharge{Charge} {};

char* ParticleType::GetName() const { return fName; }

double ParticleType::GetMass() const { return fMass; }

double ParticleType::GetCharge() const { return fCharge; }

double ParticleType::GetWidth() const { return 0; }

void ParticleType::Print() const {
  std::cout << "particle name: " << *fName << '\n'
            << "particle mass: " << fMass << '\n'
            << "particle charge: " << fCharge << '\n';
}
