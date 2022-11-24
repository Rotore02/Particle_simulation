#include "ResonanceType.hpp"

#include "iostream"

ResonanceType::ResonanceType(char* Name, double Mass, double Charge,
                             double Width)
    : ParticleType{Name, Mass, Charge}, fWidth{Width} {};

double ResonanceType::GetWidth() const { return fWidth; }

void ResonanceType::Print() const {
  Print();
  std::cout << "particle width: " << fWidth << '\n';
}
