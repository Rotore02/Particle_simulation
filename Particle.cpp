#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

Particle::Particle(char name, Impulse Imp) : fImp{Imp} {
  fIndex = FindParticle(name);
}

std::vector<ParticleType*> Particle::fParticleType;

int Particle::fNParticleType{};

Impulse Particle::GetImp() const { return fImp; };

int Particle::GetIndex() { return fIndex; };

void Particle::SetIndex(char Name_) { fIndex = FindParticle(Name_); }

void Particle::SetImp(Impulse Imp) {
  fImp.fPx = Imp.fPx;
  fImp.fPy = Imp.fPy;
  fImp.fPz = Imp.fPz;
}

int Particle::FindParticle(const char name) {
  auto i = find_if(fParticleType.begin(), fParticleType.end(),
                   [&](auto p) { return *p->GetName() == name; });
  if (i == fParticleType.end()) {
    return 11;
  } else {
    return std::distance(fParticleType.begin(), i);
  }
}

void Particle::AddParticleType(char* Name, double Mass, double Charge) {
  ParticleType* PartTp = new ParticleType{Name, Mass, Charge};
  fNParticleType = fParticleType.size();
  int Index = FindParticle(*Name);
  if (Index == 11) {
    fParticleType.push_back(PartTp);
  }
  if (Index < 11) {
    std::cout << "particle " << *PartTp->GetName() << " already exists" << '\n';
  }
  if (fParticleType.size() > fMaxNumParticleType) {
    std::cout << "number of particle exceeded" << '\n';
  }
};
void Particle::AddParticleType(char* Name, double Mass, double Charge,
                               double Width) {
  ResonanceType* PartTp = new ResonanceType{Name, Mass, Charge, Width};
  fNParticleType = fParticleType.size();
  int Index = FindParticle(*Name);
  if (Index == 11) {
    fParticleType.push_back(PartTp);
  }
  if (Index < fParticleType.size()) {
    std::cout << "particle " << *PartTp->GetName() << " already exists" << '\n';
  }
  if (fParticleType.size() > fMaxNumParticleType) {
    std::cout << "number of particle exceeded" << '\n';
  }
};

double Particle::GetCharge() const {
  return fParticleType[fIndex]->GetCharge();
};

double Particle::GetMass() const { return fParticleType[fIndex]->GetMass(); }

void Particle::PrintParticleType() {
  for (; fParticleType.begin() != fParticleType.end();
       fParticleType.begin()++) {
    (*fParticleType.begin())->Print();
  }
}

double Particle::GetEnergy() const {
  return sqrt(
      (GetMass() * GetMass()) +
      (fImp.fPx * fImp.fPx + fImp.fPy * fImp.fPy + fImp.fPz * fImp.fPz));
}

double Particle::MassInv(Particle& p) const {
  Impulse Imp{fImp.fPx + p.GetImp().fPx, fImp.fPy + p.GetImp().fPy,
              fImp.fPz + p.GetImp().fPz};
  return sqrt((GetEnergy() + p.GetEnergy()) * (GetEnergy() + p.GetEnergy()) -
              (Imp.fPx * Imp.fPx + Imp.fPy * Imp.fPy + Imp.fPz * Imp.fPz));
}

int Particle::Decay2body(Particle& dau1, Particle& dau2) const {
  if (GetMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIndex > -1) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += fParticleType[fIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf(
        "Decayment cannot be preformed because mass is too low in this "
        "channel\n");
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetImp(Impulse{pout * sin(theta) * cos(phi),
                      pout * sin(theta) * sin(phi), pout * cos(theta)});
  dau2.SetImp(Impulse{-pout * sin(theta) * cos(phi),
                      -pout * sin(theta) * sin(phi), -pout * cos(theta)});

  double energy = sqrt(fImp.fPx * fImp.fPx + fImp.fPy * fImp.fPy +
                       fImp.fPz * fImp.fPz + massMot * massMot);

  double bx = fImp.fPx / energy;
  double by = fImp.fPy / energy;
  double bz = fImp.fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz) {
  double energy = GetEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fImp.fPx + by * fImp.fPy + bz * fImp.fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fImp.fPx += gamma2 * bp * bx + gamma * bx * energy;
  fImp.fPy += gamma2 * bp * by + gamma * by * energy;
  fImp.fPz += gamma2 * bp * bz + gamma * bz * energy;
}

double Impulse::Norm() { return fPx * fPx + fPy * fPy + fPz * fPz; }
