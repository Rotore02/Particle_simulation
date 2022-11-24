#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"

void Simulation() {
  char* pPion = new char{'A'};
  char* nPion = new char{'B'};
  char* pKaon = new char{'C'};
  char* nKaon = new char{'D'};
  char* pProton = new char{'E'};
  char* nProton = new char{'F'};
  char* K = new char{'K'};

  Particle::AddParticleType(pPion, 0.13957, 1);
  Particle::AddParticleType(nPion, 0.13957, -1);
  Particle::AddParticleType(pKaon, 0.49367, 1);
  Particle::AddParticleType(nKaon, 0.49367, -1);
  Particle::AddParticleType(pProton, 0.93827, 1);
  Particle::AddParticleType(nProton, 0.93827, -1);
  Particle::AddParticleType(K, 0.89166, 0, 0.050);

  gRandom->SetSeed();

  std::vector<Particle> EventParticles;
  std::vector<Particle> Daughters;
  std::vector<Particle> Generated;

  TH1F* h1 = new TH1F{"h1", "Generated particles distribution", 7, 0, 7};
  TH2F* h2 = new TH2F{
      "h2", "Generated angles distribution", 2000, 0, 2 * M_PI, 1000, 0, M_PI};
  TH1F* h3 = new TH1F{"h3", "Generated impulse distribution", 1000, 0, 3};
  TH1F* h4 =
      new TH1F{"h4", "Generated trasverse impulse distribution", 1000, 0, 3};
  TH1F* h5 = new TH1F{"h5", "Particles energy distribution", 1000, 0, 3};
  TH1F* h6 = new TH1F{"h6", "Invariant Mass", 1000, 0, 3};
  h6->Sumw2();
  TH1F* h7 = new TH1F{"h7", "Invariant mass for particles with opposite charge",
                      1000, 0, 3};
  h7->Sumw2();
  TH1F* h8 = new TH1F{"h8", "Invariant mass for particles with same charge",
                      1000, 0, 3};
  h8->Sumw2();
  TH1F* h9 =
      new TH1F{"h9", "Invariant mass for Pions and Kaons of different charge",
               1000, 0, 3};
  h9->Sumw2();
  TH1F* h10 = new TH1F{
      "h10", "Invariant mass for Pions and Kaons of same charge", 1000, 0, 3};
  h10->Sumw2();
  TH1F* h11 = new TH1F{
      "h11", "Invariant mass for Pions and Kaons derived from a K* particle",
      1000, 0, 3};

  std::cout << "Simulating..." << '\n';

  for (int i = 0; i < 100000; i++) {
    for (int j = 0; j < 100; j++) {
      double x = gRandom->Rndm();
      double phi = gRandom->Uniform(0, 2 * M_PI);
      double theta = gRandom->Uniform(0, M_PI);
      Impulse Imp_;
      double ImpNorm = gRandom->Exp(1);
      Imp_.fPx = ImpNorm * sin(theta) * sin(phi);
      Imp_.fPy = ImpNorm * sin(theta) * cos(phi);
      Imp_.fPz = ImpNorm * cos(theta);
      double Energy;
      if (x < 0.4) {
        Particle Part{*pPion, Imp_};
        EventParticles.push_back(Part);
        Generated.push_back(Part);
        h1->Fill(Part.GetIndex());
        Energy = Part.GetEnergy();
      } else {
        if (x < 0.8) {
          Particle Part{*nPion, Imp_};
          EventParticles.push_back(Part);
          Generated.push_back(Part);
          h1->Fill(Part.GetIndex());
          Energy = Part.GetEnergy();
        } else {
          if (x < 0.85) {
            Particle Part{*pKaon, Imp_};
            EventParticles.push_back(Part);
            Generated.push_back(Part);
            h1->Fill(Part.GetIndex());
            Energy = Part.GetEnergy();
          } else {
            if (x < 0.9) {
              Particle Part{*nKaon, Imp_};
              EventParticles.push_back(Part);
              Generated.push_back(Part);
              h1->Fill(Part.GetIndex());
              Energy = Part.GetEnergy();
            } else {
              if (x < 0.945) {
                Particle Part{*pProton, Imp_};
                EventParticles.push_back(Part);
                Generated.push_back(Part);
                h1->Fill(Part.GetIndex());
                Energy = Part.GetEnergy();
              } else {
                if (x < 0.99) {
                  Particle Part{*nProton, Imp_};
                  EventParticles.push_back(Part);
                  Generated.push_back(Part);
                  h1->Fill(Part.GetIndex());
                  Energy = Part.GetEnergy();
                } else {
                  Particle Part{*K, Imp_};
                  h1->Fill(Part.GetIndex());
                  Energy = Part.GetEnergy();
                  if (x < 0.995) {
                    Particle PosPionDec{*pPion};
                    Particle NegKaonDec{*nKaon};
                    Part.Decay2body(PosPionDec, NegKaonDec);
                    EventParticles.push_back(PosPionDec);
                    EventParticles.push_back(NegKaonDec);
                    Daughters.push_back(PosPionDec);
                    Daughters.push_back(NegKaonDec);
                    h11->Fill(PosPionDec.MassInv(NegKaonDec));
                  } else {
                    Particle NegPionDec{*nPion};
                    Particle PosKaonDec{*pKaon};
                    Part.Decay2body(NegPionDec, PosKaonDec);
                    EventParticles.push_back(NegPionDec);
                    EventParticles.push_back(PosKaonDec);
                    Daughters.push_back(NegPionDec);
                    Daughters.push_back(PosKaonDec);
                    h11->Fill(NegPionDec.MassInv(PosKaonDec));
                  }
                }
              }
            }
          }
        }
      }
      h2->Fill(phi, theta);
      h3->Fill(ImpNorm);
      h4->Fill(sqrt(Imp_.fPx * Imp_.fPx + Imp_.fPy * Imp_.fPy));
      h5->Fill(Energy);
    }

    for (long unsigned int k = 0; k < EventParticles.size(); k++) {
      for (long unsigned int f = k + 1; f < EventParticles.size(); f++) {
        h6->Fill(EventParticles[k].MassInv(EventParticles[f]));
        if ((EventParticles[k].GetCharge()) * (EventParticles[f].GetCharge()) <
            0) {
          h7->Fill(EventParticles[k].MassInv(EventParticles[f]));
        }
        if ((EventParticles[k].GetCharge()) * (EventParticles[f].GetCharge()) >
            0) {
          h8->Fill(EventParticles[k].MassInv(EventParticles[f]));
        }
        if (((EventParticles[k].GetIndex() == 0) &&
             ((EventParticles[f].GetIndex()) == 3)) ||
            (((EventParticles[k].GetIndex()) == 1) &&
             ((EventParticles[f].GetIndex()) == 2))) {
          h9->Fill(EventParticles[k].MassInv(EventParticles[f]));
        }
        if ((((EventParticles[k].GetIndex()) == 0) &&
             ((EventParticles[f].GetIndex()) == 2)) ||
            (((EventParticles[k].GetIndex()) == 1) &&
             ((EventParticles[f].GetIndex()) == 3))) {
          h10->Fill(EventParticles[k].MassInv(EventParticles[f]));
        }
      }
    }
    EventParticles.clear();
    Generated.clear();
    Daughters.clear();
  }

  TFile* Simulation = new TFile{"Simulation.root", "RECREATE"};
  Simulation->Write();
  Simulation->Close();

  std::cout << "Simulation completed. Created file 'Simulation.root'" << '\n';

  delete pPion;
  delete nPion;
  delete pKaon;
  delete nKaon;
  delete pProton;
  delete nProton;
  delete K;
}
