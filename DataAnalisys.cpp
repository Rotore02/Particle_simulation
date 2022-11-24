#include "Particle.hpp"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "cmath"
#include "fstream"
#include "iostream"
#include "vector"

void SetStyle(std::vector<TH1*>& Histos) {
  TString Labels[7] = {"Generated particles",
                       "Generated azimuthal angles (rad)",
                       "Generated zenithal angles (rad)",
                       "Generated impulse (GeV)",
                       "Generated trasverse impulse (GeV)",
                       "Particles energy (GeV)",
                       "Invariant mass (GeV/(c^2))"};
  int i = 0;
  for (; i < 5; i++) {
    if (i == 0) {
      Histos[i]->GetXaxis()->SetTitle(Labels[0]);
      Histos[i]->GetYaxis()->SetTitle("Occurrencies");
      Histos[i]->GetXaxis()->SetTitleOffset(1);
      Histos[i]->GetYaxis()->SetTitleOffset(1);
      Histos[i]->SetMarkerStyle(1);
      Histos[i]->SetFillColor(30);
      Histos[i]->SetLineStyle(1);
      Histos[i]->SetLineColor(32);
    }
    if (i == 1) {
      Histos[i]->GetXaxis()->SetTitle(Labels[1]);
      Histos[i]->GetYaxis()->SetTitle(Labels[2]);
      Histos[i]->GetZaxis()->SetTitle("Occurrencies");
      Histos[i]->GetXaxis()->SetTitleOffset(1.6);
      Histos[i]->GetYaxis()->SetTitleOffset(1.6);
      Histos[i]->GetZaxis()->SetTitleOffset(0.8);
      Histos[i]->SetMarkerStyle(1);
      Histos[i]->SetFillColor(30);
      Histos[i]->SetLineStyle(1);
      Histos[i]->SetLineColor(32);
    } else {
      if (i < 5) {
        Histos[i]->GetXaxis()->SetTitle(Labels[i + 1]);
        Histos[i]->GetYaxis()->SetTitle("Occurrencies");
        Histos[i]->GetXaxis()->SetTitleOffset(1);
        Histos[i]->GetYaxis()->SetTitleOffset(0.6);
        Histos[i]->SetMarkerStyle(1);
        Histos[i]->SetFillColor(30);
        Histos[i]->SetLineStyle(1);
        Histos[i]->SetLineColor(32);
      }
    }
  }
  for (; i < 11; i++) {
    Histos[i]->GetXaxis()->SetTitle(Labels[6]);
    Histos[i]->GetYaxis()->SetTitle("Occurrencies");
    Histos[i]->GetXaxis()->SetTitleOffset(1);
    Histos[i]->GetYaxis()->SetTitleOffset(1.6);
    Histos[i]->SetMarkerStyle(1);
    Histos[i]->SetFillColor(38);
    Histos[i]->SetLineStyle(1);
    Histos[i]->SetLineColor(37);
  }
  for (; i < 13; i++) {
    Histos[i]->GetXaxis()->SetTitle(Labels[6]);
    Histos[i]->GetYaxis()->SetTitle("Occurrencies");
    Histos[i]->GetXaxis()->SetTitleOffset(1);
    Histos[i]->GetYaxis()->SetTitleOffset(1.6);
    Histos[i]->SetMarkerStyle(1);
    Histos[i]->SetFillColor(42);
    Histos[i]->SetLineStyle(1);
    Histos[i]->SetLineColor(46);
  }
};

void DataAnalisys() {
  std::vector<TH1*> Histos;
  TString s[11] = {"h1", "h2", "h3", "h4",  "h5", "h6",
                   "h7", "h8", "h9", "h10", "h11"};
  TString n[7] = {"Pion+",   "Pion-",   "Kaon+", "Kaon-",
                  "Proton+", "Proton-", "K*"};

  TFile* Simulation = new TFile{"Simulation.root"};
  for (auto i : s) {
    Histos.push_back((TH1*)Simulation->Get(i));
  };

  TF1* LinearX = new TF1{"LinearX", "[0]*x + [1]", 0, 2 * M_PI};
  TF1* LinearY = new TF1{"LinearY", "[0]*x + [1]", 0, M_PI};
  LinearX->SetParameter(0, 0);
  LinearX->SetParameter(1, 5000);
  LinearY->SetParameter(0, 0);
  LinearY->SetParameter(1, 10000);
  TH1D* ProjectionX = ((TH2F*)Simulation->Get("h2"))->ProjectionX();
  TH1D* ProjectionY = ((TH2F*)Simulation->Get("h2"))->ProjectionY();
  ProjectionX->Fit("LinearX");
  ProjectionY->Fit("LinearY");
  TF1* LinearXFit = ProjectionX->GetFunction("LinearX");
  TF1* LinearYFit = ProjectionY->GetFunction("LinearY");

  TF1* Exp = new TF1{"Exp", "[0]*exp([1]*x)", 0, 2};
  Exp->SetParameter(0, 30000);
  Exp->SetParameter(1, -1);
  Histos[2]->Fit("Exp");
  TF1* ExpFit = Histos[2]->GetFunction("Exp");

  TH1F* h7_8 = new TH1F{"h7_8", "h7 minus h8", 1000, 0, 3};
  TH1F* h9_10 = new TH1F{"h9_10", "h9 minus h10", 1000, 0, 3};
  h7_8->Add(Histos[6], Histos[7], 1, -1);
  h9_10->Add(Histos[8], Histos[9], 1, -1);
  Histos.push_back(h7_8);
  Histos.push_back(h9_10);

  TF1* Gaus1 = new TF1{"Gaus1", "gaus", 0.6, 1.2};
  TF1* Gaus2 = new TF1{"Gaus2", "gaus", 0.6, 1.2};
  Gaus1->SetParameter(0, 2500);
  Gaus1->SetParameter(1, 0.89166);
  Gaus1->SetParameter(2, 0.05);
  Gaus2->SetParameter(0, 2500);
  Gaus2->SetParameter(1, 0.89166);
  Gaus2->SetParameter(2, 0.05);
  h7_8->Fit(Gaus1, "R");
  h9_10->Fit(Gaus2, "R");

  std::ofstream os{"Data.txt"};

  if (!os) {
    std::cout << "couldn't find 'Data.txt' " << '\n';
  };

  if (os.good()) {
    for (int i = 0; i < 13; i++) {
      os << "Histo " << i + 1 << "entries: " << Histos[i]->GetEntries() << '\n'
         << "-------------------------------" << '\n'
         << '\n';
    }
    for (int i = 0; i < 7; i++) {
      os << "Histo 1 generated particles:" << '\n'
         << n[i] << " generated: " << Histos[0]->GetBinContent(i + 1) << " +- "
         << Histos[0]->GetBinError(i + 1) << '\n'
         << "-----------------------------" << '\n'
         << '\n';
    }
    os << "Histo 2 linear fit (A*x + B):" << '\n'
       << "Proj x:" << '\n'
       << "A = " << LinearXFit->GetParameter(0) << " +- "
       << LinearXFit->GetParError(0) << " , "
       << "B = " << LinearXFit->GetParameter(1) << " +- "
       << LinearXFit->GetParError(1) << '\n'
       << "ChiSquare/NDF = "
       << (LinearXFit->GetChisquare()) / (LinearXFit->GetNDF()) << " , "
       << "Fit probability: " << LinearXFit->GetProb() << '\n'
       << '\n'
       << "Proj y:" << '\n'
       << "A = " << LinearYFit->GetParameter(0) << " +- "
       << LinearYFit->GetParError(0) << " , "
       << "B = " << LinearYFit->GetParameter(1) << " +- "
       << LinearYFit->GetParError(1) << '\n'
       << "ChiSquare/NDF = "
       << (LinearYFit->GetChisquare()) / (LinearYFit->GetNDF()) << " , "
       << "Fit probability: " << LinearYFit->GetProb() << '\n'
       << "-------------------------------" << '\n'
       << '\n'
       << "Histo 3 exponential fit (A*exp(B*x)):" << '\n'
       << "A = " << ExpFit->GetParameter(0) << " +- " << ExpFit->GetParError(0)
       << " , "
       << "B = " << ExpFit->GetParameter(1) << " +- " << ExpFit->GetParError(1)
       << '\n'
       << "ChiSquare/NDF = " << (ExpFit->GetChisquare()) / (ExpFit->GetNDF())
       << " , "
       << "Fit probability: " << ExpFit->GetProb() << '\n'
       << "Fit mean: " << 1 / ExpFit->GetParameter(1) << " +- "
       << sqrt((1 / (ExpFit->GetParameter(1) * ExpFit->GetParameter(1))) *
               ExpFit->GetParError(1))
       << '\n'
       << "---------------------------------" << '\n'
       << '\n'
       << "Histo 7_8 gaussian fit (A*exp(-0.5*((x-B)/C)*((x-B)/C)):" << '\n'
       << "A = " << Gaus1->GetParameter(0) << " +- " << Gaus1->GetParError(0)
       << " , "
       << "B (K* mass) = " << Gaus1->GetParameter(1) << " +- "
       << Gaus1->GetParError(1) << " , "
       << "C (K* width) = " << Gaus1->GetParameter(2) << " +- "
       << Gaus1->GetParError(2) << '\n'
       << "ChiSquare/NDF = " << (Gaus1->GetChisquare()) / (Gaus1->GetNDF())
       << " , "
       << "Fit Probability = " << Gaus1->GetProb() << '\n'
       << '\n'
       << "Histo 9_10 gaussian fit (A*exp(-0.5*((x-B)/C)*((x-B)/C)):" << '\n'
       << "A = " << Gaus2->GetParameter(0) << " +- " << Gaus2->GetParError(0)
       << " , "
       << "B (K* mass) = " << Gaus1->GetParameter(1) << " +- "
       << Gaus2->GetParError(1) << " , "
       << "C (K* width) = " << Gaus2->GetParameter(2) << " +- "
       << Gaus2->GetParError(2) << '\n'
       << "ChiSquare/NDF = " << (Gaus2->GetChisquare()) / (Gaus2->GetNDF())
       << " , "
       << "Fit Probability = " << Gaus2->GetProb() << '\n';
  };

  TCanvas* c1 = new TCanvas{"c1", "Particles properties c1", 1000, 700};
  TCanvas* c2 = new TCanvas{"c2", "Particles properties c2", 1000, 700};
  TCanvas* c3 = new TCanvas{"c3", "Invariant mass c1", 1000, 700};
  TCanvas* c4 = new TCanvas{"c4", "Invariant mass c2", 1000, 700};
  SetStyle(Histos);
  c1->Divide(1, 2);
  c2->Divide(1, 3);
  c3->Divide(2, 2);
  c4->Divide(2, 2);
  c1->cd(1);
  Histos[0]->Draw("HIST");
  c1->cd(2);
  Histos[1]->Draw("LEGO");
  for (int i = 2; i < 13; i++) {
    if (i < 5) {
      c2->cd(i - 1);
      Histos[i]->Draw("HIST");
    } else {
      if (i < 9) {
        c3->cd(i - 4);
        Histos[i]->Draw("HIST");
      } else {
        if (i < 11) {
          c4->cd(i - 8);
          Histos[i]->Draw("HIST");
        } else {
          c4->cd(3);
          Histos[11]->Draw();
          c4->cd(4);
          Histos[12]->Draw();
        }
      }
    }
  }
  c1->Print("ParticlesProperties_c1.pdf");
  c2->Print("ParticlesProperties_c2.pdf");
  c3->Print("InvariantMass_c1.pdf");
  c4->Print("InvariantMass_c2.pdf");
}
