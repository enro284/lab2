#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TMath.h"
#include "TRandom.h"
#include "particle.hpp"

#define NBINS 100

void Main() {
  Particle::AddParticleType("pione+", 0.13957, 1.);
  Particle::AddParticleType("pione-", 0.13957, -1.);
  Particle::AddParticleType("kaone+", 0.49367, 1.);
  Particle::AddParticleType("kaone-", 0.49367, -1.);
  Particle::AddParticleType("protone+", 0.93827, 1.);
  Particle::AddParticleType("protone-", 0.93827, -1.);
  Particle::AddParticleType("k*", 0.89166, 0., 0.05);

  gRandom->SetSeed();

  const int nEvents = 10000;
  const int nParticlePerEvents = 100;

  Particle *EventParticles[120];
  for (int i{0}; i < 120; ++i) {
    EventParticles[i] = new Particle("pione+");
  }

  auto hParticleType = new TH1I("hParticleType", "tipi di particelle", 7, 0, 7);

  auto hDistributionTheta =
      new TH1F("hDistributionTheta", "distribuzione angoli theta", NBINS, 0.,
               TMath::Pi());
  auto hDistributionPhi =
      new TH1F("hDistributionPHI", "distribuzione angoli phi", NBINS, 0.,
               2 * TMath::Pi());
  auto hImpulse = new TH1F("hImpulse", "impulso", NBINS, 0., 8.);
  auto hImpulseT = new TH1F("hImpulseT", "impulso trasverso", NBINS, 0., 8.);
  auto hEnergy =
      new TH1F("hEnergy", "istogramma energia della particella", NBINS, 0., 8.);
  auto hInvMass = new TH1F("hInvMass", "massa invariante", NBINS, 0., 8.);

  auto hInvMassDisc = new TH1F("hInvMassDisc",
                               "massa invariante in combinazione "
                               "con carica di segno discorde",
                               NBINS, 0., 8.);
  hInvMassDisc->Sumw2();
  auto hInvMassConc = new TH1F("hInvMassConc",
                               "massa invariante in combinazione "
                               "con carica di segno concorde",
                               NBINS, 0., 8.);
  hInvMassConc->Sumw2();
  auto hInvMassComb1 = new TH1F("hInvMassComb1",
                                "massa invariante con combinazioni "
                                "di tipo pione+/Kaone- e pione-/kaone+",
                                NBINS, 0., 8.);
  // hInvMassDisc->Sumw2();
  auto hInvMassComb2 = new TH1F("hInvMassComb2",
                                "massa invariante con combinazioni "
                                "di tipo pione+/Kaone+ e pione-/kaone-",
                                NBINS, 0., 8.);
  // hInvMassConc->Sumw2();
  auto hInvMassK = new TH1F("hInvMassK",
                            "istogramma massa invariante fra le particelle che "
                            "derivano dal decadimento di k*",
                            NBINS, 0.5, 1.3);

  for (int event{0}; event < nEvents; ++event) {
    double x;
    int nFiglie{0};
    for (int j{0}; j < nParticlePerEvents; ++j) {
      x = gRandom->Rndm();

      if (x < 0.4) {
        EventParticles[j]->SetIndex("pione+");
        hParticleType->Fill(0.5);
      } else if (x < 0.8) {
        EventParticles[j]->SetIndex("pione-");
        hParticleType->Fill(1.5);
      } else if (x < 0.85) {
        EventParticles[j]->SetIndex("kaone+");
        hParticleType->Fill(2.5);
      } else if (x < 0.9) {
        EventParticles[j]->SetIndex("kaone-");
        hParticleType->Fill(3.5);
      } else if (x < 0.945) {
        EventParticles[j]->SetIndex("protone+");
        hParticleType->Fill(4.5);
      } else if (x < 0.99) {
        EventParticles[j]->SetIndex("protone-");
        hParticleType->Fill(5.5);
      } else {
        EventParticles[j]->SetIndex("k*");
        hParticleType->Fill(6.5);
      }

      double phi = gRandom->Rndm() * 2 * TMath::Pi();
      hDistributionPhi->Fill(phi);
      double theta = gRandom->Rndm() * TMath::Pi();
      hDistributionTheta->Fill(theta);
      double p = gRandom->Exp(1.);
      hImpulse->Fill(p);

      double pX = p * sin(theta) * cos(phi);
      double pY = p * sin(theta) * sin(phi);
      double pZ = p * cos(theta);

      hImpulseT->Fill(sqrt(pX * pX + p * p));

      double energy = EventParticles[j]->EnergyTot();

      hEnergy->Fill(energy);

      EventParticles[j]->SetP(pX, pY, pZ);

      if (EventParticles[j]->GetIndex() == 6) {
        EventParticles[j]->Decay2body(
            *EventParticles[nParticlePerEvents + nFiglie],
            *EventParticles[nParticlePerEvents + nFiglie + 1]);

        if (gRandom->Rndm() < 0.5) {
          EventParticles[nParticlePerEvents + nFiglie]->SetIndex("pione+");
          EventParticles[nParticlePerEvents + nFiglie + 1]->SetIndex("kaone-");
        } else {
          EventParticles[nParticlePerEvents + nFiglie]->SetIndex("pione-");
          EventParticles[nParticlePerEvents + nFiglie + 1]->SetIndex("kaone+");
        }
        nFiglie += 2;
      }
    }
    // invariant mass
    for (int i{0}; i < nParticlePerEvents; ++i) {
      for (int j{0}; j < nParticlePerEvents; ++j) {
        if (i != j) {
          hInvMass->Fill(EventParticles[i]->InvMass(*EventParticles[j]));
        }
      }
    }
    // invariant mass discord
    for (int i{0}; i < nParticlePerEvents + nFiglie; ++i) {
      for (int j{0}; j < nParticlePerEvents + nFiglie; ++j) {
        if (i != j &&
            EventParticles[i]->GetCharge() != EventParticles[j]->GetCharge()) {
          hInvMassDisc->Fill(EventParticles[i]->InvMass(*EventParticles[j]));
        }
      }
    }
    // invariant mass concord
    for (int i{0}; i < nParticlePerEvents + nFiglie; ++i) {
      for (int j{0}; j < nParticlePerEvents + nFiglie; ++j) {
        if (i != j &&
            EventParticles[i]->GetCharge() == EventParticles[j]->GetCharge()) {
          hInvMassConc->Fill(EventParticles[i]->InvMass(*EventParticles[j]));
        }
      }
    }
    // massa invariante fra tutte le particelle con p+/k- e p-/k+
    for (int i{0}; i < nParticlePerEvents + nFiglie; ++i) {
      for (int j{0}; j < nParticlePerEvents + nFiglie; ++j) {
        if (i != j && ((EventParticles[i]->GetIndex() == 0 &&
                        EventParticles[j]->GetIndex() == 3) ||
                       (EventParticles[i]->GetIndex() == 1 &&
                        EventParticles[j]->GetIndex() == 2))) {
          hInvMassComb1->Fill(EventParticles[i]->InvMass(*EventParticles[j]));
        }
      }
    }
    // massa invariante fra tutte le particelle con combinazione p+/k+ e p-/k-
    for (int i{0}; i < nParticlePerEvents + nFiglie; ++i) {
      for (int j{0}; j < nParticlePerEvents + nFiglie; ++j) {
        if (i != j && ((EventParticles[i]->GetIndex() == 0 &&
                        EventParticles[j]->GetIndex() == 2) ||
                       (EventParticles[i]->GetIndex() == 1 &&
                        EventParticles[j]->GetIndex() == 3))) {
          hInvMassComb2->Fill(EventParticles[i]->InvMass(*EventParticles[j]));
        }
      }
    }
    // massa invariante fra le particelle generate che derivano dal decadimento
    // di k*
    for (int i{nParticlePerEvents}; i < nParticlePerEvents + nFiglie - 1;
         i += 2) {
      hInvMassK->Fill(EventParticles[i]->InvMass(*EventParticles[i + 1]));
    }
  }

  auto c1 = new TCanvas("c1", "generazione particelle", 200, 10, 1200, 800);
  c1->Divide(3, 2);
  c1->cd(1);
  hParticleType->Draw();
  c1->cd(2);
  hDistributionTheta->Draw();
  c1->cd(3);
  hDistributionPhi->Draw();
  c1->cd(4);
  hImpulse->Draw();
  c1->cd(5);
  hImpulseT->Draw();
  c1->cd(6);
  hEnergy->Draw();

  auto c2 = new TCanvas("c2", "massa invariante con varie combinazioni", 200,
                        10, 1200, 800);
  c2->Divide(3, 2);
  c2->cd(1);
  hInvMass->Draw();
  c2->cd(2);
  hInvMassDisc->Draw("h");
  c2->cd(3);
  hInvMassConc->Draw("h");
  c2->cd(4);
  hInvMassComb1->Draw();
  c2->cd(5);
  hInvMassComb2->Draw();
  c2->cd(6);
  hInvMassK->Draw();

  TFile *file = new TFile("output_file.root", "RECREATE");

  hParticleType->Write();
  hDistributionTheta->Write();
  hDistributionPhi->Write();
  hImpulse->Write();
  hImpulseT->Write();
  hEnergy->Write();
  hInvMass->Write();
  hInvMassDisc->Write();
  hInvMassConc->Write();
  hInvMassComb1->Write();
  hInvMassComb2->Write();
  hInvMassK->Write();

  file->Close();
  delete file;
}