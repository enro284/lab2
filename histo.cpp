#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"

using namespace std;

#define NBINS 100

void histo() {
  TFile *file = TFile::Open("output_file.root", "READ");

  auto hParticleType = file->Get<TH1I>("hParticleType");
  auto hDistributionTheta = file->Get<TH1D>("hDistributionTheta");
  auto hDistributionPhi = file->Get<TH1D>("hDistributionPhi");
  auto hImpulse = file->Get<TH1D>("hImpulse");
  auto hImpulseT = file->Get<TH1D>("hImpulseT");
  auto hEnergy = file->Get<TH1D>("hEnergy");
  auto hInvMass = file->Get<TH1D>("hInvMass");
  auto hInvMassDisc = file->Get<TH1D>("hInvMassDisc");
  auto hInvMassConc = file->Get<TH1D>("hInvMassConc");
  auto hInvMassComb1 = file->Get<TH1D>("hInvMassComb1");
  auto hInvMassComb2 = file->Get<TH1D>("hInvMassComb2");
  auto hInvMassK = file->Get<TH1D>("hInvMassK");

  auto c1 = new TCanvas("c1", "generazione particelle", 200, 10, 1200, 800);
  c1->Divide(2, 2);
  c1->cd(1);
  hParticleType->LabelsOption(">");
  hParticleType->Draw();
  c1->cd(2);
  hDistributionTheta->Draw();
  c1->cd(3);
  hDistributionPhi->Draw();
  c1->cd(4);
  hImpulse->Draw();

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

  cout << " Contenuto dei Bin hPartycleType : "
       << hParticleType->GetBinContent(0) << endl;
  cout << " Errori bin hPartycleType : " << hParticleType->GetBinError(0)
       << endl;

  TF1 *fTheta = new TF1("fTheta", "[0]", 0., TMath::Pi());
  hDistributionTheta->Fit(fTheta);
  cout << " Valore parametro fit fTheta: " << fTheta->GetParameter(0) << '\n'
       << endl;
  cout << " Chiquadro ridotto fit fTheta: "
       << fTheta->GetChisquare() / fTheta->GetNDF() << endl;
  cout << " Probabilità del fit fTheta: " << fTheta->GetProb() << endl;

  TF1 *fPhi = new TF1("fPhi", "[0]", 0., 2 * TMath::Pi());
  hDistributionPhi->Fit(fPhi);
  cout << " Valore parametro fit fPhi: " << fPhi->GetParameter(0) << '\n'
       << endl;
  cout << " Chiquadro ridotto fit fPhi: "
       << fPhi->GetChisquare() / fPhi->GetNDF() << endl;
  cout << " Probabilità del fit fPhi: " << fPhi->GetProb() << endl;

  TF1 *fImpulse = new TF1("fImpulse", "[0]*exp([1]*x)", 0., 8.);
  fImpulse->SetParameter(0, 80000);
  fImpulse->SetParameter(1, 1.5);

  hImpulse->Fit(fImpulse);

  cout << " Valore parametro 1 fit Impulso: " << fImpulse->GetParameter(0)
       << endl;
  cout << " Valore parametro 2 fit Impulso: " << fImpulse->GetParameter(1)
       << endl;
  cout << " Chiquadro ridotto fit Impulso: "
       << fImpulse->GetChisquare() / fImpulse->GetNDF() << endl;
  cout << " Probabilità del fit Impulso: " << fImpulse->GetProb() << endl;
  cout << "Media dell'istogramma : "<< hImpulse->GetMean() << "+/-" << hImpulse->GetMeanError()<< endl;

  auto hSottr12 = (TH1D *)hInvMassDisc->Clone("hSottr12");
  hSottr12->Add(hInvMassConc, -1.);
  // cosmetica
  hSottr12->SetTitle("Massa invariante, sottrazione 1");
  hSottr12->GetXaxis()->SetTitle("Massa Invariante [GeV/c^2]");
  hSottr12->GetYaxis()->SetTitle("Occorrenze");
  /*
    hSottr12->SetAxisRange(0., 5., "X");
    hSottr12->SetAxisRange(-20000., 200000., "Y");
  */

  // grafico sottrazione 34
  auto hSottr34 = (TH1D *)hInvMassComb1->Clone("hSottr34");
  hSottr34->Sumw2();
  hSottr34->Add(hInvMassComb2, -1.);
  // cosmetica
  hSottr34->SetTitle("Massa invariante, sottrazione 2");
  hSottr34->GetXaxis()->SetTitle("Massa Invariante [GeV/c^2]");
  hSottr34->GetYaxis()->SetTitle("Occorrenze");
  /* hSottr34->SetAxisRange(xmin,xmax,*axis="X");
  hSottr34->SetAxisRange(ymin,ymax,*axis="Y"); */

  auto c3 = new TCanvas("c3", "Istogramma sottrazione 1-2 e 3-4 e K*", 200, 10,
                        1200, 800);
  c3->Divide(2, 2);
  c3->cd(1);
  hSottr12->Draw("h");
  c3->cd(2);
  hSottr34->Draw("h");
  c3->cd(3);
  hInvMassK->Draw();

  auto fGaus12 = new TF1("fGaus12", "[0]*exp(-(x-[1])^2/(2*[2]^2))", 0., 8.);
  fGaus12->SetParameter(0, 1E5);
  fGaus12->SetParameter(1, 0.89);
  fGaus12->SetParameter(2, 0.5);

  hSottr12->Fit(fGaus12);

  cout << " Valore ampiezza fit Gaussiana 12: " << fGaus12->GetParameter(0) << "+/-" << fGaus12->GetParError(0) 
       << '\n'
       << endl;
  cout << " Valore media fit Gaussiana 12: " << fGaus12->GetParameter(1) << "+/-" << fGaus12->GetParError(1) << '\n'
       << endl;
  cout << " Valore deviazione standard fit Gaussiana 12: "
       << fGaus12->GetParameter(2) << "+/-" << fGaus12->GetParError(2) << '\n'
       << endl;
  cout << " Chiquadro ridotto fit Gaussiana 12 : "
       << fGaus12->GetChisquare() / fGaus12->GetNDF() << endl;
  cout << " Probabilità del fit Gaussiana 12 : " << fGaus12->GetProb() << endl;

  auto fGaus34 = new TF1("fGaus34", "[0]*exp(-(x-[1])^2/(2*[2]^2))", 0., 8.);
  fGaus34->SetParameter(0, 1E5);
  fGaus34->SetParameter(1, 0.8);
  fGaus34->SetParameter(2, 0.5);
  hSottr34->Fit(fGaus34);

  cout << " Valore ampiezza fit Gaussiana 34 : " << fGaus34->GetParameter(0) << "+/-" << fGaus34->GetParError(0) 
       << '\n'
       << endl;
  cout << " Valore media fit Gaussiana 34 : " << fGaus34->GetParameter(1) << "+/-" << fGaus34->GetParError(1) 
       << '\n'
       << endl;
  cout << " Valore deviazione standard fit Gaussiana 34 : "
       << fGaus34->GetParameter(2) << "+/-" << fGaus34->GetParError(2) 
       << '\n'
       << endl;
  cout << " Chiquadro ridotto fit Gaussiana 34 : "
       << fGaus34->GetChisquare() / fGaus34->GetNDF() << endl;
  cout << " Probabilità del fit Gaussiana 34: " << fGaus34->GetProb() << endl;

  auto fK = new TF1("fK", "[0]*exp(-(x-[1])^2/(2*[2]^2))", 0., 8.);
  fK->SetParameter(0, 6000);
  fK->SetParameter(1, 0.8);
  fK->SetParameter(2, 0.05);
  hInvMassK->Fit(fK);
  cout << " Valore ampiezza fit K*: " << fK->GetParameter(0) << "+/-" << fK->GetParError(0) << '\n' << endl;
  cout << " Valore media fit K* : " << fK->GetParameter(1) << "+/-" << fK->GetParError(1) << '\n' << endl;
  cout << " Valore deviazione standard fit K* : " << fK->GetParameter(2) << "+/-" << fK->GetParError(2)  << '\n'
       << endl;
  cout << " Chiquadro ridotto fit K*: " << fK->GetChisquare() / fK->GetNDF()
       << endl;
  cout << " Probabilità del fit K*: " << fK->GetProb() << endl;

  cout << "pigreco + :" << hParticleType->GetBinContent(1) << "+/-" << hParticleType->GetBinError(1)<<endl;
  cout << "pigreco - :" << hParticleType->GetBinContent(2) << "+/-" << hParticleType->GetBinError(2)<<endl;
  cout << "k+ :" << hParticleType->GetBinContent(3) << "+/-" << hParticleType->GetBinError(3)<<endl;
  cout << "k- :" << hParticleType->GetBinContent(4) << "+/-" << hParticleType->GetBinError(4)<<endl;
  cout << "p+ :" << hParticleType->GetBinContent(5) << "+/-" << hParticleType->GetBinError(5)<<endl;
  cout << "p- :" << hParticleType->GetBinContent(6) << "+/-" << hParticleType->GetBinError(6)<<endl;
  cout << "k* : " << hParticleType->GetBinContent(7) << "+/-" << hParticleType->GetBinError(7)<<endl;
}