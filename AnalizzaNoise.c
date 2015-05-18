#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Fato.h"
#include "Direction.h"
#include <vector>
#include "TH1.h"
#include <TCanvas.h>
#include <TH2D.h>
#include "TNtuple.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAttMarker.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TPad.h"

void AnalizzaNoise() {
  
  
  TStopwatch tempo;
  
  tempo.Start(kTRUE);
  
  TFile *analizeNoise = new TFile("ricostruzione_histo_Class.root"); //apro file del salvataggio  
  TNtuple *NtuplaNoise=(TNtuple*)analizeNoise->Get("NtupleVertex"); //Apertura entupla
  // TNtuple *Ntupla = new TNtuple ("NtpleVertex","Vertici","molt:moltRhum:Zvero:Zric:Evento:yn");
  
  
  cout<<endl<<"Sto analizzando la ricostruzione degli eventi con lo studio sul noise ... "<<endl;
  
  ///////////////////////////////////////////////Efficienza(Rumore) & RMS(Rumore)///////////////////////////////////////////////
  
  Double_t Rhum[11];
  Double_t Eff3[11];
  Double_t Noise[11];
  Double_t RMS3[11];
  Double_t deltaRhum[11];
  Double_t deltaEff3[11];
  Double_t deltaRMS3[11];
  Double_t num3[11],den3[11];
  TGraphErrors *Eff_Rhum;
  TGraphErrors *RMS_Noise;
  TH1F *differenze =new TH1F("differenze","Residui",500,-0.1,0.1);
  
  for(Int_t k = 0; k<11; k++){
    Rhum[k] = 1+10*k;
    Noise[k] = 1+10*k;
    num3[k] = NtuplaNoise->GetEntries(Form("Zvero>-5 && Zvero<5 && yn==1 && moltRhum>=1 + 10*(%d) && moltRhum<11+10*(%d)",k,k));    
    den3[k] = NtuplaNoise->GetEntries(Form("Zvero>-5 && Zvero<5 && moltRhum>=1 + 10*(%d) && moltRhum<11+10*(%d)",k,k));
    NtuplaNoise->Draw("Zvero-Zric>>differenze",Form("molt>10 && moltRhum >= 1+10*(%d) && moltRhum < 11+10*(%d) ",k,k));
    RMS3[k]=differenze->GetRMS();
    if(den3[k] != 0){
      Eff3[k] = num3[k]/den3[k];
    }else{
      Eff3[k] = num3[k];
    }
    deltaRhum[k] = 0.;
    deltaEff3[k] = (TMath::Sqrt(num3[k]*(1-Eff3[k]))/den3[k]);
    deltaRMS3[k]=differenze->GetRMSError();
    differenze->Reset();
}
  Eff_Rhum=new TGraphErrors(11,Rhum,Eff3,deltaRhum,deltaEff3);
  Eff_Rhum->SetTitle("Efficienza(Noise)");
  Eff_Rhum->GetXaxis()->SetTitle("Noise (#)");
  Eff_Rhum->GetYaxis()->SetTitle("Eff");
  RMS_Noise=new TGraphErrors(11,Noise,RMS3,deltaRhum,deltaRMS3);
  RMS_Noise->SetTitle("RMS(Noise)");
  RMS_Noise->GetXaxis()->SetTitle("Noise (#)");
  RMS_Noise->GetYaxis()->SetTitle("RMS (cm)");
  delete differenze;
  
  
  TCanvas *analisinoise = new TCanvas("analisiNoise","analisiNoise");
  analisinoise->Divide(2);

  
  analisinoise->cd(1);
  Eff_Rhum->Draw("AP");
  
  analisinoise->cd(2);
  RMS_Noise->Draw("AP");
  
 
  TFile lfile("analisiNoise.root","RECREATE");
  
  
  Eff_Rhum->Write("Eff(Noise)");
  RMS_Noise->Write("RMS(Noise)");
  
  
  analisinoise->Write();
  
  
  lfile.Close();
  
  
  tempo.Stop();
  
  cout<<"L'analisi è durata "<<endl;
  tempo.Print();
    

}
