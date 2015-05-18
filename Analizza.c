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

void Analizza(Int_t N, Int_t Z, Int_t Noise) {
  
  
  TStopwatch tempo;
  
  tempo.Start(kTRUE);
  
  
  TFile *analize = new TFile("ricostruzione_histo_Class.root"); //apro file del salvataggio  
  TNtuple *Ntupla=(TNtuple*)analize->Get("NtupleVertex"); //Apertura entupla
  // TNtuple *Ntupla = new TNtuple ("NtpleVertex","Vertici","molt:moltRhum:Zvero:Zric:Evento:yn");
  
  
  cout<<endl<<"Sto analizzando la ricostruzione degli eventi ... "<<endl;
  
  /////////////////////////////////////Efficienza(Molteplicità) & RMS(Molteplicità)////////////////////////////////////////////////
  cout<<endl<<"Sto analizzando in funzione della molteplicità ... "<<endl;
  
  Double_t Moltep[10];
  Double_t Eff[10];
  Double_t RMS[10];
  Double_t deltaM[10];
  Double_t deltaE[10];
  Double_t errRMS[10];
  Double_t num,den;
  TGraphErrors *Eff_Molt;
  TGraphErrors *RMS_Molt;
 
  TH1F *differenze =new TH1F("differenze","Residui",500,-0.1,0.1);
  //Molteplicità da kinem root
  if(N==1){
	  for(Int_t w=0;w<10;w++){
	    Moltep[w]=4+5*w;
	    num = Ntupla->GetEntries(Form("molt>=5*(%d)+4 && molt<9+5*(%d) && Zvero>-8 && Zvero<8 && yn == 1",w,w));
	    den = Ntupla->GetEntries(Form("molt>=5*(%d)+4 && molt<9+5*(%d) && Zvero>-8 && Zvero<8",w,w));
	    Ntupla->Draw("Zvero-Zric>>differenze",Form("molt >= 5+4*(%d) && molt < 9+4*(%d)",w,w));
	    RMS[w]=differenze->GetRMS();
	    if(den != 0){    
	      Eff[w]=num/den;
	    }else{
	      cout<<"No data found "<<endl;
	      Eff[w]=0.;
	    }
	    deltaM[w]=1.;
	    deltaE[w]=(TMath::Sqrt(num*(1-Eff[w]))/den);
	    errRMS[w]=differenze->GetRMSError();
	    differenze->Reset();
	  }
	  Eff_Molt=new TGraphErrors(10,Moltep,Eff,deltaM,deltaE);
	  Eff_Molt->SetTitle("Efficienza(Molteplicita')");
	  Eff_Molt->GetXaxis()->SetTitle("Molt (#)");
	  Eff_Molt->GetYaxis()->SetTitle("Eff");
	  
	  RMS_Molt=new TGraphErrors(10,Moltep,RMS,deltaM,errRMS);
	  RMS_Molt->SetTitle("RMS(Molteplicita')");
	  RMS_Molt->GetXaxis()->SetTitle("Molt (#)");
	  RMS_Molt->GetYaxis()->SetTitle("RMS (cm)");
	  differenze->Reset();
  }
  ////////////////////////////////////////////Efficienza(Zeta) & RMS(Zeta)/////////////////////////////////////////////////////
  cout<<endl<<"Sto analizzando in funzione della coordinata zeta ... "<<endl;
  Double_t Zeta[15];
  Double_t Eff2[15];
  Double_t RMS2[15];
  Double_t deltaZeta[15];
  Double_t deltaEff2[15];
  Double_t deltaRMS[15];
  Double_t num2[15],den2[15];
  TGraphErrors *Eff_Zeta;
  TGraphErrors *RMS_Zeta;
  
  if(Z==-1){
  for(Int_t x=0;x<15;x++){
    Zeta[x]=2*x-15;
    num2[x]=Ntupla->GetEntries(Form("Zvero>=2*(%d)-15 && Zvero<2*(%d)-13 && molt>10 && yn == 1",x,x));
    den2[x]=Ntupla->GetEntries(Form("Zvero>=2*(%d)-15 && Zvero<2*(%d)-13 && molt>10",x,x)); 
    Ntupla->Draw("Zvero-Zric>>differenze",Form("molt>15 && Zvero >= 2*(%d) -15 && Zvero < 2*(%d)-13 ",x,x));
    RMS2[x]=differenze->GetRMS();
    if(den2[x] != 0){
      Eff2[x]=num2[x]/den2[x];
    }else{
      cout<<"No data found "<<endl;
      Eff2[x]=0.;
    }
    deltaZeta[x]=0.1;
    deltaEff2[x]=(TMath::Sqrt(num2[x]*(1-Eff2[x]))/den2[x]);
    deltaRMS[x]=differenze->GetRMSError();
    differenze->Reset();
  }
  Eff_Zeta=new TGraphErrors(15,Zeta,Eff2,deltaZeta,deltaEff2);
  Eff_Zeta->SetTitle("Efficienza(Z)");
  Eff_Zeta->GetXaxis()->SetTitle("Z (cm)");
  Eff_Zeta->GetYaxis()->SetTitle("Eff");
  RMS_Zeta=new TGraphErrors(15,Zeta,RMS2,deltaZeta,deltaRMS);
  RMS_Zeta->SetTitle("RMS(Zeta)");
  RMS_Zeta->GetXaxis()->SetTitle("Zeta (cm)");
  RMS_Zeta->GetYaxis()->SetTitle("RMS (cm)");
  differenze->Reset();
  }
  
  ///////////////////////////////////////////////Efficienza(Rumore) & RMS(Rumore)///////////////////////////////////////////////
  cout<<endl<<"Sto analizzando in funzione della molteplicità di rumore ... "<<endl;
  Double_t Rhum[10];
  Double_t Eff3[10];
  Double_t RMS3[10];
  Double_t deltaRhum[10];
  Double_t deltaEff3[10];
  Double_t deltaRMS3[10];
  Double_t num3[10],den3[10];
  TGraphErrors *Eff_Rhum;
  TGraphErrors *RMS_Noise;
  
  if(Noise==1){
  for(Int_t k = 0; k<10; k++){
    Rhum[k] = 1+2*k;
    num3[k] = Ntupla->GetEntries(Form("Zvero>-5 && Zvero<5 && yn==1 && moltRhum>=1 + 2*(%d) && moltRhum<3+2*(%d)",k,k));
    den3[k] = Ntupla->GetEntries(Form("Zvero>-5 && Zvero<5 && moltRhum>=1 + 2*(%d) && moltRhum<3+2*(%d)",k,k));
    Ntupla->Draw("Zvero-Zric>>differenze",Form("molt>10 && moltRhum >= 1+2*(%d) && moltRhum < 3+2*(%d) ",k,k));
    RMS3[k]=differenze->GetRMS();
    if(den3[k] != 0){
      Eff3[k] = num3[k]/den3[k];
    }else{
      cout<<"No data found "<<endl;
      Eff3[k] = 0.;
    }
    deltaRhum[k] = 0.;
    deltaEff3[k] = (TMath::Sqrt(num3[k]*(1-Eff3[k]))/den3[k]);
    deltaRMS3[k]=differenze->GetRMSError();
    differenze->Reset();
  }
  
  Eff_Rhum=new TGraphErrors(10,Rhum,Eff3,deltaRhum,deltaEff3);
  Eff_Rhum->SetTitle("Efficienza(Noise)");
  Eff_Rhum->GetXaxis()->SetTitle("Molt_Noise (#)");
  Eff_Rhum->GetYaxis()->SetTitle("Eff");
  RMS_Noise=new TGraphErrors(10,Rhum,RMS3,deltaRhum,deltaRMS3);
  RMS_Noise->SetTitle("RMS(Noise)");
  RMS_Noise->GetXaxis()->SetTitle("Noise (#)");
  RMS_Noise->GetYaxis()->SetTitle("RMS (cm)");
  differenze->Reset();
  }
  delete differenze;
  
  ///////////////////////////////Efficienza(Molteplicità) & Efficienza(Zeta) a diversi rumori//////////////////////////////////
  
  Double_t Xmol[18];
  Double_t Yeff[18];
  Double_t Xzeta[15];
  Double_t Yeffzeta[15];
  Double_t deltaX[18];
  Double_t deltaY[18];
  Double_t deltaXzeta[15];
  Double_t deltaYzeta[15];
  Double_t num4,den4,num5,den5;
  TGraphErrors *Eff_Molt_Rhum[4];
  TGraphErrors *Eff_Zeta_Rhum[4];
  TMultiGraph *EffClassiNoise = new TMultiGraph();
  TMultiGraph *EffZetaClassiNoise = new TMultiGraph();
  if(N==1 && Z==-1 && Noise==1){
  for(Int_t r = 0; r<4; r++){
    for(Int_t a=0;a<18;a++){
      Xmol[a]=2+3*a;
      num4=Ntupla->GetEntries(Form("molt>=3*(%d)+2 && molt<5+3*(%d) && Zvero>-8 && Zvero<8 && yn == 1 && moltRhum>=1+5*(%d) && moltRhum<6+5*(%d)",a,a,r,r));
      den4=Ntupla->GetEntries(Form("molt>=3*(%d)+2 && molt<5+3*(%d) && Zvero>-8 && Zvero<8 && moltRhum>=1+5*(%d) && moltRhum<6+5*(%d)",a,a,r,r));
      if(den4 != 0){    
	Yeff[a]=num4/den4;
      }else{
	cout<<"No data found "<<endl;
	Yeff[a]=0.;
      }
      deltaX[a]=1.;
      deltaY[a]=(TMath::Sqrt(den4*Yeff[a]*(1-Yeff[a]))/den4);
    }
    Eff_Molt_Rhum[r] = new TGraphErrors(18,Xmol,Yeff,deltaX,deltaY);
    Eff_Molt_Rhum[r]->GetXaxis()->SetTitle("Molt (#)");
    Eff_Molt_Rhum[r]->GetYaxis()->SetTitle("Eff");
    Eff_Molt_Rhum[r]->SetMarkerColor(r+2);
    Eff_Molt_Rhum[r]->SetMarkerStyle(33);
    Eff_Molt_Rhum[r]->SetTitle(Form("Efficienza(Molteplicita') a noise tra %d e %d punti per layer", 1+5*r,6+5*r));
    EffClassiNoise->Add(Eff_Molt_Rhum[r],"p");
    
    for(Int_t b=0;b<15;b++){
      Xzeta[b]=2*b-14;
      num5=Ntupla->GetEntries(Form("molt>10 && Zvero>=2*(%d)-15 && Zvero<2*(%d)-13 && yn == 1 && moltRhum>=1+5*(%d) && moltRhum<6+5*(%d)",b,b,r,r));
      den5=Ntupla->GetEntries(Form("molt>10 && Zvero>=2*(%d)-15 && Zvero<2*(%d)-13 && moltRhum>=1+5*(%d) && moltRhum<6+5*(%d)",b,b,r,r));
      if(den5 != 0){    
	Yeffzeta[b]=num5/den5;
      }else{
	cout<<"No data found "<<endl;
	Yeffzeta[b]=0.;
      }
      deltaXzeta[b]=0.1;
      deltaYzeta[b]=(TMath::Sqrt(den5*Yeffzeta[b]*(1-Yeffzeta[b]))/den5);
    }
    Eff_Zeta_Rhum[r] = new TGraphErrors(15,Xzeta,Yeffzeta,deltaXzeta,deltaYzeta);
    Eff_Zeta_Rhum[r]->GetXaxis()->SetTitle("Zeta (cm)");
    Eff_Zeta_Rhum[r]->GetYaxis()->SetTitle("Eff");
    Eff_Zeta_Rhum[r]->SetMarkerColor(r+2);
    Eff_Zeta_Rhum[r]->SetMarkerStyle(33);
    Eff_Zeta_Rhum[r]->SetTitle(Form("Efficienza(Zeta) a noise tra %d e %d punti per layer", 1+5*r,6+5*r));
    EffZetaClassiNoise->Add(Eff_Zeta_Rhum[r],"p");
    
  }
  
  }
  
  TCanvas *finale = new TCanvas("analisi","analisi");
  finale->Divide(3,2);
  
  finale->cd(1);
  if(N==1)Eff_Molt->Draw("AP");
  finale->cd(2);
  if(Z==-1)Eff_Zeta->Draw("AP");
  finale->cd(3);
  if(Noise==1)Eff_Rhum->Draw("AP");
  finale->cd(4);
  if(N==1)RMS_Molt->Draw("AP");
  finale->cd(5);
  if(Z==-1)RMS_Zeta->Draw("AP");
  finale->cd(6);
  if(Noise==1)RMS_Noise->Draw("AP");
  

  if(N==1 && Z==-1 && Noise==1){
  TCanvas *analisiNoise = new TCanvas("Analisi Noise","Analisi Noise");
  analisiNoise->Divide(2);
  analisiNoise->cd(1);
  
  //EffClassiNoise->GetXaxis()->SetTitle("Molt(#)");
  //EffClassiNoise->GetYaxis()->SetTitle("Eff");
  EffClassiNoise->SetTitle("Efficienza(Molt) a intervalli di rumore");
  EffClassiNoise->Draw("AP");
  
  analisiNoise->cd(2);
  EffZetaClassiNoise->Draw("AP");
  EffZetaClassiNoise->SetTitle("Efficienza(Zeta) a intervalli di rumore");
  //EffZetaClassiNoise->GetXaxis()->SetTitle("Zeta(cm)");
  //EffZetaClassiNoise->GetYaxis()->SetTitle("Eff");
  }
  TFile sfile("analisi.root","RECREATE");
  
  if(N==1)Eff_Molt->Write("Eff(Molt)");
  if(Z==-1)Eff_Zeta->Write("Eff(Zeta)");
  if(Noise==1)Eff_Rhum->Write("Eff(Noise)");
  if(N==1)RMS_Molt->Write("RMS(Molt)");
  if(Z==-1)RMS_Zeta->Write("RMS(Zeta)");
  if(Noise==1)RMS_Noise->Write("RMS(Noise)");
  
  if(N==1 && Z==-1 && Noise==1){
  Eff_Molt_Rhum[0]-> Write("Eff(Molt) a noise da 1 a 5");
  Eff_Molt_Rhum[1]-> Write("Eff(Molt) a noise da 6 a 10");
  Eff_Molt_Rhum[2]-> Write("Eff(Molt) a noise da 11 a 15");
  Eff_Molt_Rhum[3]-> Write("Eff(Molt) a noise da 16 a 20");
  Eff_Zeta_Rhum[0]-> Write("Eff(Zeta) a noise da 1 a 5");
  Eff_Zeta_Rhum[1]-> Write("Eff(Zeta) a noise da 6 a 10");
  Eff_Zeta_Rhum[2]-> Write("Eff(Zeta) a noise da 11 a 15");
  Eff_Zeta_Rhum[3]-> Write("Eff(Zeta) a noise da 16 a 20");
  EffClassiNoise->Write("Eff(Molt) a intervalli di noise");
  EffZetaClassiNoise->Write("Eff(Zeta) a intervalli di noise");
  }
  finale->Write();
  
  sfile.Close();
  
  
  tempo.Stop();
  
  cout<<"L'analisi è durata "<<endl;
  tempo.Print();
  
  
 

}
