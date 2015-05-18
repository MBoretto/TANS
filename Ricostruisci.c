#include <TH1F.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TFile.h>
#include "TMath.h"
#include "TNtuple.h"
#include <TStopwatch.h>
#include <vector> 
#include "TStyle.h"
#include "MVertex.h"


void Ricostruisci(Int_t Smearing, TRandom *GeneratoreEsterno=0){
	
	
	TNtuple *Ntupla = new TNtuple ("NtupleVertex","Vertici","molt:moltRhum:Zvero:Zric:Evento:yn");
	Float_t entupla[6];

	//Distribuzione delle differenze
	TH1F *Differenze = new TH1F("FrequenzeDifferenze","FrequenzeDifferenze",1000,-1.001,0.999);

	TStopwatch tempo;
	tempo.Start(kTRUE);

	
	
	//invoco la classe per ricostruire
	MVertex *calcola = new MVertex("event_tree.root", "trasporto_tree.root", Smearing, GeneratoreEsterno);
	if(GeneratoreEsterno == 0){		
		cout<<"Generatore Interno";
	}else{
		cout<<"Generatore Esterno";
	}
	GeneratoreEsterno = calcola->GetInternalGenarator();
	cout<<" FirstRNDM: "<< GeneratoreEsterno->Rndm()<<endl;
		
	cout<<endl<<"Sto ricostruendo "<<calcola->GetNevent()<<" eventi: attendere... "<<endl;
	Int_t i = 0;
	for(Int_t ev = 0; ev < calcola->GetNevent(); ev++){
	
		calcola->FindVertex(ev);
		//prendo i dati che la classe ha calcolato e i relativi presenti nei file
		entupla[0] = calcola->GetMolteplicita();
		entupla[1] = calcola->GetNoiselay1() + calcola->GetNoiselay2();		
		entupla[2] = calcola->GetVerticeVero();
		entupla[3] = calcola->GetVertex();
		entupla[4] = ev;
		//flag che dice se ci sono stati dei problemi o se la ricostruzione è andata a buon fine
		entupla[5] = calcola->GetSoFarSoGood();
	
		if(calcola->GetSoFarSoGood()){
		  Differenze->Fill(calcola->GetVertex() - calcola->GetVerticeVero());
		  //cout<<ev<<" M: "<<calcola->GetMolteplicita()<<" V: "<<calcola->GetVertex()<<" invece che: "<<calcola->GetVerticeVero()<<endl;
		  
		}else{	
			
		  //cout<<ev<<" M: "<<calcola->GetMolteplicita()<<" V: "<<calcola->GetVertex()<<" invece che:"<<calcola->GetVerticeVero();
		  i++;
		  //entupla[2] = null;
		  //cout<<ev<<" ricostruzone fallita a causa di:";
		  for(Int_t e = 1; e <= 5; e++){
		    //cout<<" "<<calcola->GetError(e);
		  }
		  //cout<<" Vedi funzione debug per chiarimenti in MVertex.cxx"<<endl;
		  

		}

		

		Ntupla->Fill(entupla);
	}
	cout<<"Esclusi: "<<i<<endl;
	
	//delete calcola;
	tempo.Stop();

	TCanvas *ricostruzione = new TCanvas("Ricostruzione","Ricostruzione");
	
	//Fit Differenze
	ricostruzione->cd(1);
	ricostruzione->SetLogy();
	Differenze->Fit("gaus");
	gStyle->SetOptFit(kTRUE);

	Differenze->SetDirectory(0);
	Differenze->Draw();

	

	//Differenze->Write();
	//Ntupla->Write();
	
	//devo passare i grafici nella classe
	calcola->SaveHisto(Ntupla,Differenze);
	//file->Close();
	
	calcola->CloseFile();

	cout<<endl<<endl<<endl<<"//////////////////////////////////////"<<endl<<endl;
	cout<<"Completato!"<<endl<<endl;
	cout<<"La ricostruzione è durata "<<endl;
	tempo.Print();
	cout<<endl<<endl;
	cout<<"PARAMETRI RICOSTRUZIONE: "<<endl;
	cout<<"\t"<<"Smearing:     ";
	if(Smearing == 1)cout<<"Smearing attivo "<<endl;
	if(Smearing == 0)cout<<"Smearing non attivo "<<endl;
	cout<<endl<<"//////////////////////////////////////"<<endl;
	
	
}
