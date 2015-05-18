#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>
#include "Fato.h"
#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "Direction.h"
#include "TStopwatch.h"

/*	USAGE: Genera(n,esp,k,w,Gen)

		1) n = numero eventi * 10^(esp)
		
                2) esp = esponente numero eventi
	
		3) k = numero -> molteplicita' fissa numero
		   k = 1 usa la distribuzione di molteplicita' da kinem.root
		
		4) w = posizione fissa di z dei vertici
		   w = -1 usa distribuzione casuale per z dei vertici
		5) puntatore al generatore esterno
*/
 
void Genera(Int_t n, Int_t esp, Int_t k,Int_t w, TRandom *GeneratoreEsterno=0){

  TStopwatch tempo;

  tempo.Start(kTRUE);

  Int_t Nevents = n*TMath::Power(10,esp);
  
  cout<<endl<<"Sto generando "<<Nevents<<" eventi: attendere... "<<endl;
	
	
	Fato *simulazione = new Fato(k,GeneratoreEsterno);
	if(GeneratoreEsterno == 0){
		cout<<"Generatore Interno";
	}else{
		cout<<"Generatore Esterno";
	}
	cout<<" FirstRNDM: "<<simulazione->Rndm()<<endl;

	TFile MMFB("event_tree.root","recreate");

	TTree *Born = new TTree ("T","Albero di un evento");  //struct con molteplicita' e vertice e un array di direzioni


	TClonesArray *direction = new TClonesArray("Direction",100);
  	TClonesArray &dir = *direction;	

	
	typedef struct {
		Double_t X,Y,Z;
		Int_t N;
		}SINGLE_EVENT;
	static SINGLE_EVENT event;    //struct con molteplicita' e vertice di un singolo evento
	
	
	
	Born->Branch("Event",&event.X,"X/D:Y:Z:N/I"); //dichiaro il ramo dell'albero che contiene la struct	

	Born->Branch("Direzioni",&direction); //dichiaro il ramo che contiene le direzioni
	
	
	for(Int_t i = 0;i < Nevents; i++){
		
		Int_t molt=0;
		simulazione->NParticelle(molt);
		
		event.N=molt;
	

		if ((w!=-1) && (w!=0)){
			event.Z=w;
			event.X=simulazione->VertX_Y(); 
			event.Y=simulazione->VertX_Y();
			
		}

		if(w==0){ //utile per i test sulla ricostruzione
		  event.Z=0;
		  event.X=0; 
		  event.Y=0;
		  
		}
		if(w==-1){
			event.Z=simulazione->VertZ();
			event.X=simulazione->VertX_Y();
			event.Y=simulazione->VertX_Y();
			
		}
	

	
		Double_t phi, theta;
	
		for(Int_t j = 0;j < event.N; j++){

			simulazione->Direzioni(phi,theta);
	
			new (dir[j]) Direction(phi,theta);
		
		}

		// Debug
		/* printf("Evento %d - molteplicità: %d\n",i,molt);
		   printf("		Vertice x= %f ; y= %f; z= %f \n",event.X,event.Y,event.Z);
		   printf("Entries nel TClonesArray: %d\n",direction->GetEntries());
		   cout<<"      "<<endl;		
		   //for (Int_t j=0; j<molt; j++){
		   //Direction *tst=(Direction*)direction->At(j);
		   //cout<<"Direzione "<<") theta, phi = "<<tst->Gettheta()<<"; "<<tst->Getphi()<<endl;
		   //}*/
		// fine del debug
		

		Born->Fill();	
		direction->Clear();

	}


	MMFB.Write();
	MMFB.Close();

	tempo.Stop();
	
	cout<<endl<<endl<<endl<<"//////////////////////////////////////"<<endl<<endl;

	cout<<"Completato ! "<<endl<<endl;
	cout<<"La generazione degli eventi è durata "<<endl;
	tempo.Print();
	cout<<endl<<endl;
	cout<<"PARAMETRI GENERAZIONE: "<<endl;
	cout<<"Numero eventi:            "<<Nevents<<endl;
	cout<<"Molteplicita':            ";
	if(k==1){cout<<" da distribuzione in kinem.root"<<endl;}else{cout<<k<<endl;}
	cout<<"Posizione Z vertici:      ";
	if(w==-1){cout<<" da distribuzione casuale"<<endl;}else{cout<<w<<endl;}

	cout<<endl<<"//////////////////////////////////////"<<endl;
	

}

