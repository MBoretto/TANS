#include <TH1F.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TFile.h>
#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "Direction.h"
#include "TH2D.h"
#include <TStopwatch.h>
#include "TRandom3.h"
#include "TRandom.h"


/*	USAGE: Trasporta(s,Rhum)
		1) s =1 scattering attivo
		   s =0 scattering non attivo
		
		2) Rhum = 0 no rumore
		   Rhum = 1 molteplicita' gaussiana
		   Rhum = n molteplicita' N
		3) Noise_Medio = 20
		se non si inserisce il terzo parametro è settato di default a 20
		altrimenti varia in base al parametro che si inserisce; 100 in caso di analisi con il noise da executeExamNoise
*/

void Trasporta(Int_t s,Int_t Rhum, TRandom *GeneratoreEsterno=0,Int_t Noise_Medio = 20) {

	TStopwatch tempo;

	tempo.Start(kTRUE);

	cout<<endl<<"Sto trasportando attraverso i rivelatori: attendere... "<<endl;
	TRandom *smear;
	if(GeneratoreEsterno == 0){
		smear = new TRandom3();
		cout<<"Generatore Interno";
	}else{
		smear = GeneratoreEsterno;
		cout<<"Generatore Esterno";
	}
	
	cout<<" FirstRNDM: "<<smear->Rndm()<<endl;
	
	//////////////////////////////////////////////////////
	//Creo un nuovo file e
	//Definisco Struct per salvare i nuovi dati x y z 
	//////////////////////////////////////////////////////

	//Definisco il nuovo albero per salvare i punti di hit	
	TFile sfile("trasporto_tree.root","RECREATE");
  
	TTree *trasporto = new TTree("Ttrasporto","TTree con 3 branches");
	 //Punti sul layer
	TTree *Rel_Lay1 = new TTree("Layer1","TTree con 1 branch");
	TTree *Rel_Lay2 = new TTree("Layer2","TTree con 1 branch");
	 //rumore
	TTree *Noise = new TTree("Rumore","TTree con 1 branch");

	typedef struct {
		Double_t X,Y,Z;
		Int_t Flag;		
	} HIT; 
	static HIT beam;  	
	static HIT lay1;  
	static HIT lay2;

	typedef struct {
		Int_t event;
		Int_t tipo;
		Int_t Noiselay1;
		Int_t Noiselay2;
	} infoRumore;
	static infoRumore InfoR;


	//Dichiaro i rami dei tree
	trasporto->Branch("BeamPipe",&beam.X,"X/D:Y:Z:Flag/I");  
	trasporto->Branch("Layer1",&lay1.X,"X/D:Y:Z:Flag/I"); 
	trasporto->Branch("Layer2",&lay2.X,"X/D:Y:Z:Flag/I");  

	Rel_Lay1->Branch("RealLayer1",&lay1.X,"X/D:Y:Z:Flag/I"); 
	Rel_Lay2->Branch("RealLayer2",&lay2.X,"X/D:Y:Z:Flag/I"); 
	
	Noise->Branch("Rumore",&InfoR,"event/I:tipo:Noiselay1:Noiselay2"); 
	Double_t temp_phi = 0;
	Int_t Nnoise=0;	

  ////////////////////////////////
  //Acquisizione Vertici
  ///////////////////////////////
  TClonesArray *dir = new TClonesArray("Direction",100);	
  typedef struct {
    Double_t X,Y,Z;
    Int_t N;
  }SINGLE_EVENT;
  static SINGLE_EVENT event;    //struct con molteplicita' e vertice di un singolo evento
	
  TFile hfile("event_tree.root");


  TTree *Born = (TTree*)hfile.Get("T");       
  TBranch *b1=Born->GetBranch("Event"); 
  TBranch *b2=Born->GetBranch("Direzioni");  //acquisisco i due branches


  b1->SetAddress(&event.X); //passo l'indirizzo del primo oggetto della struct e assegno tutto a b1
  b2->SetAddress(&dir); // lo stesso per il vettore 




  /////////////////////////
  //Geometria del rivelatore
  /////////////////////////
  Double_t R1=3;	//raggio 3 cm beam pipe
  Double_t R2=4;	//raggio 4 cm primo layer
  Double_t R3=7;	//raggio 7 cm secondo layer

  Double_t limit = 8.; //lunghezza layer su z-> z in [-8,8]

  //Variabili Varie
	Double_t Xo=0.;Double_t Yo=0.;Double_t Zo=0.;
  	Double_t X1=0.;Double_t Y1=0.;Double_t Z1=0.;
	Double_t X2=0.;Double_t Y2=0.;Double_t Z2=0.;

	Int_t N=0; //molteplicita'

	Int_t yes = 0;
	Int_t no = 0;	       
	
	for(Int_t e=0; e < Born->GetEntries(); e++){
	
		Born->GetEvent(e);
		Xo=event.X;
		Yo=event.Y;
		Zo=event.Z;		
		N=event.N;	    
		
		for(Int_t i=0; i<N; i++){
			
			//Cast dell'elemenento i di TClones a Direction
			Direction *angolacci=(Direction*)dir->At(i);
			angolacci->SetRNDGenerator(smear);//uso lo stesso generatore anche nella classe
			
			//primo hit beam pipe
			angolacci->GeneraHit(Xo,Yo,Zo,R1);//genero il punto di impatto sul beam pipe		
			
			beam.X=angolacci->GetNewX(); //recupero le coordinate del punto d'impatto sul BP
			beam.Y=angolacci->GetNewY();					
			beam.Z=angolacci->GetNewZ();
			beam.Flag=1;
		

			///////////////////////////////////////////////////
			/////////////scattering sul beam pipe//////////////
			///////////////////////////////////////////////////
			if(s==1){
				//dipende dal tipo di materiale
				angolacci->Scattering(0.08,35.28);
			}
					
				
			//secondo hit layer 1			
			angolacci->GeneraHit(beam.X,beam.Y,beam.Z,R2);

			X1 = angolacci->GetNewX();	
			Y1 = angolacci->GetNewY();
			Z1 = angolacci->GetNewZ();
			
			lay1.X=X1;
			lay1.Y=Y1;							
			lay1.Z=Z1;	
			
			//verifico che la particella colpisca il layer
			if(TMath::Abs(Z1) < limit){

				lay1.Flag = e;
				Rel_Lay1->Fill();					       
				
				///////////////////////////////////////////////
				/////////////scattering sul layer//////////////	
				///////////////////////////////////////////////
				if(s==1){
					angolacci->Scattering(0.02,9.37);
				}				
							
				yes++;
				
			}else no++;	      
				

		      //terzo hit layer 2			
		      angolacci->GeneraHit(X1,Y1,Z1,R3);

		      X2 = angolacci->GetNewX();	
		      Y2 = angolacci->GetNewY();
		      Z2 = angolacci->GetNewZ();
		      lay2.X=X2;
		      lay2.Y=Y2;							
		      lay2.Z=Z2;
			

		      //verifico che la particella colpisca il layer
		      if(TMath::Abs(Z2) < limit){

			lay2.Flag = e;	
			Rel_Lay2->Fill();			
	
			yes++;

		      }else{
			no++;	
		      }

			angolacci->RemoveGenerator();
			trasporto->Fill(); //riempie tutto con quello che ho definito sopra

		      // Debug
		      /*printf("Evento %d : part %d \n",e,i+1);
			printf("x beam= %f ; y beam= %f; z beam= %f \n",beam.X,beam.Y,beam.Z);

			if(lay1.Flag){
			printf("x lay1= %f ; y lay1= %f; z lay1= %f \n",lay1.X,lay1.Y,lay1.Z);
			}else{
			printf("Non urta sul layer 1 \n");
			}


			if(lay2.Flag){
			printf("x lay2= %f ; y lay2= %f; z lay2= %f \n",lay2.X,lay2.Y,lay2.Z);
			}else{
			printf("Non urta sul layer 2 \n");
			}*/
		}

		////////////////////////////////////////////////////////////////////////////
		//////////////////////////RUMORE////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////
		InfoR.event = e;
		InfoR.tipo = Rhum;
		if(Rhum != 0){
		//genero rumore lay 1
			if(Rhum == 1){
				//Nnoise= TMath::Abs(smear->Gaus(20,5));
				//Nnoise = 1+(Int_t)(20*smear->Rndm());
				Nnoise = 1+(Int_t)(Noise_Medio*smear->Rndm());
			}else{
				Nnoise= Rhum;
			}
		
		
			InfoR.Noiselay1 = Nnoise;
		
			for(Int_t y =0; y < Nnoise; y++){
				temp_phi = smear->Uniform(0,2*TMath::Pi());
				lay1.X = R2*TMath::Cos(temp_phi);
				lay1.Y = R2*TMath::Sin(temp_phi);							
				lay1.Z = smear->Uniform(-limit,limit);
	
				lay1.Flag=e;

				Rel_Lay1->Fill();		
		
			}

		      //genero rumore lay 2					
			if(Rhum == 1){
				//Nnoise= TMath::Abs(smear->Gaus(20,5));
				//Nnoise = 1+(Int_t)(20*smear->Rndm());
				Nnoise = 1+(Int_t)(Noise_Medio*smear->Rndm());
			}else{
				Nnoise= Rhum;
			}

			InfoR.Noiselay2 = Nnoise;

			for(Int_t w =0; w < Nnoise; w++){
				temp_phi = smear->Uniform(0,2*TMath::Pi());

				lay2.X = R3*TMath::Cos(temp_phi);
				lay2.Y = R3*TMath::Sin(temp_phi);							
				lay2.Z = smear->Uniform(-limit,limit);
	
				lay2.Flag=e;

				Rel_Lay2->Fill();	
			}
		}else{
			InfoR.Noiselay1 = 0;
			InfoR.Noiselay2 = 0;
		}
	
	//fill per il rumore
	Noise->Fill();
	}


	sfile.Write(); 
 
	sfile.Close();
	//ho il file con tutti gli eventi

  

 	tempo.Stop();
  
	cout<<endl<<endl<<endl<<"//////////////////////////////////////"<<endl<<endl;
	cout<<"Completato!"<<endl<<endl;
	cout<<"Il trasporto è durato "<<endl;
	tempo.Print();
	cout<<endl<<endl;
	cout<<"PARAMETRI TRASPORTO: "<<endl;
	cout<<"\t"<<"Scattering:     "<<s;
	if(s==1)cout<<"  Scattering attivo"<<endl;
	if(s==0)cout<<"  Scattering non attivo"<<endl;
	cout<<"\t"<<"Rumore:         ";
	if(Rhum==1)cout<<"  Rumore gaussiano  "<<endl;
	if(Rhum==0)cout<<"  Nessun rumore"<<endl;
	if((Rhum!=0) & (Rhum!=1))cout<<"  Rumore con molteplicita' fissa "<<Rhum<<endl;
  	cout<<endl<<"//////////////////////////////////////"<<endl;

}
