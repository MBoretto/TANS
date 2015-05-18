#include <Fato.h>
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>


ClassImp(Fato)

//----------COSTRUTTORE DI DEFAULT-------------------------------------------------------------------

Fato::Fato (TRandom * RNDM): TRandom3(),

  fdistribuzione(1)

  
{
//link al generatore
GeneratoreEsterno = RNDM;


}



//----------------COSTRUTTORE STANDARD----------------------------------------------------------------

Fato::Fato (Int_t k,TRandom * RNDM): TRandom3(),
   
  fdistribuzione(k) 
  	
    
{
	//link al generatore
	GeneratoreEsterno = RNDM;

	TFile F("kinem.root");
	rap = (TH1F*) F.Get("heta");
	mol = (TH1F*) F.Get("hmul");

	rap->SetDirectory(0);
	mol->SetDirectory(0);

	F.Close(); 

}  


//----------------DISTRUTTORE STANDARD-----------------------------------------------------------------

Fato:: ~Fato() {

delete rap;
delete mol;

}

//---------------FUNZIONE PER COORDINATE VERTICE----------------------------------------------------------------


Double_t Fato::VertZ() {
	
	Double_t Temp = 0;
		if(GeneratoreEsterno != 0){
			Temp = GeneratoreEsterno->Gaus(0.,5.3) ;
			
		}else{
			Temp = this->Gaus(0.,5.3) ;
		}
	return Temp;
}

Double_t Fato::VertX_Y() {

	Double_t Temp = 0;
	if(GeneratoreEsterno != 0){
		Temp = GeneratoreEsterno->Gaus(0.,0.01);
		
	}else{
		Temp = this->Gaus(0.,0.01);
	}
	return Temp;

}




//---------------FUNZIONE MOLTEPLICITA'----------------------------------------------------------------


void Fato::NParticelle(Int_t &N) {
	
	N = 0 ;


        if (fdistribuzione == 1){
	  N = 2 + mol->GetRandom(); //molteplicità minima = 2, distribuzione da kinem.root
	}else{
	  N=fdistribuzione; //molteplicità fissa 
	}
	
		
	       
}


//---------------FUNZIONE GENERATORE DIREZIONE---------------------------------


void Fato::Direzioni(Double_t &phi,Double_t &theta) {

	//phi = 2*TMath::Pi()*gRandom->Rndm();
	if(GeneratoreEsterno != 0){
		phi = 2*TMath::Pi()*GeneratoreEsterno->Rndm();
	}else{
		phi = 2*TMath::Pi()*this->Rndm();
	}
	Double_t etha;
	do{
	etha = rap->GetRandom();

		}while (etha>2||etha<-2); 

	theta = 2*TMath::ATan(TMath::Exp(-etha));


}




