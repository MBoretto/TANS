#ifndef FATO_H
#define FATO_H

#include "TRandom.h"
#include "TRandom3.h"
#include "TH1F.h"


//-------------------------Classe per la generazione di numeri casuali----------------------//

class Fato : public TRandom3 {


public: 
	Fato(TRandom * RNDM);
	Fato(Int_t k, TRandom * RNDM);
	virtual ~Fato();

	
	void NParticelle(Int_t &N);	
	void Direzioni(Double_t &phi, Double_t &theta);
         //Double_t VertZ() {return this->Gaus(0.,5.3) ;}
	//Double_t VertX_Y() {return this->Gaus(0.,0.01) ;}
	
	Double_t VertZ();
	Double_t VertX_Y();
	
private:

	Int_t fdistribuzione;
	TH1F *rap;
	TH1F *mol;
	TRandom *GeneratoreEsterno;


ClassDef(Fato,1);
};


#endif
