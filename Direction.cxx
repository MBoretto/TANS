#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Direction.h"
#include "TRandom.h"


ClassImp(Direction)

//________________________________________________________________________
Direction::Direction():TObject(),
 fPhi(0.),
 fTheta(0.){
   // default constructor
	GeneratoreEsterno=0;
 }


//___________________________________________________________________________
Direction::Direction(Double_t Phi,Double_t Theta):TObject(),
fPhi(Phi),
fTheta(Theta){

// Standard constructor
GeneratoreEsterno=0;

}

//___________________________________________________________________________
Direction::~Direction()	 {
  // destructor
}


//_____________________________FUNZIONE GENERAZIONE HIT SUI LAYER_______________________________________________

void Direction::GeneraHit(Double_t X,Double_t Y,Double_t Z,Double_t R) {


Double_t t = -2*TMath::Sin(fTheta)*(
					X*TMath::Cos(fPhi)	+	Y*TMath::Sin(fPhi)
					);

	Double_t tt = TMath::Sqrt(
						(2*TMath::Sin(fTheta)*(
							X*TMath::Cos(fPhi)	+	Y*TMath::Sin(fPhi)))*

						(2*TMath::Sin(fTheta)*(
							X*TMath::Cos(fPhi)	+	Y*TMath::Sin(fPhi)))

						-	4*((TMath::Sin(fTheta)*TMath::Cos(fPhi))*(TMath::Sin(fTheta)*TMath::Cos(fPhi))
						+	(TMath::Sin(fTheta)*TMath::Sin(fPhi))*
							(TMath::Sin(fTheta)*TMath::Sin(fPhi)))*

							(X*X+Y*Y-R*R)
						
				); 



	Double_t den = 2*(TMath::Sin(fTheta)*TMath::Sin(fTheta)
				);


			Double_t t1=(t+tt)/den;
			
		
	fX=X + t1*(TMath::Sin(fTheta)*TMath::Cos(fPhi));
	fY=Y + t1*(TMath::Sin(fTheta)*TMath::Sin(fPhi));
	fZ=Z + t1*TMath::Cos(fTheta);

} 



//___________________________FUNZIONE SCATTERING_________________________________________________
void Direction::Scattering(Double_t x,Double_t Xo){


	//      x 		spessore del materiale 
	//      Xo 		lunghezza di radiazione
	
	Double_t Beta = 1.;
	//Double_t c = 3*10^8;	//Velocità luce -> 1 in unità naturali
	Double_t p = 750.;	//impulso particella
	Int_t nZ = 1;		//carica della particella


	
	//calcolo RMS della distribuzione del nuovo angolo theta
	Double_t Theta0 = (13.6/(Beta*p))	
				*	nZ	
				*	TMath::Sqrt(x	/	Xo)
				*(1	+ 0.038	*	TMath::Log(x/Xo));




	//Angoli dal vertice
	//fTheta
	//fPhi
	
	//angoli dal multiple scattering
	Double_t php = GeneratoreEsterno->Uniform(0.,2*TMath::Pi());
	Double_t thp = GeneratoreEsterno->Gaus(0.,Theta0);
	



	//Matrici per il prodotto Matriciale
	Double_t mr[3][3];	//matrice rotazione
	Double_t cdp[3];	//matrice nuova direzione
	Double_t cd[3];		//componenti nuovo vettore nel SR del lab

	//definizione delle matrici
	mr[0][0]=-TMath::Sin(fPhi);
	mr[1][0]=TMath::Cos(fPhi);
	mr[2][0]=0.;
	mr[0][1]=-TMath::Cos(fPhi)*TMath::Cos(fTheta);
	mr[1][1]=-TMath::Cos(fTheta)*TMath::Sin(fPhi);
	mr[2][1]=TMath::Sin(fTheta);
	mr[0][2]=TMath::Sin(fTheta)*TMath::Cos(fPhi);
	mr[1][2]=TMath::Sin(fTheta)*TMath::Sin(fPhi);
	mr[2][2]=TMath::Cos(fTheta);

	
	cdp[0] = TMath::Sin(thp)*TMath::Cos(php);
	cdp[1] = TMath::Sin(thp)*TMath::Sin(php);
	cdp[2] = TMath::Cos(thp);

	
	//prodotto matriciale
	for(Int_t i=0; i<3; i++){
		cd[i]=0.;
		for(Int_t j=0; j<3; j++){
			cd[i]+=mr[i][j]*cdp[j];
		}

	}

	//output


	fTheta = TMath::ACos(cd[2]);
	
	Double_t ang = TMath::ATan2(cd[1],cd[0]);

	fPhi=ang;


}





