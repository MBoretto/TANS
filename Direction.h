#ifndef DIRECTION_H
#define DIRECTION_H

#include "TObject.h"
#include "TRandom3.h"


class Direction : public TObject
{

public:

	Direction();
	Direction(Double_t phi,Double_t theta);

	virtual ~Direction();
	
	void Setphi(Double_t phi){fPhi = phi;}
	void Settheta(Double_t theta){fTheta = theta;} 

	

	Double_t Getphi()  {return fPhi;} 
	Double_t Gettheta(){return fTheta;}
	Double_t GetNewX() {return fX;}
	Double_t GetNewY() {return fY;}
	Double_t GetNewZ() {return fZ;}


	

	void GeneraHit(Double_t X,Double_t Y,Double_t Z,Double_t R);

	void Scattering(Double_t x,Double_t Xo);

	void SetRNDGenerator(TRandom *gen){
				if(GeneratoreEsterno == 0){
					GeneratoreEsterno = gen;
				}
			};
	void RemoveGenerator(){	
				GeneratoreEsterno = 0;
			}


private:

	TRandom *GeneratoreEsterno;

	Double_t fPhi;
	Double_t fTheta;

	Double_t fX;
	Double_t fY;
	Double_t fZ;

	ClassDef(Direction,1)

};


#endif 


