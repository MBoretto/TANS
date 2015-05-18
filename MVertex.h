#ifndef FATO_H
#define FATO_H

#include "TObject.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector> 
#include "TNtuple.h"

typedef struct {
		Double_t X,Y,Z;
		Int_t N;
}SINGLE_EVENT;

typedef struct {
    Double_t X,Y,Z;
    Int_t Flag;		
} HIT;

typedef struct {
  Int_t event;
  Int_t tipo;
  Int_t Noiselay1;
  Int_t Noiselay2;
} infoRumore;



class MVertex : public TObject {

public: 	
	MVertex(const char* fname1, const char* fname2,Int_t S, TRandom *Rndm);
	virtual ~MVertex();
	//funzioni pubbliche	
	void FindVertex(Int_t id_ev);
	
	//get inline per analisi dati
	Int_t GetNevent(){return Neventi;};
	Int_t GetNpuntiLay1(){return Npuntilay1;};
	Int_t GetNpuntiLay2(){return Npuntilay2;};
	Int_t GetError(Int_t cpt);
	Double_t GetVertex(){return Vertex;};
	Int_t GetMolteplicita(){return Molteplicita;};
	Int_t GetNoiselay1(){return Noiselay1;};	
	Int_t GetNoiselay2(){return Noiselay2;};
	Double_t GetVerticeVero(){return VerticeVero;};
	Int_t GetSoFarSoGood(){return SoFarSoGood;};

	TRandom * GetInternalGenarator(){return GeneratoreEsterno;};
	void SaveHisto(TNtuple *Ntupla,TH1F *Differenze);


private:
	///////////////
	//funzioni private
	///////////////
	void LinkGeneraData(const char* fname);
	void LinkTrasportaData(const char* fname);
	void GetData(Int_t id_ev);
	void Associazioni();
	void ClearStartVector();
	void Metodo1(Int_t id_ev);
	void NewGraf(Int_t id_ev);
	void Fill();
	void Cluster();
	Bool_t Indicizza(const vector<Int_t> &vec, vector<Int_t> &index, const Bool_t directionDown);
	void ClearClusterVector();
	void EvaluateVertex(Int_t min_bin, Int_t max_bin);
	void Debug(Bool_t yn, Int_t checkpoint);	
	void CloseFile();

	///////////////
	//Variabili Condivise nella classe
	///////////////
	Bool_t Smearing;
	Int_t Neventi;
	Int_t Npuntilay1;//Nb3
	Int_t Npuntilay2;//Nb4
	Double_t srad;
	Double_t Extensione;
	
	//Generatore Esterno
	TRandom *GeneratoreEsterno;

	//Dati da passare alla Ntupla
	Int_t Molteplicita;
	Int_t Noiselay1;
	Int_t Noiselay2;
	Double_t VerticeVero;
	
	//Debug
	Bool_t SoFarSoGood; //Flag per sancire se l'evento è stato ricostruiro o ci sono stati dei problemi 1 ok 0 non ricostruito
	Int_t CheckPoint[10];	     

	//Geometria dei rivelatori
	Double_t R1;
	Double_t R2;

	//Link Dati Genera	
	TFile *Generafile;//k
		TTree *Origine;
		TBranch *b1;
		SINGLE_EVENT event; 
	
	//Link Dati Trasporta
	TFile *Trasportafile;//h
	
		//urti
		TTree *Colpi_1;
		TTree *Colpi_2;
		TBranch *b3;
		TBranch *b4;
		HIT lay1;  
		HIT lay2;

		//rumore
		TTree *Rumore;
		TBranch *b5;
		infoRumore InfoR;

	//doppio indice per i layer indispensabile per leggere tutti gli eventi
	Int_t i_1;//2 contatori diversi per scorrere gli hit sui layer e fare gli abbinamenti
	Int_t i_2;

	//Acquisizione Dati
	vector<Double_t> temp_lay1_z;
	vector<Double_t> temp_lay1_rPhi;
	vector<Double_t> temp_lay2_z;
	vector<Double_t> temp_lay2_rPhi;
	
	
	//Associazione
	//z
	vector<Double_t>  Vertice_costruito;

	//FindVertex
	TH1F *Zvertex;
	Int_t Nbin;
	Double_t ZminROI;
	Double_t ZmaxROI;
	Double_t Vertex;
	
	//Cluster
	vector<Int_t> nmax_bin;
	vector<Int_t> f_peak;
	vector<Int_t> xmin_bin;
	vector<Int_t> xmax_bin;
	vector<Int_t> peso_cluster;

	//ordine Cluster
	vector<Int_t> Indice_peak;
	vector<Int_t> Indice_peso;
	//info sui cluster
	Bool_t max_peso;
	Int_t multiPeak;
	Bool_t sameCluster;

	//Per memorizzare i grafici
	TFile *Histofile;
	

ClassDef(MVertex,1);
};


#endif
