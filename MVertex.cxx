#include <MVertex.h>
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <vector> 
#include "TRandom.h"
#include "TRandom3.h"
#include "TObject.h"
#include "TNtuple.h"

ClassImp(MVertex)

//----------COSTRUTTORE DI DEFAULT-------------------------------------------------------------------

void MVertex::SaveHisto(TNtuple *Ntupla,TH1F *Differenze){

	Ntupla->Write();
	Differenze->Write();

}


MVertex::MVertex(const char* fname1, const char* fname2, Int_t S, TRandom *Rndm): TObject() {
	//Smearing
	if(S == 0){
		Smearing=0;
	}else{
		Smearing=1;	
	}

	this->LinkGeneraData(fname1);
	this->LinkTrasportaData(fname2);

	Histofile = new TFile("ricostruzione_histo_Class.root","RECREATE");
	
	if(Rndm == 0){	
		GeneratoreEsterno = new TRandom3();
	}else{
		GeneratoreEsterno = Rndm;
	}



	//definisco la geometria dei rivelatori
	R1=4.0;
	R2=7.0;
	srad = 0.015;//radianti finestra
	Extensione = 15.;

	i_1=0;//2 contatori diversi per scorrere gli hit sui layer e fare gli abbinamenti
	i_2=0;
	
	//Debug
	for(UInt_t i = 0; i < 10; i++) CheckPoint[i] = 3;	

	Zvertex = 0;
	
}


//----------------DISTRUTTORE STANDARD-----------------------------------------------------------------

MVertex:: ~MVertex() {	
	delete Generafile;
	delete Origine;
	delete b1;
	delete Trasportafile;
	delete Colpi_1;
	delete Colpi_2;
	delete b3;
	delete b4;
	delete Rumore;
	delete b5;
	delete Zvertex;

	
	
}

void MVertex::LinkGeneraData(const char* fname){
	
	Generafile = new TFile(fname);	 
	Generafile->GetObject("T",Origine);	
	b1=Origine->GetBranch("Event");
	b1->SetAddress(&event.X); //passo l'indirizzo del primo oggetto della struct e assegno tutto a b1
	Neventi = Origine->GetEntries();
}

void MVertex::LinkTrasportaData(const char* fname){

	Trasportafile = new TFile(fname);	

	Trasportafile->GetObject("Layer1",Colpi_1);	
	Trasportafile->GetObject("Layer2",Colpi_2);	

	b3 = Colpi_1->GetBranch("RealLayer1"); 
	b4 = Colpi_2->GetBranch("RealLayer2");  //acquisisco i due branches	

	b3->SetAddress(&lay1.X); //passo l'indirizzo del primo oggetto della struct e assegno tutto a b3
	b4->SetAddress(&lay2.X); // lo stesso per il vettore 

	Npuntilay1 = b3->GetEntries();
	Npuntilay2 = b4->GetEntries();	
	
	Trasportafile->GetObject("Rumore",Rumore);
	b5 = Rumore->GetBranch("Rumore");
	b5->SetAddress(&InfoR.event);
}

void MVertex::GetData(Int_t id_ev){
	//pulisco prima di riempire
	this->ClearStartVector();	

		if(id_ev == 0){
			Colpi_1->GetEvent(i_1);
			Colpi_2->GetEvent(i_2);
		}
	
		while((lay1.Flag == id_ev) & (i_1 < Npuntilay1)){
		
			temp_lay1_z.push_back(lay1.Z);
			temp_lay1_rPhi.push_back(R1*TMath::ATan(lay1.Y/lay1.X));

			++i_1;
			Colpi_1->GetEvent(i_1);
		}

		while((lay2.Flag == id_ev) & (i_2 < Npuntilay2)){
			
			temp_lay2_z.push_back(lay2.Z);
			temp_lay2_rPhi.push_back(R2*TMath::ATan(lay2.Y/lay2.X));

			++i_2;
			Colpi_2->GetEvent(i_2);
		}
}

void MVertex::Associazioni(){
	for(UInt_t i = 0;i < temp_lay1_rPhi.size(); i++){
			
			/////////Smearing sul layer 1	
			if(Smearing){
				temp_lay1_rPhi[i] = temp_lay1_rPhi[i]+GeneratoreEsterno->Gaus(0.,0.003);
				temp_lay1_z[i] = temp_lay1_z[i]+GeneratoreEsterno->Gaus(0.,0.012);				
			}
			
			//Loop sulle possibilità di associazione in base alla finestra azimutale 			
			for(UInt_t u = 0; u < temp_lay2_rPhi.size(); u++){				
				
				///////////Smearing sul layer 2	
				if(Smearing){					
					temp_lay2_rPhi[u] = temp_lay2_rPhi[u]+GeneratoreEsterno->Gaus(0.,0.003);
					temp_lay2_z[u] = temp_lay2_z[u]+GeneratoreEsterno->Gaus(0.,0.012);									
				}			
				
				if(
						(TMath::Abs((temp_lay1_rPhi[i]/R1)	-	(temp_lay2_rPhi[u]/R2)) < srad)
					|	(
							(TMath::Abs((temp_lay1_rPhi[i]/R1		-	temp_lay2_rPhi[u]/R2)) > (2*TMath::Pi()		-	srad/2))
						&	(TMath::Abs((temp_lay1_rPhi[i]/R1		-	temp_lay2_rPhi[u]/R2)) < (2*TMath::Pi()		+	srad/2))
						)
						

					){
						Double_t m = (temp_lay1_rPhi[i]-temp_lay2_rPhi[u])/(temp_lay1_z[i]-temp_lay2_z[u]);
						//trovo ora l'intersezione con (x,y)=0, cioè con l'asse z
						Double_t z = temp_lay1_z[i]-temp_lay1_rPhi[i]/m;
						
						Vertice_costruito.push_back(z);						
				}				
			}		
		}
}

void MVertex::NewGraf(Int_t id_ev){
	if(Zvertex != 0){//evito di distruggere qualcosa prima di averlo creato 
		delete Zvertex;
	}

	this->ClearClusterVector();	
	Indice_peso.clear();
	Indice_peak.clear();

	max_peso = 0;
	multiPeak = 0;
	sameCluster = 0;

	Double_t Zmezzobin =TMath::Abs((ZmaxROI-ZminROI))/(2*Nbin);
	Zvertex = new TH1F(Form("Zricostruito%i",id_ev),"FrequenzeZricostruito",Nbin, ZminROI-Zmezzobin,ZmaxROI-Zmezzobin);
	//Zvertex = new TH1F(Form("Zricostruito%i",id_ev),"FrequenzeZricostruito",Nbin, ZminROI,ZmaxROI);
	
	this->Fill();	
	//Zvertex->Write();
	this->Cluster();

	///////////////indicizzazione e ordinamento///////////////////		
	if(!this->Indicizza(peso_cluster,Indice_peso,kTRUE)){ this->Debug(0, 1); }
	if(!this->Indicizza(f_peak,Indice_peak,kTRUE)) { this->Debug(0, 2); }
	
	//tutte le info sulla clusterizzazione del grafico
        
	//il picco con il cluster piu' grande ha anche l'altezza massima	
	if(Indice_peso[0] == Indice_peak[0] || f_peak[Indice_peak[0]] == f_peak[Indice_peso[0]]) max_peso=1;

	//quanti picchi ci sono con la stessa altezza 
	//multipeak=n°dipicchi-1
	//multipeak=0 picco singolo
	if(f_peak.size() > 1){
		UInt_t hh = 0;
		while(hh < (f_peak.size()-1)){
			if(f_peak[Indice_peak[hh]] == f_peak[Indice_peak[hh+1]]){				
				++multiPeak;
				++hh;
			}else{
				break;
			}					
		}			
	}
	
	//ci sono delle strutture con lo stesso numero di bin
	//sameCluster=0 non ci sono
	//sameCluster=1 ci sono
	if(Indice_peso.size()> 1){
		if(peso_cluster[Indice_peso[0]] == peso_cluster[Indice_peso[1]]) sameCluster=1;
	}

}

void MVertex::Fill(){
	//inserisce i valori di Vertice_costruito nel Grafico nell'intervallo ZminROI e ZMazRoi
	for(UInt_t y=0; y<Vertice_costruito.size();y++){
		//estremi del grafico: [Zvertex->GetBinLowEdge(1) , Zvertex->GetBinLowEdge(Zvertex->GetNbinsX()) + Zvertex->GetBinWidth(1))
		if( (Vertice_costruito[y] >= Zvertex->GetBinLowEdge(1))
			& (Vertice_costruito[y] < Zvertex->GetBinLowEdge(Zvertex->GetNbinsX()) + Zvertex->GetBinWidth(1))
			){					
			Zvertex->Fill(Vertice_costruito[y]);
		}
	}	
}

void MVertex::EvaluateVertex(Int_t min_bin, Int_t max_bin){
	//potrei mettere un controllo che i bin siano >=1 è più piccoli del massimo di bin presenti
	Double_t Num = 0;
	Double_t Den = 0;
	
	Vertex=0;//azzero il valore

	for (Int_t q = min_bin; q <= max_bin; q++ ) {
		Num += Zvertex->GetBinContent(q)*Zvertex->GetBinCenter(q);
		Den += Zvertex->GetBinContent(q);				
	}
	if(Den != 0){
		Vertex = Num/Den;			
	}else{	
		//c'è un problema..		  
		Vertex = Num/0.0000001;
		this->Debug(0, 3);
		}

	
}

void MVertex::Debug(Bool_t yn, Int_t cpt){
	//nuova chance	
	if((yn == 1) & (cpt == 0)){
		//da regolare se si aggiungono altri punti di debug
		for(UInt_t i = 1; i <= 5; i++) CheckPoint[i] = 1;
		SoFarSoGood = 1;
	}
	//checkpoint
	//3 non sono ancora passato
	//1 passato con successo
	//0 passato senza successo

	//id check point
	//0 private azzera il debug
	//1 indicizza peso
	//2 indicizza peak
	//3 evaluateVertex
	//4 errori vari
	//5 meno di due entrate nel grafico
	//6
	//7
	//8
	//9

	if(yn==1){	
		CheckPoint[cpt] = 1;
	}else if(yn==0){
		CheckPoint[cpt] = 0;
		SoFarSoGood = 0;//chance giocata
	}
}
Int_t MVertex::GetError(Int_t cpt){
	if(CheckPoint[cpt] == 0){
		return cpt;
	}else{
		return 0;
	}
}


void MVertex::FindVertex(Int_t id_ev){
	this->GetData(id_ev);
	this->Associazioni();

	//linko al nuovo evento per la Ntuple
	Origine->GetEvent(id_ev);
	Rumore->GetEvent(id_ev);
	
	Molteplicita = event.N;
	Noiselay1 = InfoR.Noiselay1;
	Noiselay2 = InfoR.Noiselay2;
	VerticeVero = event.Z;

	//qua si possono cambiare i metodi e eventualmente impilarli
	this->Debug(1, 0);//nuova chanche
	this->Metodo1(id_ev);
	//this->Metodo2(id_ev);
	
}

void MVertex::Metodo1(Int_t id_ev){

	Int_t ottimi = 0;	
	Int_t no = 0;

	Int_t new_min_roi = 0;
	Int_t new_max_roi = 0;

	
	Nbin = 60;
	ZminROI= -Extensione;
	ZmaxROI= Extensione;
	this->NewGraf(id_ev);

	if(Zvertex->GetEntries() <= 1) this->Debug(0, 5);
	
	for(UInt_t l=0; l<2; l++){

		new_min_roi = 0;
		new_max_roi = 0;

		//il picco con il cluster piu' grande ha anche l'altezza massima e non ci sono altri picchi con la stessa altezza o cluster con la stessa dimensione
		if(max_peso & !multiPeak & !sameCluster) {
			new_min_roi = xmin_bin[Indice_peso[0]]-1;
			new_max_roi = xmax_bin[Indice_peso[0]]+1;

			if(l==0)ottimi++;

		}else if(!max_peso & !multiPeak & (f_peak[Indice_peak[0]] > 1 + f_peak[Indice_peso[0]])){
			
			new_min_roi = xmin_bin[Indice_peak[0]]-1;
			new_max_roi = xmax_bin[Indice_peak[0]]+1;
		
		}else if(multiPeak == 1){
		  //cout<<"MULTIPEAK";

			if(xmin_bin[Indice_peak[1]]<= 2 + xmax_bin[Indice_peak[0]]) {
				
				new_min_roi = xmin_bin[Indice_peak[0]]-1;
				new_max_roi = xmax_bin[Indice_peak[1]]+1;
			}			
		}else{
			//cout<<"STOP NON RICOSTRUITO "<<ev<<" livello "<<l<<" Cause:";						
			//if(max_peso) cout<<" MAX_PESO ";
			//if(multiPeak) cout<<" DOPPIO PICCO ";
			//if(sameCluster) cout<<" DOPPIO CLUSTER ";
			//cout<<endl;
			//cout<<"-----------------------------------------------"<<endl;	
			
			//qualcosa non va..
			this->Debug(0, 4);
			no++;			
		}
		
		if(l == 1){//calcolo il vertice		
		  this->EvaluateVertex(new_min_roi,new_max_roi);			
					
		}else{  //zoom della ROI			
			ZminROI = Zvertex->GetBinLowEdge(new_min_roi);
			ZmaxROI = Zvertex->GetBinLowEdge(new_max_roi) + Zvertex->GetBinWidth(new_max_roi);		
			Nbin = (UInt_t)(ZmaxROI-ZminROI)*20;			
			this->NewGraf(id_ev);
		}		
	}	
}

	
Bool_t MVertex::Indicizza(const vector<Int_t> &vec, vector<Int_t> &index, const Bool_t directionDown){

	if(vec.size()==0) return kFALSE;
	Bool_t ord;
	for(UInt_t i=0;i<vec.size();++i) index.push_back(i);

		do{
		ord=kFALSE;
			for(UInt_t i=0;i<vec.size()-1;++i) {
			Bool_t condizione = directionDown ? vec[index[i]] < vec[index[i+1]] : vec[index[i]] > vec[index[i+1]]; 
				if(condizione){
				ord=kTRUE;
				Int_t tmp = index.at(i+1);
				index[i+1] = index[i];
				index[i] = tmp;
				}
			}
		} while(ord);
	return kTRUE;	
}

void MVertex::Cluster(){
	
	UInt_t temp_nmax_bin=0;	
	UInt_t temp_f_peak=0;
	UInt_t temp_xmin_bin=0;
	UInt_t temp_peso_cluster=0;


	for(Int_t g = 1; g <= Zvertex->GetNbinsX(); g++){			

		//start cluster
		if(g !=1){
			if((Zvertex->GetBinContent(g-1) == 0) & (Zvertex->GetBinContent(g) != 0)){
				
				temp_peso_cluster = 0;
				temp_f_peak = 0;
				temp_nmax_bin = 0;
				temp_xmin_bin = g;//memorizza inizio cluster
			}
			
		}else if(Zvertex->GetBinContent(g) != 0){
		
			temp_peso_cluster = 0;
			temp_f_peak = 0;
			temp_nmax_bin = 0;
			temp_xmin_bin = g;//memorizza inizio cluster									
		}

		//cluster
		if((Zvertex->GetBinContent(g) != 0)){
		temp_peso_cluster += Zvertex->GetBinContent(g);
			if(Zvertex->GetBinContent(g) > temp_f_peak){
				temp_f_peak = Zvertex->GetBinContent(g);
				temp_nmax_bin = g;
			}
		}
		
		//stop cluster
		if(g != Zvertex->GetNbinsX()){

			if((Zvertex->GetBinContent(g) != 0) & (Zvertex->GetBinContent(g+1) == 0)){
				
				//if(temp_f_peak != 1){

			 	//if(temp_peso_cluster != 1){
					//memorizzo
					nmax_bin.push_back(temp_nmax_bin);
					f_peak.push_back(temp_f_peak);
					xmin_bin.push_back(temp_xmin_bin);
					xmax_bin.push_back(g);
					peso_cluster.push_back(temp_peso_cluster);
					
					//azzero le variabili
					temp_nmax_bin=0;
					temp_xmin_bin=0;						
					temp_f_peak=0;
					temp_peso_cluster=0;						

				//}
			}
		
		}else{
			if((temp_f_peak != 0) & (temp_f_peak != 1)){
				//memorizzo
				nmax_bin.push_back(temp_nmax_bin);
				f_peak.push_back(temp_f_peak);
				xmin_bin.push_back(temp_xmin_bin);
				xmax_bin.push_back(g);
				peso_cluster.push_back(temp_peso_cluster);
				
				//azzero le variabili
				temp_nmax_bin=0;
				temp_xmin_bin=0;						
				temp_f_peak=0;
				temp_peso_cluster=0;
			}		
		}
	}
}
void MVertex::ClearStartVector(){
	temp_lay1_z.clear();
	temp_lay1_rPhi.clear();				
	temp_lay2_z.clear();
	temp_lay2_rPhi.clear();
	Vertice_costruito.clear();

}
void MVertex::ClearClusterVector(){
	nmax_bin.clear();
	f_peak.clear();
	xmin_bin.clear();
	xmax_bin.clear();
	peso_cluster.clear();	
}
void MVertex::CloseFile(){
	Generafile->Close();
	Trasportafile->Close();
	Histofile->Close();
	
}

