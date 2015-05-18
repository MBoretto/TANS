void executeExam(TString myopt="fast"){
  Int_t n = 0;
  Int_t esp = 0;
  Int_t molt = 0;
  Int_t zeta = 0;
  Bool_t scatt = 0;
  Int_t rhum = 0;
  Bool_t smear = 0;

  ////PARAMETRI DELLA SIMULAZIONE////

  cout<<"INSERIRE I PARAMETRI PER LA GENERAZIONE : "<<endl;
  cout<<"            Inserire il numero di eventi da simulare nella forma n*10^(esp) : "<<endl;
  cout<<" -> Scegliere n: "<<endl;
  cin>>n;
  cout<<" -> Scegliere esp: "<<endl;
  cin>>esp;
  cout<<" -> Scegliere la molteplicità degli eventi : 1->molteplicità da kinem.root, n->molteplicità fissa n"<<endl;
  cin>>molt;
  cout<<" -> Scegliere la posizione della z dei vertici : -1->zeta con distribuzione gaussiana, n->z in posizione fissa a n cm"<<endl;
  cin>>zeta;
  cout<<endl;
  cout<<"INSERIRE I PARAMETRI PER IL TRASPORTO : "<<endl;
  cout<<" -> Attivare lo scattering multiplo? 1->attivo, 0->non attivo :"<<endl;
  cin>>scatt;
  cout<<" -> Attivare il rumore? 1->attivo con molteplicità gaussiana, 0->non attivo, n->attivo con molteplicità fissa n su ogni layer : "<<endl;
  cin>>rhum;
  cout<<endl;
  cout<<"INSERIRE I PARAMETRI PER LA RICOSTRUZIONE : "<<endl;
  cout<<" -> Attivare lo smearing sui layer? 1->attivo, 0->non attivo : "<<endl;
  cin>>smear;

  
  TStopwatch totale;
  totale.Start(kTRUE);

  char * opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
	gSystem->CompileMacro("Fato.cxx",opt);
	gSystem->CompileMacro("Direction.cxx",opt);
	gSystem->CompileMacro("MVertex.cxx",opt);
	gSystem->CompileMacro("Genera.c",opt);
	gSystem->CompileMacro("Trasporta.c",opt);
	gSystem->CompileMacro("Ricostruisci.c",opt);
	if( (molt==1) & (zeta==-1) & (scatt==1) & (rhum==1) & (smear==1)) gSystem->CompileMacro("Analizza.c",opt);
	

	delete gRandom;
	gRandom = new TRandom3();
	gRandom->SetSeed(45623);

		
	Genera(n,esp,molt,zeta,gRandom);
	Trasporta(scatt,rhum,gRandom);
	Ricostruisci(smear,gRandom);
	//delete Generatore;
	if( (molt==1) & (zeta==-1) & (scatt==1) & (rhum==1) & (smear==1)) Analizza(molt,zeta,rhum);

	totale.Stop();

	cout<<endl;
	cout<<"SIMULAZIONE COMPLETATA CON SUCCESSO!!!"<<endl;
	cout<<"L'intera simulazione è durata : "<<endl;
	totale.Print();
	cout<<endl<<endl;



}
