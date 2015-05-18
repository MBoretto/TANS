void compileExam(TString myopt="fast"){
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
	gSystem->CompileMacro("Analizza.c",opt);
}
