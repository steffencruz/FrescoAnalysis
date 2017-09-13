#ifndef TFRESCOANALYSIS_H
#define TFRESCOANALYSIS_H

#include <stdlib.h>
#include <vector>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TList.h>
#include <Rtypes.h>
#include <TROOT.h>

#ifndef FADDIR
#define FADDIR "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/FrescoAnalysis"
#endif

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
//  ____TFrescoAnalysis____
//  
//  Uses FRESCO template files and neutron state information to 
//  generate a DWBA calculation and fit it to data to extract SF
//  
//  Can also be used to generate elastic scattering files to 
//  normalize data to DWBA
//
//  _____SearchMode_____
//
//  Calculates SF for all possible neutron states and returns
//  lowest chi2 fit
//
//  _____ScanMode_____
//
//  Calculates SF & Sig for all possible neutron states and 
//  produces a graphs of these values as a function of ExcHi and
//  fits this to a cumulative function which indicates the best SF.
//
//
//  ___Description of variables ___
/*

  REAC - "dp", "pp", "dd" selects the reaction
  ** neutron state information **  
  EXC - Excitation energy [MeV]
  BE  - Binding energy of neutron [MeV]
  JF  - Total angular momentum of final state, coupling [s1/2][nlj]
  J0  - Total angular momentum of orbital that neutron is placed in
  L0  - Orbital angular momentum of orbital that neutron is placed in
  N0  - Principal quantum number of orbital that neutron is placed in
  SF  - Spectroscopic Factor
  CS - Total cross section [with SF appplied]

  ** required to determine spectroscopic factor error **
  NORM - Normalization Used
  GAM  - Gamma gate used 
  ** these are taken from the Results_*Exc*.root file **

  ExcHi - Upper excitation energy window, see ScanEnergyWindow

  Verbose : -
 -> 0   =   No Printing 
 -> 1   =   Only High level summary eg. CheckFeeding
 -> 2   =   Standard printing
 -> 3   =   Print everything
*/
//
// ___Notes___
/*
 Further work : -
   -> Add (d,t) reaction
   -> Copy files and make directories?
   -> Speed up the process of doing calculations. Scripts.
   
 Known bugs :-
*/  
//
// Made By Steffen Cruz, 2016
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

class TFrescoAnalysis 	{
	
	public:

		TFrescoAnalysis(void);
		~TFrescoAnalysis();		
		
		static void Print(Option_t * = "");			
		static void Clear(Option_t * = "");	
		
		static TGraphErrors *GetData(){return gdata;}	
		static TGraph *GetDWBA(){return gcalc;}
   
    static Double_t Norm(){ return NORM;}
    static Double_t NormErr(){ return NORM_ERR;}
        
    static Bool_t CalculateSF(Bool_t usechi2=true, Double_t br=1.0, Double_t br_err=0.0);   
    static TGraph *ChiSquaredGraph(Double_t val=1.0, Double_t dval=0.01, TGraphErrors *gd=NULL, TGraph *gc=NULL);	
    static Double_t ChiSquared(Double_t scale=1.0, TGraphErrors *gd=NULL, TGraph *gc=NULL);
    static Double_t ChiSquaredError(Double_t dval=0.01, TGraphErrors *gd=NULL, TGraph *gc=NULL);
    		
    static Bool_t SetInfo(UInt_t A, std::string reac="dp", std::string om="PP", 
        Double_t exc=0.0, UInt_t l=0, Double_t jo=-1, Double_t jf=-1,std::string fname="");
        
    static Bool_t SetFileName(std::string fname);
    static std::string GetFileName(){ return FNAME; }
    static void SetVerbose(UInt_t verb) { verbose = verb; }
      
    // use reaction info to make fresco files                			   
    static Bool_t MakeFile(std::string file_opt = "frin+search+min");
    
    static TGraph *CalculateAngDist(Double_t sf=1.0, UInt_t colour=1, Bool_t sine=false);
    static TList *FitData(Bool_t draw=true, Bool_t write=false);
    
//    static Bool_t TestOpticalModels();
    static TList *StateSearch(UInt_t A, Double_t exc, std::string dname, std::string om="PP",Bool_t all=false);

    static Bool_t SaveResults(std::string rootfile_name, std::string dir);
    
    static Bool_t DoFrescoAnalysis();

    // get OM
    static Double_t GetPotVal(std::string fname);
    static UInt_t  ReadPotVals(std::string fname); // parse from a file (in case of OM fit)
    
    static UInt_t  WritePotVals(std::string fname); 
    // set OM
    static UInt_t SetPotVals(UInt_t kbpot=0);

    static Bool_t UserSetPotVal(std::string parname, Double_t val);       
    static Bool_t UserSetPotVal(UInt_t kbpot, UInt_t type, UInt_t parnum, Double_t val);   
    // scan OM
    static TCanvas *UserScanPotVal(std::string parname, Double_t valmin, Double_t valmax, UInt_t nsteps);    
    static TCanvas *UserScanPotVal(UInt_t kbpot, UInt_t kbtype, UInt_t parnum, 
                                  Double_t valmin, Double_t valmax, UInt_t nsteps);
    static Bool_t UserScanPotVal2D(std::string var1="R0", Double_t val1lo=0.8, Double_t val1hi=1.2, UInt_t ns1=5,   
                                   std::string var2="RD", Double_t val2lo=1.2, Double_t val2hi=1.6, UInt_t ns2=5);                                  
    
    // control the sfresco fit
    static UInt_t AddFitData(std::string dfile, Double_t thmin=0, Double_t thmax=0);
    static UInt_t GetNFitData(){ return NFITDATA; }             
    
    static Bool_t AddFitPar(UInt_t kind, UInt_t kbpot=1, UInt_t type=1, UInt_t parnum=1, 
                            Double_t val=-1, Double_t valmin=-1, Double_t valmax=-1, Double_t vstep=-1);
    static Bool_t AddFitPar(std::string,Double_t val=-1, Double_t valmin=-1, Double_t valmax=-1, Double_t vstep=-1);
    static UInt_t GetNFitPars(){ return NFITPARS; }
   
    // get results from root file                         
    static Bool_t ExtractAngDistInfo(std::string dfile);
    // get total cross section from fort.13
    static Bool_t ExtractTotalCrossSection(const char *fort_name);

    static std::string GetParName(UInt_t kbpot, UInt_t type, UInt_t parnum, Bool_t holdername=true);
    static Bool_t GetParInfo(std::string parname, UInt_t &kind, UInt_t &kbpot, UInt_t &kbtype, UInt_t &parnum);    
    static UInt_t GetAllowedTransfers(UInt_t A, std::vector<int> &lo, std::vector<double> &jo,
           std::vector<double> &jf, Bool_t all=true);

    static TList *MakePotentials(Int_t kbpot);
    
  private:

    static Bool_t CheckInfo();   
  
    static Bool_t SetParticles(UInt_t A, Double_t exc, std::string reac);
    static Bool_t SetOM(std::string reac, std::string om);
    static Bool_t SetTransfer(Double_t exc, UInt_t l, Double_t jo, Double_t jf);  
  
    static UInt_t BuildReplacementStrings();
    
    static Bool_t CreateFile(const char *template_name, const char *file_name);  
    static void ReplaceHolders(std::string &line);
    static Int_t GetHolderIndex(std::string hname);
  
    static Bool_t ExtractFitResult(std::string line);    
     
    static Double_t Coulomb(Double_t x, Double_t R);
    static Double_t WS_Volume(Double_t x, Double_t V, Double_t R, Double_t A);
    static Double_t WS_Surface(Double_t x, Double_t V, Double_t R, Double_t A);  
  
    static void GetRid(const char *name, Bool_t delete_all=false);
    static Bool_t CheckState(TGraph *g);  
    
      
  private: 
  
    static UInt_t verbose;  // how much stuff gets printed

    static Bool_t set_info; // flag which indicates that everything is set up
    static Int_t TYPE;      // reaction type, same as Easy Trees
    
    static std::vector<std::string> HOLDERNAME;
    static std::vector<std::string> HOLDERVAL;  
    
    // file info
    static std::string FNAME; // base name for all fresco files
    static std::string DNAME; // experimental ang. dist. data file
    static std::string RNAME; // experimental ang. dist. root file

    // beam info
    static UInt_t BEAMA;
    static Double_t BMASS;
    static Double_t BSPIN;
  
    static Double_t EPERU;
    static Double_t ELAB;
  
    // reaction 
    static std::string REAC;    
    static Double_t QVAL;
    
    // optical model
    static std::string OM;  // full input 
    static std::string OM1; // entrance elastic channel   
    static std::string OM2; // exit elastic channel           
    
    // light target info
    static UInt_t TARGA;
    static Double_t TMASS;
    static Double_t TSPIN;
    
    // heavy recoil info
    static UInt_t RECOA;
    static Double_t RMASS;
    static Double_t REXC; // excited state energy [MeV]
    static Double_t RSN; // one neutron separation energy [MeV]

    // transferred neutron info
    static Double_t JO;
    static Double_t JF;
    static UInt_t L;
    static UInt_t NO;
    static Double_t NBE; // neutron binding energy [MeV]
    
    // how many potentials, how many parameters and data sets to fit
    static UInt_t NFITPARS;
    static UInt_t NFITDATA;    
    static UInt_t NKBPOT;        

    // experimental angular distribution information
    static Double_t NORM;
    static Double_t NORM_ERR;
    static Double_t GAM;
    static Double_t GAMEFF;
    static Double_t GAMEFF_ERR;
    static Double_t SA;
    static Double_t SA_ERR;     
    static Double_t SF;
    static Double_t SF_ERR;
    static Double_t CS;
    static Double_t CS_ERR;  
    static Double_t CHI2;
  
    static Double_t EXCHI;  

    static TGraphErrors *gdata;
    static TGraph *gcalc;
    
	ClassDef(TFrescoAnalysis,0)
};


#endif