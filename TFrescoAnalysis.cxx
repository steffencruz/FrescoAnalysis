#include "TFrescoAnalysis.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <Rtypes.h>
#include <TF1.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2F.h>
#include <TROOT.h>

ClassImp(TFrescoAnalysis)

Bool_t   TFrescoAnalysis::set_info = false;
Int_t    TFrescoAnalysis::TYPE     = -1;
UInt_t   TFrescoAnalysis::verbose=2;

UInt_t   TFrescoAnalysis::BEAMA=0;
UInt_t   TFrescoAnalysis::RECOA=0;
UInt_t   TFrescoAnalysis::TARGA=0;

Double_t TFrescoAnalysis::BMASS=0.0;
Double_t TFrescoAnalysis::BSPIN=0.0;

Double_t TFrescoAnalysis::EPERU=0.0;
Double_t TFrescoAnalysis::ELAB=0.0;
Double_t TFrescoAnalysis::QVAL=0.0;

Double_t TFrescoAnalysis::TMASS=0.0;
Double_t TFrescoAnalysis::TSPIN=0.0;
Double_t TFrescoAnalysis::RMASS=0.0;
Double_t TFrescoAnalysis::REXC=0.0;
Double_t TFrescoAnalysis::RSN=0.0;
Double_t TFrescoAnalysis::NBE=0.0;  
Double_t TFrescoAnalysis::JF=0.0;
UInt_t   TFrescoAnalysis::L=0;
Double_t TFrescoAnalysis::JO=0.0;    
UInt_t   TFrescoAnalysis::NO=0.0;    
Double_t TFrescoAnalysis::EXCHI=0.0;   
Double_t TFrescoAnalysis::NORM=0.0;
Double_t TFrescoAnalysis::NORM_ERR=0.0;
Double_t TFrescoAnalysis::GAM=0.0;
Double_t TFrescoAnalysis::GAMEFF=0.0;
Double_t TFrescoAnalysis::GAMEFF_ERR=0.0;     
Double_t TFrescoAnalysis::SA=0.0;
Double_t TFrescoAnalysis::SA_ERR=0.0;
Double_t TFrescoAnalysis::SF=0.0;
Double_t TFrescoAnalysis::SF_ERR=0.0;  
Double_t TFrescoAnalysis::CS=0.0;
Double_t TFrescoAnalysis::CS_ERR=0.0;
Double_t TFrescoAnalysis::CHI2=0.0;  

std::string TFrescoAnalysis::OM   = "";
std::string TFrescoAnalysis::OM1  = "";
std::string TFrescoAnalysis::OM2  = "";
std::string TFrescoAnalysis::REAC = "";  

std::vector<std::string> TFrescoAnalysis::HOLDERNAME;
std::vector<std::string> TFrescoAnalysis::HOLDERVAL;
UInt_t   TFrescoAnalysis::NFITDATA = 0;
UInt_t   TFrescoAnalysis::NFITPARS = 0;
UInt_t   TFrescoAnalysis::NKBPOT   = 0;
  
std::string TFrescoAnalysis::FNAME = "";    
std::string TFrescoAnalysis::DNAME = "";    
std::string TFrescoAnalysis::RNAME = "";    

TGraphErrors *TFrescoAnalysis::gdata;
TGraph *TFrescoAnalysis::gcalc;

/////////////   /////////////   /////////////   /////////////   /////////////   /////////////

static const unsigned long npos = std::string::npos;	

// data directories
std::string BASEDIR = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools";
std::string FDIR = BASEDIR+"/FrescoAnalysis";
std::string TDIR = FDIR+"/templates";
std::string MDIR = FDIR+"/opticalmodel";
std::string DDIR = BASEDIR+"/AngularDistribution";

// exectuables
std::string FREX = BASEDIR+"/fresco/fres/source/fresco";
std::string SFREX = BASEDIR+"/fresco/fres/source/sfresco";

/////////////   /////////////   /////////////   /////////////   /////////////   /////////////

TFrescoAnalysis::TFrescoAnalysis()	{	}

TFrescoAnalysis::~TFrescoAnalysis()	{	}

TList *TFrescoAnalysis::StateSearch(UInt_t A, Double_t exc, std::string dname, std::string om, Bool_t all){

  std::vector<int> lo;
  std::vector<double> jo, jf; 
  
  UInt_t n = GetAllowedTransfers(A,lo,jo,jf,all);
  TList *list;
  
  if(!n)
    return list;
    
  list = new TList();    
    
  std::vector<double> sa, sa_err, sf, sf_err, sig, sig_err, chi2;    
  Double_t chi2_best = 1000.0;
  UInt_t indx = 0;
  
  TGraph *gthry[n];

  for(int i=0; i<n; i++){
    if(verbose) printf("\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    // make fresco files, locate fit data and root file with gamma and norm info
    if(!SetInfo(A,"dp",om,exc,lo[i],jo[i],jf[i]) || !AddFitData(dname))
      return list;
      
    // now fit it!!
    TList *fitlist = FitData(false);
    if(!fitlist || fitlist->GetEntries()!=2)
      return list;

    if(!i){ // only add exp data once
      gdata = (TGraphErrors*) fitlist->FindObject(Form("%s_ExpData",FNAME.c_str()));    
      list->Add(gdata);
    }
    gthry[i] = (TGraph*) fitlist->FindObject(Form("%s_Thry",FNAME.c_str()));
    list->Add(gthry[i]);

    sa.push_back(SF);  
    sa_err.push_back(SF_ERR);      
    sf.push_back(SF);  
    sf_err.push_back(SF_ERR);  
    sig.push_back(CS);
    sig_err.push_back(CS_ERR);
    chi2.push_back(CHI2);
    
    if(CHI2<chi2_best){
      indx = i;
      chi2_best = CHI2;
    }
  }
  
  std::string cname = Form("Sr%i_dp%.0f",A,exc*1e3);
  std::string ctitle = Form("%iSr(d,p) : %iSr E_{exc} = %.0f keV",A,A+1,exc*1e3);

  GetRid(cname.c_str());
  TCanvas *c = new TCanvas(cname.c_str(),ctitle.c_str(),800,600);
  TLegend *leg = new TLegend(0.62,0.63,0.97,0.9);
  
  // draw invisible graph to make sure we can change axes
  /*
  TGraph *gtmp = (TGraphErrors *)gdata->Clone("tmp");

  Double_t gmin = 0.5*TMath::MinElement(gtmp->GetN(),gtmp->GetY());
  Double_t gmax = 1.5*TMath::MaxElement(gtmp->GetN(),gtmp->GetY());
  gtmp->SetPoint(gtmp->GetN()-1,180,gmax);  
  gtmp->GetXaxis()->SetRangeUser(0,110);  
  gtmp->GetYaxis()->SetRangeUser(gmin,gmax); 
  gtmp->SetLineColor(0);
  gtmp->Draw("AL");
  */
  gdata->SetTitle(Form("^{%i}Sr(%c,%c): ^{%i}Sr E_{exc} = %.0f keV; #theta_{CM} [#circ]; #frac{d#sigma}{d#Omega} [mb/sr]",A,REAC.at(0),REAC.at(1),A+1,exc*1e3));
  
  leg->AddEntry(gdata,"Experimental Data","lp");    
  gdata->Draw("A P");          
  
  std::string fname_base = Form("StateSearch_Sr%i_dp%.0f",A,exc*1e3);
  TFile *fout = new TFile(Form("%s.root",fname_base.c_str()),"RECREATE");
  gdata->Write();
  
  std::string msgstr;
  std::ofstream outfile(Form("%s.txt",fname_base.c_str()));
  if(verbose){
    printf("\n\n\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
    printf("   Fit summary for %iSr %.3f MeV State: -\n",RECOA,exc);
  }
  for(int i=0; i<n; i++){

    // set nice display styles
    if(lo[i]==0){ // s1/2 is red
      gthry[i]->SetMarkerColor(2);
      gthry[i]->SetLineColor(2);     
      gthry[i]->SetLineStyle(2);     
    } else if(lo[i]==2){
      if(jo[i]==1.5){ // d3/2 is nice green
        gthry[i]->SetMarkerColor(8);
        gthry[i]->SetLineColor(8);          
        if(jf[i]==1 || jf[i]==jo[i]) // J=1+ is dashed line
          gthry[i]->SetLineStyle(2); 
        if(jf[i]==2) // J=2+ is dotted line
          gthry[i]->SetLineStyle(3);    
      } else if(jo[i]==2.5){// d5/2 is cyan
        gthry[i]->SetMarkerColor(7);
        gthry[i]->SetLineColor(7);          
        if(jf[i]==2 || jf[i]==jo[i]) // J=2+ is dashed line
          gthry[i]->SetLineStyle(2); 
        if(jf[i]==3) // J=3+ is dotted line
          gthry[i]->SetLineStyle(3);    
      }            
    } else if(lo[i]==4){ // g7/2 is blue
      gthry[i]->SetMarkerColor(4);
      gthry[i]->SetLineColor(4);                            
      if(jf[i]==3 || jf[i]==3.5) // J=3+ is dashed line
        gthry[i]->SetLineStyle(2); 
      if(jf[i]==4) // J=4+ is dotted line
        gthry[i]->SetLineStyle(3);         
    } 
    // set the best fit result curve to have a solid line    
    if(i==indx)
      gthry[i]->SetLineStyle(1);
      
    if(A==94 || A==96)
      leg->AddEntry(gthry[i],Form("L=%i Jo=#frac{%.0f}{2}^{+} Jf=#frac{%.0f}{2}^{+} #chi^{2}/N=%.2f ",lo[i],2*jo[i],2*jf[i],chi2[i]),"l");        
    else if(A==95)
      leg->AddEntry(gthry[i],Form("L=%i Jo=#frac{%.0f}{2}^{+} Jf=%.0f^{+} #chi^{2}/N=%.2f ",lo[i],2*jo[i],jf[i],chi2[i]),"l");        
    gthry[i]->Draw("same L");
  
    msgstr.assign(Form("%i. I=%.1f x [ L=%i , j=%.0f/2 ]-> J=%.1f\t SF = %.4f +/- %.4f   Sig(tot) = %.3f +/- %.3f   Chi2 = %.2f",
    i,BSPIN,lo[i],2*jo[i],jf[i],sf.at(i),sf_err.at(i),sig.at(i),sig_err.at(i),chi2.at(i)));
    if(verbose) {printf("\n%s",msgstr.c_str()); fflush(stdout); }
    outfile << msgstr.c_str() << "\n";
    
    gthry[i]->Write();
  }  
  leg->Draw();    
  outfile.close();  
  if(verbose) printf("\n\n   Saved results to ' %s '.",fout->GetName());  
  if(verbose) printf("\n\n * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");

  gPad->SetLogy(true);
  c->Print(Form("%s.pdf",cname.c_str()));  
  c->Write();
 // fout->Close();
  
  // use the best fit to set the current status
  SetInfo(A,"dp",om,exc,lo.at(indx),jo.at(indx),jf.at(indx));
  SetVerbose(0);  
  AddFitData(dname);
  SetVerbose(2);

  SA = sa.at(indx);
  SA_ERR = sa_err.at(indx);  
  SF = sf.at(indx);
  SF_ERR = sf_err.at(indx);
  CS = sig.at(indx);
  CS_ERR = sig_err.at(indx);  
  CHI2 = chi2_best;
  
  gcalc = gthry[indx];
  Print();
  
  return list;
}

Double_t TFrescoAnalysis::ChiSquaredError(Double_t dval, TGraphErrors *gd, TGraph *gc){

  TGraph *g = ChiSquaredGraph(1.0,dval,gd,gc);
  
  Double_t *x=g->GetX();
  Double_t *y=g->GetY();  
  Double_t chi2_min = TMath::MinElement(g->GetN(),y);
  Double_t xmin=100, xmax=0;
  
  for(int i=0; i<g->GetN(); i++){

    if(y[i]<=chi2_min+1){
      if(x[i]<xmin)
        xmin=x[i];
      if(x[i]>xmax)
        xmax=x[i];
    }  
  }
  Double_t xerr = 0.5*(xmax-xmin);
  printf("\n\n Chi2 Fit : Minimum chi2 = %.3f, estimated rel error = %.3f\n",chi2_min,xerr);
  return xerr;
}

TGraph *TFrescoAnalysis::ChiSquaredGraph(Double_t val, Double_t dval, TGraphErrors *gd, TGraph *gc){

  if(ChiSquared(val,gd,gc)==0)
    return 0;
    
  TGraph *gchi2 = new TGraph();
  
  Double_t chi2, vtmp;
  UInt_t imax=100; // scan the chi squared value for val±50%
  for(int i=0; i<imax; i++){
    
    vtmp = val+dval*((double)i-(double)imax/2);
    chi2 = ChiSquared(vtmp,gd,gc);
    gchi2->SetPoint(i,vtmp,chi2);  
  }
  
  return gchi2;
}

Double_t TFrescoAnalysis::ChiSquared(Double_t scale, TGraphErrors *gd, TGraph *gc){

  if(!gd && gdata)
    gd = gdata;
  if(!gc && gcalc)
    gc = gcalc;
    
  if(!gd || !gc){
    printf("\n\n\t Error :  Graphs must first be made. gdata=%p, gcalc=%p\n\n",gd,gc);
    return 0;
  }
  Double_t chi2=0;
  Double_t *x = gdata->GetX();
  Double_t *y = gdata->GetY(); 
  Double_t *ey = gdata->GetEY();   
  
  for(int i=0; i<gdata->GetN(); i++){ 
    chi2 += pow(y[i]-scale*gcalc->Eval(x[i]),2)/pow(ey[i],2);
  }

  return chi2/gdata->GetN();
} 

Bool_t TFrescoAnalysis::SetInfo(UInt_t A, std::string reac, std::string om, Double_t exc, UInt_t l, Double_t jo, Double_t jf, std::string fname){

  Clear();
  
  if(!SetParticles(A,exc,reac))
    return false;

  if(!SetOM(reac,om))
    return false;
    
  if(TYPE==0 || TYPE==4){
  
    if(!SetTransfer(exc,l,jo,jf))
      return false;
        
    if(verbose>1) 
      printf("\n %iSr(%c,%c) : Exc = %.3f MeV,  I=%.1f x [ L=%i , j=%.0f/2 ]-> J=%.1f, OM = '%s'[d]+'%s'[p]: -\n",
            A,REAC.at(0),REAC.at(1),REXC,BSPIN,L,2*JO,JF,OM1.c_str(),OM2.c_str());  

    if(!fname.length())
      FNAME = Form("sr%i_%s%.0f_L%ij%.1fJ%.1f_p%s_d%s",A,reac.c_str(),REXC*1e3,L,JO,JF,OM1.c_str(),OM2.c_str());
    else SetFileName(fname);

    if(EXCHI){
      DNAME = Form("AngDist_%.0f_ExcHi%.0f.txt",REXC*1e3,EXCHI);
      RNAME = Form("Results_%.0f_ExcHi%.0f.root",REXC*1e3,EXCHI);
    } else {
      DNAME = Form("AngDist_%s%.0f.txt",REAC.c_str(),REXC*1e3);
      RNAME = Form("Results_%s%.0f.root",REAC.c_str(),REXC*1e3);    
    }

    NBE = RSN - REXC; // neutron binding energy is separation energy minus excitation energy
    
  } else {
    
    DNAME = Form("AngDist_%s.txt",REAC.c_str());
    RNAME = Form("Results_%s.root",REAC.c_str());   
  
    if(verbose>1) printf(" ' %iSr %s ' Elastic Scattering, OM = ' %s ': -\n",A,reac.c_str(),om.c_str());  
  
    if(!fname.length())
      FNAME = Form("sr%i_%s_%s",A,REAC.c_str(),OM1.c_str());
    else SetFileName(fname);

  }
  
  // if everything worked well the set-up is complete  
  set_info = true;  

  if(!BuildReplacementStrings()){
    set_info = false;    
    return false;
  }
  
  // I don't like this at all.  
  DDIR=BASEDIR + Form("/AngularDistribution/sr%i",BEAMA);
    
  return true;
}

Bool_t TFrescoAnalysis::MakeFile(std::string file_opt){
  if(verbose>1)printf("\n\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    
  if(!CheckInfo())
    return false;
  
  if(file_opt.find("frin")==npos && file_opt.find("search")==npos && file_opt.find("min")==npos){
    printf("\n\n\t Error :  Option ' %s ' is invalid."
    "\n\t ->' frin ' creates a FRESCO input file"
    "\n\t ->' search ' creates an SFRESCO input file"
    "\n\t ->' min ' creates a script to run SFRESCO\n\n", file_opt.c_str());
    return false;
  }
  
  if(file_opt.find("frin")!=npos){
    const char *frin_temp = Form("fresco_%s_template.txt",TYPE==0?"dp":"elastic");
    if(!CreateFile(frin_temp,Form("%s.frin",FNAME.c_str())))
      return false;
  
    if(verbose>1)printf("\n -> Made ' %s.frin '\n",FNAME.c_str());
  }
  if(file_opt.find("search")!=npos && NFITPARS){
    
    if(!NFITDATA){
      printf("\n\n\t Error :  No data has been loaded and so SFRESCO cannot be used.\n\n");
      return false;
    }
    
    if(!CreateFile("sfresco_template.txt",Form("%s.search",FNAME.c_str())))
      return false;
  
    if(verbose>1)printf("\n -> Made ' %s.search '\n",FNAME.c_str());
  }
  if(file_opt.find("min")!=npos){
    if(!CreateFile("script_template.txt",Form("%s.min",FNAME.c_str())))
      return false;
  
    if(verbose>1)printf("\n -> Made ' %s.min '\n",FNAME.c_str());
  }
    
  if(verbose>1)printf("\n\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  
  return true;
}

TGraph *TFrescoAnalysis::CalculateAngDist(Double_t sf, UInt_t colour, Bool_t sine){

  TGraph *g;
  if(!MakeFile("frin"))
    return g; 
    
  system(Form("%s <  %s.frin  > %s.frout",FREX.c_str(),FNAME.c_str(),FNAME.c_str()));
  
  std::string foutname;
  if(TYPE==0)
    foutname.assign("fort.202");
  else  
   foutname.assign("fort.201");

  std::ifstream calcfile(foutname.c_str());
  if(!calcfile.is_open()){
    printf("\n\n\t Error :  Could not locate FRESCO angular distribution file ' %s '.\n\n",foutname.c_str());
    return g;
  }
  
  g = new TGraph();
  g->SetNameTitle(FNAME.c_str(),Form("%s; #theta_{CM} [#circ]; #frac{d#sigma}{d#Omega} [mb/sr]",FNAME.c_str()));
  g->SetMarkerColor(colour);
  g->SetLineColor(colour);
  
  std::string line;
  Double_t theta, sigma, d2r = TMath::DegToRad();
  Bool_t readvals = false;
  
  while (std::getline(calcfile, line)){
    if(line.find("#  Theta       sigma")!=npos)
      readvals = true;
    else if(readvals){
      std::istringstream ss(line);
      ss >> std::scientific >> theta >> sigma;
      if(theta && sigma){
        sigma*=sf;
        if(sine) sigma*=sin(theta*d2r);
        g->SetPoint(g->GetN(),theta,sigma);
      }
     }
  }
  
  return g;
}

Bool_t TFrescoAnalysis::ExtractFitResult(std::string line){

  if(line.find("ChiSq/N")!=npos){
    UInt_t chi2_beg = line.find("=");
    UInt_t chi2_end = line.find("from");
    std::string chi2val  = line.substr(chi2_beg+5,chi2_end-chi2_beg-5);
    //printf("\n chi2val = '%s'",chi2val.c_str());
    CHI2 = atof(chi2val.c_str());
    if(verbose>1) printf("\n\t -> Chi2 / N                 = %.3f\n",CHI2);
    return true;
  }
  
  if(line.find("value")==npos)
    return false;

  UInt_t vpos = line.find("value");
  UInt_t cpos = line.find(",");
  Double_t val = atof(line.substr(vpos+5,cpos-(vpos+5)).c_str());    
  
  UInt_t epos = line.find("error");
  Double_t err = atof(line.substr(epos+5,line.length()-(epos+5)).c_str());    
  
  if(line.find("Spec. Ampl")!=npos){ // spectroscopic amplitude
    SA = val;
    SA_ERR = err;
    if(verbose>1) printf("\n\t -> Spectroscopic Amplitude  = %.2e +/- %.1e  [rel. err. = %.3f %%]",SA,SA_ERR,SA_ERR/SA*100.0);    
  } else if(line.find("exptnorm")!=npos){ // normalization
    NORM = val;
    NORM_ERR = err;
    if(verbose>1) printf("\n\t -> Normalization            = %.2e +/- %.1e  [rel. err. = %.3f %%]",NORM,NORM_ERR,NORM_ERR/NORM*100.0);          
  } else if(line.find("=")!=npos){ // optical model parameter
    UInt_t spos = line.find("=");
    UInt_t kind, kbpot, kbtype, parnum;
    if(!GetParInfo(line.substr(spos+1,6),kind,kbpot,kbtype,parnum))
      return false;  
    std::string hname = Form("<KP%i_%iP%i>",kbpot,kbtype,parnum);
    Int_t indx = GetHolderIndex(hname.c_str());
  //  printf("\n**** SETTING HOLDERVAL %i ' %s ' to %.3f from %s ****\n",indx,hname.c_str(),val,HOLDERVAL.at(indx).c_str());  
    HOLDERVAL.at(indx) = Form("%.3f",val);
    
    if(verbose>1) printf("\n\t -> %6s                   = %.2e +/- %.1e  [rel. err. = %.3f %%]",line.substr(spos+1,6).c_str(),val,err,err/val*100.0);          
  }
  
  
  return true;
}

// Bool_t TFrescoAnalysis::SetGamEff(Double_t GamEff, Bool_t ScaleData){
// 
//   GAMEFF = GamEff;
//   if(NORM && SA && SA_ERR){
//     GAMEFF_ERR = GAMEFF*sqrt(pow(2*SA_ERR/SA,2)+pow(NORM_ERR,
//     SF_ERR = SF * sqrt(pow(2*SA_ERR/SA,2.0)+pow(GAMEFF_ERR/GAMEFF,2.0)+pow(NORM_ERR/NORM,2.0));
//   }
//  
//   if(ScaleData){
//     if(
//   }
// }

// Add the following functions which adjust the SF & SF_ERR using SA && SA_ERR

// Set Normalization
// Set TigEff

// reads root file
Bool_t TFrescoAnalysis::ExtractAngDistInfo(std::string dfile){

  if(verbose>1) printf("\n Searching ' %s ' for norm. and gam. eff. factors : \n",dfile.c_str());

  if(!dfile.length()) // look for default file name
    dfile.assign(DNAME);
  if(dfile.find("/")==npos) // use default path
    dfile = DDIR + "/" + dfile;
    
// use local copy of data ?
  TFile *f = new TFile(dfile.c_str(),"READ");
  if(!f->IsOpen()){
    printf("\n\t Error :  Couldn't find root file ' %s '\n\n",dfile.c_str());
    return false;
  }
  
  TH1D *hnorm = (TH1D*)f->Get("Results/Normalization");
  NORM = hnorm->GetBinContent(1);
  NORM_ERR = hnorm->GetBinError(1);  
  if(verbose>1) printf("\n\t -> Normalization            = %.2e +/- %.1e  [rel. err. = %.3f %%]",NORM,NORM_ERR,NORM_ERR/NORM*100.0);
  
  TH1D *hgam = (TH1D*)f->Get("Results/TigEfficiency");
  if(hgam){
    TCanvas *Options = (TCanvas*)f->Get("Results/Options");
    TPaveStats *ps = (TPaveStats*) Options->GetListOfPrimitives()->At(0);
    std::string str = ps->GetLineWith("Gam Energy")->GetTitle();
    
    GAM = atof(str.substr(str.find("=")+2,10).c_str());
    GAMEFF = 1/hgam->GetBinContent(1);
    GAMEFF_ERR = hgam->GetBinError(1)/hgam->GetBinContent(1)*GAMEFF;
    if(verbose>1) printf("\n\t -> TIG. Eff. @ %.1f keV   = %.2e +/- %.1e  [rel. err. = %.3f %%]",GAM,GAMEFF,GAMEFF_ERR,GAMEFF_ERR/GAMEFF*100.0);
  } else {    
    GAM = 0.0;
    GAMEFF = 1.0;
    GAMEFF_ERR = 0.0;
  } 
  if(verbose>1)printf("\n");
  f->Close();
  return true;
}

Bool_t TFrescoAnalysis::ExtractTotalCrossSection(const char *fort_name){

  if(verbose>1) printf("\n FRESCO output file ' %s ' will be searched for total cross section : -\n",fort_name);

  std::ifstream infile(fort_name);
  if(!infile.is_open()){
    printf("\n\t Error :  Couldn't find fort file ' %s '\n\n",fort_name);
    return false;
  }  
  
  
  std::string line;
  std::getline(infile,line);  //printf("\n%s",line.c_str());
  std::getline(infile,line);  //printf("\n%s",line.c_str());
  std::getline(infile,line);  //printf("\n%s",line.c_str());
  std::getline(infile,line);  //printf("\n%s",line.c_str());
  std::getline(infile,line);  //printf("\n%s",line.c_str());       

  Double_t v1,v2,v3,v4,v5,v6,v7;
  infile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7;
  if(v3!=REXC) {
    printf("\n\t Error :  Couldn't find total cross section in fort file ' %s '\n\n",fort_name);  
    return false;
  }
  CS = v7;
  CS_ERR = SF_ERR/SF*CS; // has same relative error as spec. factor
  if(verbose>1) printf("\n\t -> Total Cross Section      = %.2e +/- %.1e mb \n\n",CS,CS_ERR);
  
  return true;
}

TList *TFrescoAnalysis::FitData(Bool_t draw, Bool_t write){
  
  TList *list = new TList();
  if(!MakeFile("frin+search+min"))
    return list;
        
  system(Form("%s <  %s.min  > %s.min-out",SFREX.c_str(),FNAME.c_str(),FNAME.c_str()));
  
  if(verbose>1) printf("\n [ Wrote sfresco output to %s.min-out and %s-fit.plot ]\n",FNAME.c_str(),FNAME.c_str());
  
  std::string foutname = Form("%s-fit.plot",FNAME.c_str());   
  std::ifstream res_file(foutname.c_str());
  if(!res_file.is_open()){
    printf("\n\n\t Error :  Could not locate SFRESCO results file ' %s '.\n\n",foutname.c_str());
    return list;
  }  
  
  std::string line;
  Double_t theta, sigma, sigerr;
  Bool_t readvals = false;
 
   /////////////////  first look for experimental data  ///////////////// 
  while (std::getline(res_file, line)){
    if(line.find("#")!=npos) 
      ExtractFitResult(line); // SF is now calculated later in this function
    if(line.length() && line.find("#")==npos && line.find("@")==npos)
        break; // skip to experimental data  
  }

  gdata = new TGraphErrors();
  gdata->SetName(Form("%s_ExpData",FNAME.c_str()));
  std::string units, info;
  if(TYPE==0){
    units = "mb/sr";
    info = Form("SF = %.2e +/- %.1e",SF,SF_ERR);
  } else {
    units = "Ratio to Rutherford";
    info = Form("Norm = %.2e +/- %.1e",NORM,NORM_ERR);
  }
  
  gdata->SetTitle(Form("^{%i}Sr(%c,%c) @ %.3f MeV/u; #theta_{CM} [#circ]; #frac{d#sigma}{d#Omega} [%s]",
                  BEAMA,REAC.at(0),REAC.at(1),EPERU,units.c_str()));
  gdata->SetMarkerColor(1);
  gdata->SetLineColor(1);
  gdata->SetMarkerStyle(8);
  gdata->SetMarkerSize(0.7);  
  
  UInt_t n=0;  
  do{
    std::istringstream ss(line);
    ss >> std::scientific >> theta >> sigma >> sigerr;
    if(theta && sigma && sigerr){
      gdata->SetPoint(n,theta,sigma);
      gdata->SetPointError(n,0,sigerr);
    }
    n++;
  } while (std::getline(res_file, line) && line.find("&")==npos);

  Double_t thterr = 10.0, thterr_tmp;
  for(int i=1; i<n; i++){ // get smallest x distance between points
    thterr_tmp = 0.5*(gdata->GetX()[i]-gdata->GetX()[i-1]);
    if(thterr_tmp<thterr)
      thterr = thterr_tmp;
  }
  for(int i=0; i<n; i++) // set x error 
    gdata->SetPointError(i,thterr,gdata->GetEY()[i]);
    
   /////////////////  Now look for theory data  /////////////////
  
  while (std::getline(res_file, line)){
    if(line.length() && line.find("#")==npos && line.find("@")==npos)
        break; // skip to theory data  
  }  

  gcalc = new TGraph();
  gcalc->SetName(Form("%s_Thry",FNAME.c_str()));
  gcalc->SetTitle(Form("%s [DWBA]; #theta_{CM} [#circ]; #frac{d#sigma}{d#Omega} [%s]",
                  gdata->GetTitle(),units.c_str()));
  gcalc->SetMarkerSize(5);  
  
  int col[] = {2, kMagenta, 8, kCyan, 4};
  gcalc->SetMarkerColor(col[L]);
  gcalc->SetLineColor(col[L]);  

  theta = 0.0;
  sigma = 0.0;
  do{
    std::istringstream ss(line);
    ss >> std::scientific >> theta >> sigma;
    if(theta && sigma)
      gcalc->SetPoint(gcalc->GetN(),theta,sigma);
  } while (std::getline(res_file, line) && line.find("&")==npos);
     
  // calculate SF and SF_ERR
  if(TYPE==0 && SA && SA_ERR){
    // default root file name is the same as AngDist name but with .root
    int pos1 = DNAME.find("dp"), pos2 = DNAME.find(".");    
    std::string rfile = Form("Results_%s.root",DNAME.substr(pos1,pos2-pos1).c_str());      
    // get NORM and GAM from root file    
    ExtractAngDistInfo(rfile);  
    // Combine this uncertainty with NORM, GAM & BR
    CalculateSF(true);
  }
       
  if(draw){
    TCanvas *c = new TCanvas("FitData","FitData",800,600);  
   // TLegend *leg = new TLegend(0.6,0.7,0.95,0.9); 
     if(TYPE==0)
       info = Form("SF = %.2e +/- %.1e",SF,SF_ERR);
 
    TLegend *leg = new TLegend(0.6,0.65,0.96,0.9,Form("#chi^{2}/N = %.2f, %s",CHI2,info.c_str()));
    
    leg->AddEntry(gdata,"Experimental Data","lp");
    leg->AddEntry(gcalc,"FRESCO Fit","l");
    
    // drawing thry first sets range to be 0-180    
    gcalc->GetXaxis()->SetRangeUser(0,180);
    gcalc->Draw("A L");
    gdata->Draw("same P");    
    
    leg->Draw();
    list->Add(c);
    c->Print(Form("FitData_%iSr%s.pdf",BEAMA,REAC.c_str()));   
  }  
  list->Add(gcalc);
  list->Add(gdata); 
  

    
  // this should also add in the errors on tigress and normalization if transfer
  if(TYPE==0 && !ExtractTotalCrossSection("fort.13"))
    return list;   
  
  if(write){
    TFile *fout = new TFile(Form("FitData_%iSr%s_%s.root",BEAMA,REAC.c_str(),OM.c_str()),"RECREATE");
    if(verbose) printf("\n\n   Saved results to ' %s '.\n\n",fout->GetName());  
    list->Write();
  }

  return list;
}

Bool_t TFrescoAnalysis::CalculateSF(Bool_t usechi2, Double_t br, Double_t br_err){

  if(!SA || !SA_ERR){
    printf("\n\n Error :  Spectroscopic Amplitude has not been set.\n\n");
    return false;
  } 
 
  if(usechi2){    // use chi2 fit curve to estimate SF_RELERR
    Double_t SF_RELERR = ChiSquaredError();  
    // Divide by 2 to get SA_RELERR and convert to SA_ERR  
    SA_ERR = 0.5*SF_RELERR*SA;    
    if(verbose>1) printf("\n\t -> Spectroscopic Amplitude  = %.2e +/- %.1e  [rel. err. = %.3f %%]",SA,SA_ERR,SA_ERR/SA*100.0);    
  }
  
  SF = SA * SA;
  Double_t SF_RELERR2 = pow(2*SA_ERR/SA,2.0);
  
  if(NORM && NORM_ERR) // use normalization constant in total SF error
    SF_RELERR2 += pow(NORM_ERR/NORM,2.0);
  else if(verbose>1) printf("\n\n Warning :  Normalization has not been set.\n\n");

  if(GAMEFF && GAMEFF_ERR) // add tigress efficiency uncertainty in quadrature
    SF_RELERR2 += pow(GAMEFF_ERR/GAMEFF,2.0);
  if(br && br_err)  // add branching ratio uncertainty in quadrature
    SF_RELERR2 += pow(br_err/br,2.0);
  
  SF_ERR = SF * sqrt(SF_RELERR2);

  if(verbose>1) printf("\n\t -> Spectroscopic Factor     = %.2e +/- %.1e  [rel. err. = %.3f %%]\n\n",SF,SF_ERR,sqrt(SF_RELERR2)*100.0); 
  
  return true;     
}


Double_t TFrescoAnalysis::GetPotVal(std::string vname){

  UInt_t kind, kbpot, kbtype, parnum;
  Double_t val=0.00;
  // get 'real' name
  if(!GetParInfo(vname,kind,kbpot,kbtype,parnum))
    return val;

  std::string parname = GetParName(kbpot,kbtype,parnum);
  if(!parname.length())
    return val;
    
  // get index and value  
  Int_t indx = GetHolderIndex(parname);
  if(indx>=0){
    val = (double)atof(HOLDERVAL.at(indx).c_str());
  //  printf("\n %s = %s [%f]\n\n",vname.c_str(),HOLDERVAL.at(indx).c_str(),val);    
  }else if(verbose) 
    printf("\n\n Error :  Parameter ' %s ' was not recognized",vname.c_str());

  return val;
}

TCanvas *TFrescoAnalysis::UserScanPotVal(UInt_t kbpot, UInt_t kbtype, UInt_t parnum, Double_t valmin, Double_t valmax, UInt_t nsteps){

  std::string parname = GetParName(kbpot,kbtype,parnum,false);

  GetRid(parname.c_str());
  TCanvas *c = new TCanvas(parname.c_str(),parname.c_str(),800,600);

  Bool_t vtmp = verbose;
  verbose = 1;
  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  Double_t val;
  for(int i=0; i<nsteps; i++){
    
    val = valmin+(valmax-valmin)*(double)i/(double)(nsteps-1);
    
    UserSetPotVal(kbpot,kbtype,parnum,val);
    
    TGraph *g = CalculateAngDist(1);
    if(!g)
      continue;
    g->SetName(Form("PotScan_%s%.0f",parname.c_str(),val));
    g->SetLineColor(kRed+i);

    if(i)
      g->Draw("same L");
    else{
      g->SetTitle(Form("^{%i}Sr(%c,%c) : '%s' Parameter",BEAMA,REAC.at(0),REAC.at(1),parname.c_str()));
      g->GetXaxis()->SetTitle("#theta_{CM} [#circ]");
      g->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} [mb/sr]");      
      g->Draw("A L");     
    }
    leg->AddEntry(g,Form("%s = %.2f",parname.c_str(),val),"l");
  }
  leg->Draw();
  
  verbose = vtmp;
  
  return c;
}

TCanvas *TFrescoAnalysis::UserScanPotVal(std::string parname, Double_t valmin, Double_t valmax, UInt_t nsteps){

  UInt_t kind, kbpot, kbtype, parnum;
  if(!GetParInfo(parname,kind,kbpot,kbtype,parnum))
    return 0;

  return UserScanPotVal(kbpot,kbtype,parnum,valmin,valmax,nsteps);
}

Bool_t TFrescoAnalysis::UserScanPotVal2D(std::string var1, Double_t val1lo, Double_t val1hi, UInt_t ns1, 
                                         std::string var2, Double_t val2lo, Double_t val2hi, UInt_t ns2){

  if(ns1*ns2>100 || ns1*ns2==0){
    printf("\n\t Error :  ns1*ns2>100 && 0 Not allowed [val=%i]\n\n",ns1*ns2);
    return false;
  }
  TGraphErrors *ge[ns1*ns2];
  TGraph *gt[ns1*ns2];
  
  std::string fname = Form("PotScan_%s_%s",var1.c_str(),var2.c_str());
  SetFileName(fname);
  
   Bool_t vtmp = verbose;
   verbose = 1;  
  
  if(!UserSetPotVal(var1,val1lo)){ 
    printf("\n\nt Error :  Pot. val. '%s' not recognized.\n\n",var1.c_str());
    return false;
  }
  if(!UserSetPotVal(var2,val2lo)){
    printf("\n\nt Error :  Pot. val. '%s' not recognized.\n\n",var2.c_str());
    return false;
  }          
  
  TH2F *hchi2 = new TH2F("Chi2",Form("#chi^{2} Surface; %s; %s",var1.c_str(),var2.c_str()),
            ns1,val1lo,val1hi,ns2,val2lo,val2hi);

  UInt_t indx=0,col=0;
  Double_t val1=0, val2=0, xlo=0, xhi=0, ylo=0, yhi=0;
  Double_t chi2min=1e3, val1min=0, val2min=0;
  std::string info, infomin;
  
  TList *list;
  TCanvas *c1 = new TCanvas("DWBA","DWBA Fits",1200,900);
  TPad *pad[ns1*ns2];
  
  gStyle->SetTitleW(0.88);
  gStyle->SetTitleH(0.12);
  gStyle->SetTitleTextColor(kBlue);
  
  c1->Divide(ns1,ns2,0,0);
    
  for(int i=0; i<(int)ns1; i++){
    val1 = val1lo+i*(val1hi-val1lo)/(double)(ns1-1);  
    UserSetPotVal(var1.c_str(),val1);
    
    for(int j=0; j<(int)ns2; j++){ 
      val2 = val2lo+j*(val2hi-val2lo)/(double)(ns2-1);
      UserSetPotVal(var2.c_str(),val2);
      
      indx = i*ns2+j;

      list = FitData(false,false);
      if(list->GetEntries()!=2)
        return false;
        
      if(TYPE==0)
        info.assign(Form("SF = %.2e+/-%.0e",SF,SF_ERR));
      else
        info.assign(Form("Norm = %.2e+/-%.0e",NORM,NORM_ERR));
            
      ge[indx] = (TGraphErrors*)list->FindObject(Form("%s_ExpData",FNAME.c_str()));
      ge[indx]->SetMarkerStyle(1);
      ge[indx]->SetMarkerSize(0.1);
      ge[indx]->SetTitle(Form("%s = %.2f, %s = %.2f, %s",var1.c_str(),val1,var2.c_str(),val2,info.c_str()));      
      gt[indx] = (TGraph*)list->FindObject(Form("%s_Thry",FNAME.c_str()));
      if(CHI2<2.5)
        col=8;
      else if(CHI2<5)
        col=kOrange;
      else if(CHI2>5)
        col=kRed;
      gt[indx]->SetLineColor(col);
                      
      c1->cd(indx+1);    
      ge[indx]->Draw("AP");
      gt[indx]->Draw("same C");
      
      TLegend *leg = new TLegend(0.6,0.55,0.96,0.8,Form("#chi^{2}/N = %.2f",CHI2));
//      leg->SetLineColor(col);
      leg->AddEntry(ge[indx],Form("Exp. Data"),"lp");
      leg->AddEntry(gt[indx],"DWBA","l");      
      leg->Draw("same");
      
      hchi2->Fill(val1,val2,CHI2);
      if(CHI2<chi2min){
        chi2min=CHI2;
        val1min=val1;
        val2min=val2;
        infomin.assign(info);
      } 
      printf("\n-> %s = %.2f \t %s = %.2f \t chi2 = %.2f \t %s",var1.c_str(),val1,var2.c_str(),val2,CHI2,info.c_str());
 
    }
  }
  printf("\n\n Complete \n\n\n\t < - - - - - - - - - > BEST KIND < - - - - - - - - >");
  
  printf("\n\n %s = %.2f \t %s = %.2f \t chi2 = %.2f \t %s\n\n",
  var1.c_str(),val1min,var2.c_str(),val2min,chi2min,infomin.c_str());
  
  verbose = vtmp;
  return true;
}


Bool_t TFrescoAnalysis::DoFrescoAnalysis(){

  if(verbose>1)printf("\n\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  if(!CheckInfo())
    return false;
  
  std::string new_dir;
  if(REAC.find("dp")!=npos){
    new_dir = Form("%s/State@%.0fkeV",FDIR.c_str(),REXC*1e3);  
    if(EXCHI)
      new_dir+="/ExcHi";
  } else 
    new_dir = Form("%s/%s",FDIR.c_str(),REAC.c_str()); 
    
  //  copy angdist .txt data and result .root file into this directory too
  system(Form("cp %s/%s %s",DDIR.c_str(),DNAME.c_str(),new_dir.c_str()));
  system(Form("cp %s/%s %s",DDIR.c_str(),RNAME.c_str(),new_dir.c_str()));
  
  if(verbose>1)printf("\n\n\t____ PREPARING INPUT FILES _____\n");

  if(!MakeFile("frin+search+min"))
    return false;
    
  if(verbose>1)printf("\n\n\n\t____ PARSING RESULTS FILES _____\n");
    
  // if (d,p), parse various output files to extract SF, Sig, and make final results
  if(TYPE==0 && !ExtractAngDistInfo(RNAME.c_str()))
    return false;
  
  if(TYPE==0 && !ExtractTotalCrossSection("fort.13"))
    return false;
  
  const char *rootfile_name = Form("%s/%s%s.root",new_dir.c_str(),FNAME.c_str(),EXCHI?Form("_ExcHi%.01f",EXCHI):"");    
  if(!SaveResults(rootfile_name,new_dir))
    return false;
   
  if(verbose>1)printf("\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
   
  return true;
}

Bool_t TFrescoAnalysis::SaveResults(std::string rootfile_name, std::string dir){

  
  const char * rfile_name = Form("%s/%s-fit.plot",dir.c_str(),FNAME.c_str());
  std::ifstream res_file(rfile_name);  
  std::string line;
    
  if(!res_file.is_open()){
    printf("\n\t Error :  Couldn't find the file '%s'.\n",rfile_name);
    return false;
  }
  
   /////////////////  first look for experimental data  /////////////////

  while (std::getline(res_file, line)) 
    if(line.length() && line.find("#")==npos && line.find("@")==npos)
      break; // skip to experimental data
  
  Int_t n=0;
  TGraphErrors *gdata = new TGraphErrors();
  gdata->SetNameTitle(Form("%s_Data",FNAME.c_str()),Form("sr95_%s%.0f_Data; #theta_{CM} [#circ]; #frac{d #sigma}{d #Omega} [#frac{mb}{sr}]",REAC.c_str(),REXC));
  gdata->SetMarkerSize(5);
     
  Double_t theta, thterr, sigma, sigerr;
  std::string theta_str, thterr_str, sigma_str, sigerr_str;
  Int_t pos1,pos2,pos3;
  
  if(verbose==3) printf("\n\tReading anglar distribution data :-\n");   
  do{
    if(line.find("&")!=npos) // end of data
      break;
   
    if(verbose==3) printf("\n%s",line.c_str());
   // manual data extraction because of mixed numerical formats
    pos1 = line.find(".");
    pos2 = line.find(".",pos1+1);    
    pos3 = line.find(".",pos2+1);  
    
    theta_str=line.substr(pos1-3,pos2-pos1);
    sigma_str=line.substr(pos2-2,pos3-pos2); 
    sigerr_str=line.substr(pos3-2,12);

    // make expt. data graph      
    theta  = atof(theta_str.c_str());
    sigma  = atof(sigma_str.c_str());
    sigerr = atof(sigerr_str.c_str());   
         
    gdata->SetPoint(n,theta,sigma);  
    gdata->SetPointError(n,0,sigerr);  

    n++;    
  } while(std::getline(res_file, line));
      
  if(n<5){ // check expt. dat sie
    printf("\n\t Error :  Found insufficient experimental data [ %i points ] in the file '%s'.\n",n,rfile_name);
    return false;
  }
  thterr = 0.5*(gdata->GetX()[1]-gdata->GetX()[0]);
  for(int i=0; i<n; i++) // set x error 
    gdata->SetPointError(i,thterr,gdata->GetEY()[i]);


   /////////////////  Now look for theory data  /////////////////
  
  while (std::getline(res_file, line)) 
    if(line.length() && line.find("#")==npos && line.find("@")==npos)
      break; // skip to experimental data
      
  n=0;    
  TGraphErrors *gthry = new TGraphErrors();
  gthry->SetNameTitle(Form("%s_Thry",FNAME.c_str()),Form("%s_Thry; #theta_{CM} [#circ]; #frac{d #sigma}{d #Omega} [#frac{mb}{sr}]",FNAME.c_str()));
  
  gthry->SetMarkerSize(5);  
  gthry->SetMarkerColor(0.5*(double)L+2);
  gthry->SetLineColor(0.5*(double)L+2);   
       
  if(verbose==3) printf("\n\n\tReading anglar distribution thry :-\n");   
       
  do{
    if(line.find("&")!=npos) // end of data
      break;

 //   if(verbose==3) printf("\n%s",line.c_str());
    
    pos1 = line.find(".");
    pos2 = line.find(".",pos1+1);    
    
    theta_str=line.substr(pos1-3,pos2-pos1);
    sigma_str=line.substr(pos2-2,12); 

    theta  = atof(theta_str.c_str());
    sigma  = atof(sigma_str.c_str());

   // make thry graph      
    gthry->SetPoint(n,theta,sigma);  
    n++;    
  } while(std::getline(res_file, line));

  if(n<5){ // check expt. dat sie
    printf("\n\t Error :  Found insufficient theoretical data [ %i points ] in the file '%s'.\n",n,rfile_name);
    return false;
  }
   
   /////////////////  Now draw results and make residuals  /////////////////

  GetRid("Results");  
  TCanvas *c = new TCanvas("Results","Results",1000,1000);
  
  TPad *p1 = new TPad("p1","p1",0,0.5,0.5,1);
  p1->Draw(); 
  
  TPad *p2 = new TPad("p2","p2",0.5,0.5,1,1);
  p2->Draw();  
 
  TPad *p3_1 = new TPad("p3_1","p3_1",0,0,0.5,0.25);
  p3_1->SetTopMargin(0);
  p3_1->Draw();  

  TPad *p3_2 = new TPad("p3_2","p3_2",0,0.25,0.5,0.5);
  p3_2->SetBottomMargin(0);  
  p3_2->Draw();   
 
  TPad *p4 = new TPad("p4","p4",0.5,0,1,0.5);
  p4->Draw();     
  
  
  TLegend *leg = new TLegend(0.5,0.65,0.95,0.9);
  leg->AddEntry(gdata,"Experimental Data","lp"); 
  if(REAC.find("dp")!=npos)
    leg->AddEntry(gthry,Form("SF = %.3f[%.3f], Chi2 = %.3f",SF,SF_ERR,CHI2),"lp");   
  else 
    leg->AddEntry(gthry,Form("NORM = %.3e[%.3e], Chi2 = %.3f",NORM,NORM_ERR,CHI2),"lp");   

  p1->cd();
  gdata->Draw("AP");    
  gthry->Draw("C same");
  leg->Draw();  
  gPad->SetLogy();


  TGraphErrors *gres = new TGraphErrors();
  gres->SetNameTitle("Residuals","Residuals From FRESCO Fit; #theta_{CM} [#circ]; #frac{d #sigma}{d #Omega} Residual [#frac{mb}{sr}]");
  TGraphErrors *gres_sum = new TGraphErrors();
  gres_sum->SetNameTitle("ResidualsSum","Cumulative Residuals From FRESCO Fit; #theta_{CM} [#circ]; #frac{d #sigma}{d #Omega} Total Residual [#frac{mb}{sr}]");
  
  TGraphErrors *gchi2 = new TGraphErrors();
  gchi2->SetNameTitle("ResidualChi2","#chi^{2} Per Datum From FRESCO Fit; #theta_{CM} [#circ]; #chi^{2}");
  TGraphErrors *gchi2_sum = new TGraphErrors();
  gchi2_sum->SetNameTitle("ResidualChi2Sum","Cumulative #chi^{2} Per Datum From FRESCO Fit; #theta_{CM} [#circ];Total #chi^{2}");
  Double_t xdat, yres, xerr, yerr, ysum=0, yesum=0, chi2=0, chi2sum=0;    
  //(TGraphErrors*) gdata->Clone(Form("%sResiduals",gdata->GetName()))
  for(int i=0; i<gdata->GetN(); i++){
  
    xdat = gdata->GetX()[i];
    yres = gdata->GetY()[i]-gthry->Eval(gdata->GetX()[i]);
    xerr = gdata->GetEX()[i];
    yerr = gdata->GetEY()[i];    
    
    gres->SetPoint(i,xdat,yres);
    gres->SetPointError(i,xerr,yerr);
    
    ysum += yres;
    yesum += ysum;
    
    gres_sum->SetPoint(i,xdat,ysum);
   // gres_sum->SetPointError(i,0,yesum);  
    
    chi2 = pow(yres/yerr,2.0);
    gchi2->SetPoint(i,xdat,chi2);

    chi2sum += chi2/(double)(i+1);    
    gchi2_sum->SetPoint(i,xdat,chi2sum);    
  }
  
  p3_2->cd();
  gres->Draw("AP");
  gres->GetYaxis()->SetTitleSize(0.08);
  gres->GetYaxis()->SetTitleOffset(0.8); 
  gres->GetYaxis()->CenterTitle();  
  p3_2->SetGridy();

  p3_1->cd();
  gchi2->SetLineColor(kRed);
  gchi2->SetMarkerColor(kRed);
  gchi2->SetMarkerStyle(27);
  gchi2->GetYaxis()->SetTitleSize(0.08);
  gchi2->GetYaxis()->SetTitleOffset(0.8);
  gchi2->GetYaxis()->CenterTitle();
  gchi2->GetXaxis()->SetTitleSize(0.08);
  gchi2->GetXaxis()->SetTitleOffset(0.6);
  gchi2->GetXaxis()->CenterTitle();
  
  gchi2->Draw("ALP");
 // gchi2_sum->Draw("same LP");
  
  
  TList *potlist = new TList();
    
  // potential
  for(int j=0; j<2; j++){
    if(j==0)
      p2->cd();
    else 
      p4->cd();
        
    TList *list = MakePotentials(j+1);
    TIter iter(list);
    
    TGraph *g;
  
    Double_t ymin_tmp, ymax_tmp, ymin=100.0, ymax=-100.0;
    TLegend *leg1 = new TLegend(0.5,0.6,0.95,0.9);

    for(int i=0; i<2; i++){
      n=0;
      iter.Reset();
      while((g = (TGraph*)iter.Next())){  
        if(i==0){
          ymin_tmp = TMath::MinElement(g->GetN(),g->GetY());
          ymax_tmp = TMath::MaxElement(g->GetN(),g->GetY());
    
          if(ymin_tmp<ymin)
            ymin = ymin_tmp;
          if(ymax_tmp>ymax)
            ymax = ymax_tmp;
        
        } else {
          if(fabs(ymin)<fabs(ymax))
            ymin = -ymax;
          if(fabs(ymin)>fabs(ymax))
            ymax = -ymin;      
              
          if(verbose>2){        
            printf("\n\n--------------------------------------------------\n");
            printf("\n Drawing %s",g->GetName());
            printf("\n\t Setting y axis to be [%.3f - %.3f]\n",ymin,ymax);
            g->Print();
            printf("\n--------------------------------------------------\n");
          }
          g->GetYaxis()->SetRangeUser(ymin,ymax);
          leg1->AddEntry(g,g->GetTitle(),"l");
          potlist->Add(g);
        
          if(n==0)
            g->Draw("AC");
          else
            g->Draw("same C");
        }
        n++;
      }
    }
    leg1->Draw();
  }
  
  
  TFile *file = new TFile(rootfile_name.c_str(),"RECREATE");
  
  c->SaveAs(Form("%s_FrescoResults.pdf",FNAME.c_str()));
 
 // c->Write();
  gdata->Write();
  gthry->Write();
  gres->Write();
  gchi2->Write();
  
  file->mkdir("Potentials");
  file->cd("Potentials");
  potlist->Write();
  
  file->Close();
  
  return true;  
}

Bool_t TFrescoAnalysis::SetParticles(UInt_t A, Double_t exc, std::string reac){

  if(A<94 || A>96){
    printf("\n\t Error :  A=%i is not valid. A must be either 94, 95 or 96.\n\n",A);
    return false;
  }
  
  if(reac.find("dp")!=npos)
    TYPE = 0;
  else if(reac.find("pp")!=npos)
    TYPE = 1;    
  else if(reac.find("dd")!=npos)
    TYPE = 2;    
  else if(reac.find("cc")!=npos){
    printf("\n\t Error :  Carbon elastic scattering is not currently implemented.\n\n");
    return false;  
  } else if(reac.find("dt")!=npos)
    TYPE = 4;                
  else {
    printf("\n\t Error :  reaction ' %s ' not recognized.\n\n",reac.c_str());
    return false;
  }  

  BEAMA = A;
  REAC.assign(reac);
  
  static double eperu[]    = {5.341, 5.378, 5.35}; // 94Sr, 95Sr, 96sr
  static double elab[]     = {499.5, 510.9, 513.5}; // 94Sr, 95Sr, 96sr
  static double mass[]     = {92.9140, 93.9154, 94.9194, 95.9217, 96.9262};    // 94Sr, 94Sr, 95Sr, 96Sr, 97Sr
  static double spin[]     = {2.5, 0.0, 0.5, 0.0, 0.5}; // 93Sr, 94Sr, 95Sr, 96Sr, 97Sr  
  static double qval_dp[]  = {2.123, 3.654, 1.499}; // 94Sr->95Sr, 95Sr->96Sr, 96Sr->97Sr  
  static double qval_dt[]  = {-0.5739, 1.9094, 0.3783}; // 94Sr->93Sr, 95Sr->94Sr, 96Sr->95Sr  
  static double sn[]       = {5.290, 6.831, 4.348, 5.879, 3.724}; // 93Sr, 94Sr, 95Sr, 96Sr, 97Sr 
    
  int indx = A-94;
  EPERU = eperu[indx];
  ELAB  = elab[indx];
  BMASS = mass[indx+1];
  BSPIN = spin[indx+1];  
  
  if(TYPE==0){ // dp
    RECOA = A+1;
    RMASS = mass[indx+2];
    QVAL  = qval_dp[indx] - exc;
    RSN   = sn[indx+2];
    TARGA = 2;    
    TMASS = 2.0141;
    TSPIN = 1.0;
    NKBPOT = 5; // transfer requires 5 potentials    
  } else if(TYPE==4){ // dt
    RECOA = A-1;
    RMASS = mass[indx];
    QVAL  = qval_dt[indx] - exc;
    RSN    = sn[indx];    
  } else { // elastic
    RECOA = A;
    RMASS = mass[indx];
    QVAL  = 0.0;
    RSN    = sn[indx]; 
    NKBPOT = 1;  // only one potential required      
    if(TYPE==1){
      TARGA = 1;
      TMASS = 1.008;
      TSPIN = 0.5;
    }else if(TYPE==2){
      TARGA = 2;    
      TMASS = 2.0141;
      TSPIN = 1.0;
    }
  }
  
  return true;
}
  
Bool_t TFrescoAnalysis::SetOM(std::string reac, std::string om){

  if(om.find(".")!=npos && ReadPotVals(om)>15){
    OM.assign(om.substr(0,om.find(".")));
    OM1.assign(OM);
    OM2.assign(OM);        
    return true;
  }
  
  if(om=="PP")
    om.assign("PP1+PP3");
  
  if(reac[0]=='d'){
    if(om.find("PP1")!=npos){
      OM1.assign("PP1");
    } else if(om.find("PP2")!=npos){
      OM1.assign("PP2"); 
    } else if(om.find("LH")!=npos){
      OM1.assign("LH");  
    } else if(om.find("D")!=npos){
      OM1.assign("D");        
    } else{
      printf("\n Error :  Optical model option ' %s ' does not contain a deuteron potential.",om.c_str());
      printf("\n\t  For deuteron, try : 'PP1' or 'PP2' or 'LH' or 'D'");
      return false;  
    }
  } else if(reac[0]=='p'){
      if(om.find("PP3")!=npos){
      OM1.assign("PP3");
    } else if(om.find("PP4")!=npos){
      OM1.assign("PP4"); 
    } else if(om.find("BG")!=npos){
      OM1.assign("BG");  
    } else if(om.find("CH")!=npos){
      OM1.assign("CH");        
    } else{
      printf("\n Error :  Optical model option ' %s ' does not contain a deuteron potential.",om.c_str());
      printf("\n\t  For proton, try : 'PP3' or 'PP4' or 'BG' or 'CH'");
      return false;  
    }
  }
  
  if(reac[1]=='p'){
    if(om.find("PP3")!=npos){
      OM2.assign("PP3");
    } else if(om.find("PP4")!=npos){
      OM2.assign("PP4"); 
    } else if(om.find("BG")!=npos){
      OM2.assign("BG");  
    } else if(om.find("CH")!=npos){
      OM2.assign("CH");        
    } else{
      printf("\n Error :  Optical model option ' %s ' does not contain a deuteron potential.",om.c_str());
      printf("\n\t  For deuteron, try : 'PP1' or 'PP2' or 'LH' or 'D'");
      return false;  
    }
  }    
  
  OM.assign(om);
  
  UInt_t nvals = SetPotVals();
  if(!nvals){
    printf("\n Error : Potential parameters not set [nvals=%i].\n\n",nvals); 
    OM1.clear();
    OM2.clear();
    OM.clear();     
    return false;
  }    
  
  return true;
}

Bool_t TFrescoAnalysis::SetTransfer(Double_t exc, UInt_t l, Double_t jo, Double_t jf){

  if(l!=0 && l!=2 && l!=4){
    printf("\n\t Error :  L=%i is not valid. L must be either 0, 2 or 4.\n\n",l);
    return false;
  }

  if(exc>6.0 || exc<0){
    printf("\n\t Error :  Exc=%.0f is not valid.. Use units of MeV! \n\n",exc);
    return false;  
  }

  REXC = exc;
  L = l;
  
  if(TYPE==0){
    // for 95Sr(d,p) the total spin is the coupling between the 95Sr nucleus and the transferred neutron  
    if(jo>0) JO = jo;
    if(jf>=0) JF = jf;

    if(BEAMA==95){
      if(l==0){ // s1/2
        NO = 3;
        if(jo<=0) JO = 0.5;
        if(jf<0) JF = 0;    
      } else if (l==2 && (jo==1.5 || jo<0)){ // d3/2
        NO = 2;
        if(jo<=0) JO = 1.5;
        if(jf<0) JF = 2;       
      }  else if (l==2 && jo==2.5){ // d5/2
        NO = 2;
        if(jo<=0) JO = 2.5;
        if(jf<0) JF = 2;       
      } else if (l==4){ // g7/2
        NO = 1;
        if(jo<=0) JO = 3.5;
        if(jf<0) JF = 4;        
      } else {
        printf("\n\t Error :  Could not recognize orbital .. L=%i, jo=%.1f, jf=%.1f\n\n",l,jo,jf);
        return false;
      }
    } 
      // for 94Sr(d,p) and 96Sr(d,p) the total spin is just the spin of the transferred neutron
    else if(BEAMA==94 || BEAMA==96){
      if(l==0){ // s1/2
        NO = 3;
        JO = 0.5;
        JF = 0.5;    
      } else if (l==2 && (jo==1.5 || jo<=0)){ // d3/2
        NO = 2;
        JO = 1.5;
        JF = 1.5;       
      }  else if (l==2 && jo==2.5){ // d5/2
        NO = 2;
        JO = 2.5;
        JF = 2.5;       
      } else if (l==4){ // g7/2
        NO = 1;
        JO = 3.5;
        JF = 3.5;    
      }
    }  
  } else if (TYPE==4){
    // this is complicated. (d,t) is a neutron pickup reaction from the Sr beam.
    // A neutron is more likely to be found in a full orbitals close to the Fermi surface 
    // than mostly empty orbitals in the valence space.
    // The lower orbitals would also be more bound though.
    printf("\n\n\t Error :  (d,t) is under construction ...\n\n");
    return false;    
  }
  
  return true;
}

UInt_t TFrescoAnalysis::GetAllowedTransfers(UInt_t A, std::vector<int> &lo, std::vector<double> &jo, std::vector<double> &jf, Bool_t all){

  lo.clear();
  jo.clear();
  jf.clear();
  
  if(A==94 || A==96){
    lo.push_back(0); jo.push_back(0.5); jf.push_back(0.5);
    lo.push_back(2); jo.push_back(1.5); jf.push_back(1.5);
    lo.push_back(2); jo.push_back(2.5); jf.push_back(2.5);
    lo.push_back(4); jo.push_back(3.5); jf.push_back(3.5);         
  } else if(A==95){
    lo.push_back(0); jo.push_back(0.5); jf.push_back(0.0);
    if(all){lo.push_back(2); jo.push_back(1.5); jf.push_back(1.0);}
    lo.push_back(2); jo.push_back(1.5); jf.push_back(2.0);
    if(all){lo.push_back(2); jo.push_back(2.5); jf.push_back(2.0);}
    lo.push_back(2); jo.push_back(2.5); jf.push_back(3.0);
    if(all){lo.push_back(4); jo.push_back(3.5); jf.push_back(3.0);}      
    lo.push_back(4); jo.push_back(3.5); jf.push_back(4.0);         
  } else 
    printf("\n\n\t Error :  A = %i is invalid. A must be 94, 95 or 96.\n\n",A);
  
  return (UInt_t)lo.size();
}

UInt_t TFrescoAnalysis::ReadPotVals(std::string fname){

  std::ifstream infom(Form("%s/%s",MDIR.c_str(),fname.c_str()));
  Int_t kpot,kbtype,parnum;
  Double_t val;
  UInt_t n = 0, ndel = 0;
  std::string name,line;
  Int_t indx;

  while(infom.good()){
    
    infom >> kpot >> kbtype >> parnum >> val;
    
    name = Form("<KP%i_%iP%i>",kpot,kbtype,parnum);
 /*
    if(kpot>1 && (TYPE==1 || TYPE==2))
      break; // don't store potentials beyond first one if elastic
  */  
    indx = GetHolderIndex(name); 
    if(indx>=0){
      if(verbose==3) printf("\n\n - Deleted holderval [%i] %s : Value was %s",
        indx,HOLDERNAME.at(indx).c_str(),HOLDERVAL.at(indx).c_str());
      HOLDERNAME.erase(HOLDERNAME.begin()+indx);
      HOLDERVAL.erase(HOLDERVAL.begin()+indx); 
      ndel++;
    }
     
    HOLDERNAME.push_back(name.c_str());
    HOLDERVAL.push_back(Form("%.3f",val));
    if(verbose==3) printf("\n\n + Added holderval ' %s ' : Value = %.3f",name.c_str(),val);     
    n++;
  }    
  
  if(verbose>1 && ndel) printf("\n\n -- Deleted %i potential parameters",ndel);
  if(verbose>1) printf("\n\n ++ Added %i potential parameters from ' %s ' in ' %s '\n\n",n,fname.c_str(),MDIR.c_str());

  return n;
}

UInt_t TFrescoAnalysis::WritePotVals(std::string fname){

  std::ofstream outf(Form("%s/%s",MDIR.c_str(),fname.c_str()));
 
  UInt_t kind,kbpot,kbtype,parnum;
  UInt_t n=0;
    
  for(int i=0; i<(int)HOLDERNAME.size(); i++){
    
    if(HOLDERNAME.at(i).find("<KP")==npos)
      continue;
      
    if(!GetParInfo(HOLDERNAME.at(i),kind,kbpot,kbtype,parnum))
      continue;  

    outf << kbpot << "\t" << kbtype << "\t" << parnum << "\t" << HOLDERVAL.at(i).c_str() << "\n";
    n++;    
  }    

  if(verbose>1) printf("\n\n Wrote %i potential parameters to ' %s ' in ' %s '\n\n",n,fname.c_str(),MDIR.c_str());
  
  return n;
}

UInt_t TFrescoAnalysis::SetPotVals(UInt_t kbpot){
  
// calculate OM
// pot #1 : detueron + 95Sr
// pot #2 : proton   + 96Sr
// pot #3 : neutron  + 95Sr    *FIXED*
// pot #4 : proton   + neutron *FIXED*
// pot #5 : proton   + 95Sr
  
  std::string namebase, name;
  UInt_t pmin=1, pmax=5, kbtype, pn, n=0;
  if(kbpot){
    pmin=kbpot;
    pmax=kbpot;
  }
  if(REAC.find("dp")==npos){ // if elastic only get first potential
    pmin=1;
    pmax=1;    
  }
  
  Double_t par[19];
  Double_t val;
  for(int pot=pmin ; pot<=(int)pmax; pot++){
  
    for(int i=0;i<19;i++)
      par[i]=0.0;

    // if (d,d) or (d,p) potential1
    if(TYPE==0 || (TYPE==2 && pot==1)){ // deut + 95Sr
    
      if(OM1.find("PP1")!=npos){       
        // kbtype = 0 [coulomb]
        par[0] = 1.3    ;
        // kbtype = 1 [volume WS]
        par[1] = 109.45 ; par[2] = 1.05	  ; par[3] = 0.86	  ;     
        // kbtype = 2 [surface WS]               
        par[10] = 10.471; par[11] = 1.43	; par[12] = 0.771 ;
        // kbtype = 3 [Spin Orbit]
        par[13] = 7.0   ; par[14] = 0.75	; par[15] = 0.5   ;        
      } else if(OM1.find("PP2")!=npos){       
        // kbtype = 0 [coulomb]
        par[0] = 1.15   ;
        // kbtype = 1 [volume WS]
        par[1] = 95.29  ; par[2] = 1.15	  ; par[3] = 0.81	  ;      
        // kbtype = 2 [surface WS]               
        par[10] = 16.981; par[11] = 1.34	; par[12] = 0.68  ;
        // * NO SPIN ORBIT*        
      } else if(OM1.find("LH")!=npos){
        // kbtype = 0 [coulomb]
        par[0] = 1.3    ;
        // kbtype = 1 [volume WS]
        par[1] = 109.45 ; par[2] = 1.05	  ; par[3] = 0.86	  ;      
        // kbtype = 2 [surface WS]               
        par[10] = 10.471; par[11] = 1.43	; par[12] = 0.771 ;
        // kbtype = 3 [Spin Orbit]
        par[13] = 7.0   ; par[14] = 0.75	; par[15] = 0.5   ; 
      } else if(OM1.find("D")!=npos){
        // kbtype = 0 [coulomb]
        par[0] = 1.3    ;
        // kbtype = 1 [volume WS]
        par[1] = 93.032 ; par[2] = 1.17	  ; par[3] = 0.727  ;   
        // kbtype = 2 [surface WS]               
        par[10] = 12.336; par[11] = 1.325 ; par[12] = 0.849 ;
        // kbtype = 3 [Spin Orbit]
        par[13] = 7.0 ; par[14] = 1.07	; par[15] = 0.66   ;
      } else {
        printf("\n\t Error :  optical model ' %s ' not found.\n\n",OM1.c_str());
        return n;
      }
    }
    // if (p,p) or (d,p) potential 2 or 5
    if(TYPE==1 || ((pot==2 || pot==5) && TYPE==0)){   // prot + 95Sr
  
      if(OM2.find("PP3")!=npos){       
        // kbtype = 0 [coulomb]
        par[0] = 1.25   ;
        // kbtype = 1 [volume WS]
        par[1] = 58.725 ; par[2] = 1.25	  ; par[3] = 0.65	  ;     
        // kbtype = 2 [surface WS]               
        par[10] = 13.5  ; par[11] = 1.25	; par[12] = 0.47  ;
        // kbtype = 3 [Spin Orbit]
        par[13] = 7.5   ; par[14] = 1.25	; par[15] = 0.47  ;
      } else if(OM2.find("PP4")!=npos){       
        // kbtype = 0 [coulomb]
        par[0] = 1.15   ;
        // kbtype = 1 [volume WS]
        par[1] = 59.599 ; par[2] = 1.17	  ; par[3] = 0.75	  ;      
        // kbtype = 2 [surface WS]               
        par[10] = 12.956; par[11] = 1.32	; par[12] = 0.51  ;
        // kbtype = 3 [Spin Orbit]
        par[13] = 6.2   ; par[14] = 1.01	; par[15] = 0.75  ;
      } else if(OM2.find("BG")!=npos){
        // kbtype = 0 [coulomb]
        par[0] = 1.25    ;
        // kbtype = 1 [volume WS]
        par[1] = 60.41  ; par[2] = 1.17	  ; par[3] = 0.75	  ;      
        // kbtype = 2 [surface WS]               
        par[10] = 12.856; par[11] = 1.32	; par[12] = 0.65  ;
        // kbtype = 3 [Spin Orbit]
        par[13] = 6.2   ; par[14] = 1.01	; par[15] = 0.75   ;    
      } else if(OM2.find("CH")!=npos){
        // kbtype = 0 [coulomb]
        par[0] = 1.266  ;
        // kbtype = 1 [volume WS]
        par[1] = 57.314 ; par[2] = 1.201	; par[3] = 0.69 ;   
        par[4] = 0.558  ; par[5] = 1.201	; par[6] = 0.69 ;         
        // kbtype = 2 [surface WS]               
        par[10] = 10.292; par[11] = 1.238 ; par[12] = 0.69;
        // kbtype = 3 [Spin Orbit]
        par[13] = 5.9   ; par[14] = 1.077	; par[15] = 0.63;
      } else{
        printf("\n\t Error :  optical model ' %s ' not found.\n\n",OM2.c_str());
        return n;
      }
    }

    kbtype=0;
    pn=0;
    val=0.0;
    for(int j=0; j<19; j++) { // par must be 13 in size
      val=par[j]; 
      if(j==0)  kbtype=0;
      else if(j<=6)  kbtype=1;
      else if(j<=12) kbtype=2;
      else if(j>12)  kbtype=3;  
    
      if(j)
        pn = j-6*(kbtype-1);

      name = GetParName(pot,kbtype,pn);
      if(!name.length())
        return n;
        
      HOLDERNAME.push_back(name.c_str());
      HOLDERVAL.push_back(Form("%.3f",val));
      n++;   
    
      //printf("\n j=%i, %s = %.3f",j,name.c_str(),val);
    }
  }
  return n;
}

std::string TFrescoAnalysis::GetParName(UInt_t kbpot, UInt_t kbtype, UInt_t parnum, Bool_t holdername){

  std::string pname="";     
  
  if(kbpot==0 || kbpot>NKBPOT){ // check potential number
    printf("\n\n\t Error :  1< kbpot [=%i] <%i\n\n",kbpot,NKBPOT);
    return pname;
  }
  if(kbtype==0 && parnum!=0){ // check potential type
    printf("\n\n\t Error :  Potential type %i has only 1 parameter, parnum must be 0 [=%i] \n\n",kbtype,parnum);
    return pname;
  } else if((kbtype>0 && parnum==0) || kbtype>3 || parnum>6){
    printf("\n\n\t Error :  1< kbtype [=%i] <4 and 1< parnum [=%i] <7 .\n\n",kbtype,parnum);
    return pname;
  }
  
  // return holder name
  if(holdername){
    pname.assign(Form("<KP%i_%iP%i>",kbpot,kbtype,parnum));      
    return pname;
  }
  
  // return actual parameter name
  static std::vector<std::string> pnames;
  if(!pnames.size()){
    pnames.push_back("RC"); // coulomb
    pnames.push_back("V0"); // vol real
    pnames.push_back("R0");
    pnames.push_back("A0"); 
    pnames.push_back("W0"); // vol imag
    pnames.push_back("RW0");
    pnames.push_back("AW0");     
    pnames.push_back("VD0"); // surf real
    pnames.push_back("RD0");
    pnames.push_back("AD0"); 
    pnames.push_back("WD"); // surf imag
    pnames.push_back("RD");
    pnames.push_back("AD");     
    pnames.push_back("VSO"); // spin orbit real
    pnames.push_back("RSO");
    pnames.push_back("ASO"); 
    pnames.push_back("WSO"); // spin orbit imag
    pnames.push_back("RWSO");
    pnames.push_back("AWSO");     
  }
  
  if(kbtype==0)
    return pnames.at(0);
  else
    return pnames.at(1+(kbtype-1)*6+parnum-1);      
}

Bool_t TFrescoAnalysis::GetParInfo(std::string parname, UInt_t &kind, UInt_t &kbpot, UInt_t &kbtype, UInt_t &parnum){

  if(parname.find("<")!=npos && parname.find(">")!=npos){
    kind = 1;
    kbpot  = (UInt_t)parname[3]-'0';
    kbtype = (UInt_t)parname[5]-'0';
    parnum = (UInt_t)parname[7]-'0';        
    return true;
  }

  static std::string numstr = "012345";
  // return actual parameter name
  static std::vector<std::string> pnames;
  if(!pnames.size()){
    pnames.push_back("RC"); // coulomb
    pnames.push_back("V0"); // vol real
    pnames.push_back("R0");
    pnames.push_back("A0"); 
    pnames.push_back("W0"); // vol imag
    pnames.push_back("RW0");
    pnames.push_back("AW0");     
    pnames.push_back("VD0"); // surf real
    pnames.push_back("RD0");
    pnames.push_back("AD0"); 
    pnames.push_back("WD"); // surf imag
    pnames.push_back("RD");
    pnames.push_back("AD");     
    pnames.push_back("VSO"); // spin orbit real
    pnames.push_back("RSO");
    pnames.push_back("ASO"); 
    pnames.push_back("WSO"); // spin orbit imag
    pnames.push_back("RWSO");
    pnames.push_back("AWSO");     
  }  
  
  for(int i=0; i<(int)pnames.size(); i++){
    
    if(parname.find(pnames.at(i).c_str())!=npos){
      
      kind = 1;
      
      if(numstr.find(parname[0])!=npos)
        kbpot = (UInt_t)numstr.find(parname[0]);
      else 
        kbpot = 1;
      
      if(i==0){
        kbtype = 0;
        parnum = 0;
      } else if(i<=6){
        kbtype=1;
        parnum=i;
      } else if(i<=12){
        kbtype=2;
        parnum=i-6;
      }else if(i>12){
        kbtype=3;
        parnum=i-12;
      }
      return true;
    }
  }
  return false;
}

Bool_t TFrescoAnalysis::UserSetPotVal(UInt_t kbpot, UInt_t kbtype, UInt_t parnum, Double_t val){
  
  std::string vname = GetParName(kbpot,kbtype,parnum);
  if(!vname.length())
    return false;
    
  Int_t indx = GetHolderIndex(vname);
  if(indx>=0){
    if(verbose) printf("\n\n Changed parameter %s [%s] value from %.3f to %.3f\n\n",
          vname.c_str(),GetParName(kbpot,kbtype,parnum,false).c_str(),atof(HOLDERVAL.at(indx).c_str()),val);
    HOLDERVAL.at(indx) = Form("%.3f",val);
  }else{
    HOLDERNAME.push_back(vname);
    HOLDERVAL.push_back(Form("%.3f",val));
  }

  return true;
}

Bool_t TFrescoAnalysis::UserSetPotVal(std::string parname, Double_t val){

  UInt_t kind, kbpot, kbtype, parnum;
  if(!GetParInfo(parname,kind,kbpot,kbtype,parnum))
    return false;

  return UserSetPotVal(kbpot,kbtype,parnum,val);    
}

Bool_t TFrescoAnalysis::CheckInfo(){

  if(!set_info){
    printf("\n\n Error : The reaction must first be defined using the SetInfo function.\n\n");
    return false;   
  }

  return true;
}

UInt_t TFrescoAnalysis::AddFitData(std::string dfile, Double_t thmin, Double_t thmax){

  if(!CheckInfo())
    return false;

  if(!dfile.length()) // look for default file name
    dfile.assign(DNAME);
  if(dfile.find("/")==npos) // use default path
    dfile = DDIR + "/" + dfile;
    
// use local copy of data ?
  std::ifstream file_data(dfile.c_str()); 
  if(!file_data.is_open()){
    printf("\n\n\t Error :  data file '%s' could not be found.\n\n",dfile.c_str());
    return false;
  }
  
  HOLDERNAME.push_back(Form("<FITDATA%i>",NFITDATA));  
  
  NFITDATA++;  
  Int_t indx = GetHolderIndex("<NFITDATA>");
  HOLDERVAL.at(indx)=Form("%i",NFITDATA);
  
  std::string line, added_data="";
  if(TYPE==0 || TYPE==4)
    added_data+="\n&data type=0 iscale=2 idir=0 lab=F abserr=T ic=2/ \n\n";
  else 
    added_data+="\n&data type=0 iscale=-1 idir=2 lab=F abserr=T/ \n\n";
  
  UInt_t n=0;
  Double_t theta, sigma, sigerr;
  
  gdata = new TGraphErrors();
  gdata->SetName(Form("%s_ExpData",FNAME.c_str()));
  std::string units;
  if(TYPE==0)
    units = "mb/sr";
  else   
    units = "Ratio to Rutherford";
  
  gdata->SetTitle(Form("^{%i}Sr(%c,%c) @ %.3f MeV/u; #theta_{CM} [#circ]; #frac{d#sigma}{d#Omega} [%s]",
                  BEAMA,REAC.at(0),REAC.at(1),EPERU,units.c_str()));
  gdata->SetMarkerColor(1);
  gdata->SetLineColor(1);
  gdata->SetMarkerStyle(8);
  gdata->SetMarkerSize(0.7);  
    
  while(std::getline(file_data, line)){
    std::stringstream ss(line);
    ss >> theta >> sigma >> sigerr;
    if((thmin && theta<thmin) || (thmax && theta>thmax)) // ignore points
      continue;
      
    gdata->SetPoint(n,theta,sigma);
    gdata->SetPointError(n,0,sigerr);
      
    if(verbose==3) printf("\n%s",line.c_str());
    added_data+=line;
    added_data+="\n";
    n++;
  }  	   
  Double_t thterr = 10.0, thterr_tmp;
  for(int i=1; i<n; i++){ // get smallest x distance between points
    thterr_tmp = 0.5*(gdata->GetX()[i]-gdata->GetX()[i-1]);
    if(thterr_tmp<thterr)
      thterr = thterr_tmp;
  }   
  for(int i=0; i<n; i++) // set x error 
    gdata->SetPointError(i,thterr,gdata->GetEY()[i]);  
  
  added_data+="\n & \n";
  
  HOLDERVAL.push_back(added_data);  

  if(verbose>1)
    printf("\n\n Added %i lines of experimental data from ' %s '.\n\n",n,dfile.c_str());
  
  // this is now the data file
  DNAME.assign(dfile);  
  
  return n;
}

Bool_t TFrescoAnalysis::AddFitPar(UInt_t kind, UInt_t kbpot, UInt_t kbtype, UInt_t parnum, Double_t val, Double_t valmin, Double_t valmax, Double_t vstep){

  if(!CheckInfo())
    return false;
    
  HOLDERNAME.push_back(Form("<FITPARS%i>",NFITPARS));      

  NFITPARS++;  
  Int_t indx = GetHolderIndex("<NFITPARS>");
  HOLDERVAL.at(indx)=Form("%i",NFITPARS);

  std::string added_par = "";
  if(kind==1){
    std::string parname = GetParName(kbpot,kbtype,parnum,false);  
    if(val<0)
      val = std::atof(HOLDERVAL.at(GetHolderIndex(parname.c_str())).c_str());
    if(valmin<0 || valmax<0){
      valmin = 0.5*val;
      valmax = 2*val;
      vstep = 0.01; // small step size for R and A parameters
      if(parname.find("V")!=npos || parname.find("W")!=npos)
        vstep=0.5; // modify potential depth by larger increments
    }
    added_par+=Form("\n&variable kind=1 name='%s' kp=%i pline=%i col=%i potential=%.2f valmin=%.2f valmax=%.2f step=%.3f/ \n",
        parname.c_str(),kbpot,kbtype+1,parnum,val,valmin,valmax,vstep);
  }else if(kind==2){
    added_par+="\n&variable kind=2 nafrac=1 afrac=1.0 name='Spec. Ampl.' valmin=0.001 valmax=10.0/ \n";
  }else if(kind==5){
    added_par+="\n&variable kind=5 name='exptnorm'/ \n";
  }else {
    printf("\n\n\t Error :  kind must be equal to 1 [Optical Model], 2 [exptnorm] or 5 [Spec. Ampl.]\n\n");
    return false;
  }  
  
  HOLDERVAL.push_back(added_par);
  if(verbose>2)
    printf("\n\n Added the following fit parameter :-\n\t %s \n\n",added_par.c_str());

  return true;
}

Bool_t TFrescoAnalysis::AddFitPar(std::string parname, Double_t val, Double_t valmin, Double_t valmax, Double_t vstep){

  UInt_t kind, kbpot, kbtype, parnum;
  if(!GetParInfo(parname,kind,kbpot,kbtype,parnum))
    return false;
    
  return AddFitPar(kind,kbpot,kbtype,parnum,val,valmin,valmax,vstep);  
}

Int_t TFrescoAnalysis::GetHolderIndex(std::string hname){

  for(int i=0; i<(int)HOLDERNAME.size(); i++)
    if(HOLDERNAME.at(i)==hname)
      return i;

  return -1;
}

Bool_t TFrescoAnalysis::SetFileName(std::string fname){

  FNAME = fname; // set file name
  Int_t indx = GetHolderIndex("<FRESNAME>");
  if(indx>=0){
    HOLDERVAL.at(indx)=fname; // also fresco file name
    return true;
  }else
    return false;
}

UInt_t TFrescoAnalysis::BuildReplacementStrings(){
  // now we build all the replacement strings

  HOLDERNAME.push_back("<FRESNAME>");HOLDERVAL.push_back(FNAME);          
  
  HOLDERNAME.push_back("<BEAMA>");  HOLDERVAL.push_back(Form("%i",BEAMA));  
  HOLDERNAME.push_back("<BMASS>");  HOLDERVAL.push_back(Form("%.4f",BMASS));  
  HOLDERNAME.push_back("<BSPIN>");  HOLDERVAL.push_back(Form("%.1f",BSPIN));
  HOLDERNAME.push_back("<EPERU>");  HOLDERVAL.push_back(Form("%.3f",EPERU));
  HOLDERNAME.push_back("<ELAB>");   HOLDERVAL.push_back(Form("%.1f",ELAB));    
  HOLDERNAME.push_back("<TARGA>");  HOLDERVAL.push_back(Form("%i",TARGA));  
  HOLDERNAME.push_back("<TMASS>");  HOLDERVAL.push_back(Form("%.4f",TMASS));  
  HOLDERNAME.push_back("<TSPIN>");  HOLDERVAL.push_back(Form("%.1f",TSPIN));    
  HOLDERNAME.push_back("<RECOA>");  HOLDERVAL.push_back(Form("%i",RECOA));    
  HOLDERNAME.push_back("<RMASS>");  HOLDERVAL.push_back(Form("%.4f",RMASS));
  
  HOLDERNAME.push_back("<QVAL>");   HOLDERVAL.push_back(Form("%.3f",QVAL));
  HOLDERNAME.push_back("<OM>" );    HOLDERVAL.push_back(OM);  
  
  HOLDERNAME.push_back("<EXC>");    HOLDERVAL.push_back(Form("%.3f",REXC));
  HOLDERNAME.push_back("<JF>" );    HOLDERVAL.push_back(Form("%.1f",JF));
  HOLDERNAME.push_back("<SN>");     HOLDERVAL.push_back(Form("%.4f",RSN));    
  HOLDERNAME.push_back("<BE>" );    HOLDERVAL.push_back(Form("%.3f",NBE));
  HOLDERNAME.push_back("<L>"  );    HOLDERVAL.push_back(Form("%i",L));
  HOLDERNAME.push_back("<JO>" );    HOLDERVAL.push_back(Form("%.1f",JO));
  HOLDERNAME.push_back("<NO>" );    HOLDERVAL.push_back(Form("%i",NO));
  
  HOLDERNAME.push_back("<NFITPARS>");   HOLDERVAL.push_back(Form("%i",NFITPARS));
  HOLDERNAME.push_back("<FITPARS>");    HOLDERVAL.push_back("");    
  HOLDERNAME.push_back("<NFITDATA>");   HOLDERVAL.push_back(Form("%i",NFITDATA));    
  HOLDERNAME.push_back("<FITDATA>");    HOLDERVAL.push_back("");   
      
  if(TYPE==0)
    AddFitPar(2); // spectroscopic factor fit
  else 
    AddFitPar(5); // excptnorm fit

  if(verbose==3){
    printf("\n\n Built wildcard vocabulary with %i strings : -\n",(int)HOLDERNAME.size());
    Print("all");
  }
  
  return (UInt_t)HOLDERVAL.size();
}

Bool_t TFrescoAnalysis::CreateFile(const char *template_name, const char *file_name){
  if(verbose==3)printf("\n\n ------------------------------------------------------------------------------");
  if(verbose>1) printf("\n File ' %s ' will be created : ",file_name);
  if(verbose==3)printf("\n .  .  .  .  . . . .....-------------------...... . . . . . .  .  .  .  .  .  .\n");
  
  
  std::ifstream file_tmp(Form("%s/%s",TDIR.c_str(),template_name));
  if(!file_tmp.is_open()){
    printf("\n\t Error :  template file '%s' could not be found in '%s'.\n",template_name,TDIR.c_str());
    return false;
  }      
  std::ofstream file_out(file_name);

	std::string line;
  while (std::getline(file_tmp, line)) {
		ReplaceHolders(line);
    if(verbose==3)  printf("\n%s",line.c_str());    
    file_out << line << "\n";  
  }
  
  file_tmp.close();
  file_out.close();
  
  if(verbose==3){
    printf("\n .   .   .  .  .  .  . . . ..... Complete ..... . . . .  .  .  .  .   .   .   .");
    printf("\n ------------------------------------------------------------------------------\n\n");
  } else if(verbose>1) printf("    DONE!");
  return true;
}

void TFrescoAnalysis::ReplaceHolders(std::string &line){

  unsigned long pos;
  
  for(int i=0; i<(int)HOLDERNAME.size(); i++){
  	
    pos = line.find(HOLDERNAME.at(i).c_str());
    if(pos==npos)
      continue;
      
    // for fit data and parameters we push all of data/pars into a single string to be replaced all at once 
    if(HOLDERNAME.at(i).find("<FITDATA>")!=npos){ // put all fit data into this holder
      HOLDERVAL.at(i) = ""; // clear it first
      for(int j=0; j<NFITDATA; j++)
        HOLDERVAL.at(i)+=HOLDERVAL.at(GetHolderIndex(Form("<FITDATA%i>",j)));
    } else if(HOLDERNAME.at(i).find("<FITPARS>")!=npos){ // put all fit pars into this holder
      HOLDERVAL.at(i) = ""; // clear it first
      for(int j=0; j<NFITPARS; j++)
        HOLDERVAL.at(i)+=HOLDERVAL.at(GetHolderIndex(Form("<FITPARS%i>",j)));
    }
       
    while(pos!=npos){
      line.replace(pos,HOLDERNAME.at(i).length(),HOLDERVAL.at(i).c_str());
      pos = line.find(HOLDERNAME.at(i).c_str());
    }	 
    
    if(HOLDERNAME.at(i).find("<FITDATA>")!=npos || HOLDERNAME.at(i).find("<FITPARS>")!=npos)
      HOLDERVAL.at(i) = ""; // clear it
  }
}

TList *TFrescoAnalysis::MakePotentials(Int_t kbpot){

  Double_t rpars[10], ipars[10];
  // par[0] is coulomb
  // pars[1,2,3] are volume
  // pars[4,5,6] are surface
  // pars[7,8,9] are L.S
  
  // assume WS potential
  Int_t kbtype, par, indx;
  unsigned long pos;  
  Double_t val; 
  for(int i=0; i<(int)HOLDERNAME.size(); i++){ 
    std::string pname = HOLDERNAME.at(i);
    if(pname.find(Form("<KP%i",kbpot))==npos) // get to potentials
      continue;
    
    pos = pname.find("_");
    if(pos==npos)
      continue;
      
//    printf("\n\t %s \t=\t %s",holdername.at(i).c_str(),holderval.at(i).c_str());   
    // < KP %1 _ %2 P %3 >  
    // %1 - potential
    // %2 - kbtype
    // %3 - par
    kbtype = atoi(pname.substr(pos+1,1).c_str());
    par  = atoi(pname.substr(pos+3,1).c_str());  
     
    val  = atof(HOLDERVAL.at(i).c_str());
    
    if(kbtype==0){ // coulomb radius
      rpars[0] = val;
      continue;
    }
       
    if(par<=3) // REAL POTENTIAL, p1/p2/p3
      rpars[par+3*(kbtype-1)] = val;  
    else       // IMAG POTENTIAL, p4/p5/p6
      ipars[par+3*(kbtype-1)-3] = val;  
  }    
  

  
  TGraph *greal = new TGraph();
  greal->SetNameTitle(Form("real_pot_%i",kbpot),Form("RealPotential_%i; R [fm]; V [MeV]",kbpot));
  greal->SetLineColor(1);
  greal->SetLineWidth(1);
  
  TGraph *gr_co = new TGraph();
  gr_co->SetNameTitle(Form("coul_pot_%i",kbpot),Form("CoulombPotential_%i; R [fm]; V [MeV]",kbpot)); 
  gr_co->SetLineColor(5);  
  gr_co->SetLineWidth(1);  
  
  TGraph *gr_ll = new TGraph();
  gr_ll->SetNameTitle(Form("centrifugal_pot_%i",kbpot),Form("CentrifugalPotential_%i; R [fm]; V [MeV]",kbpot));
  gr_ll->SetLineColor(6);  
  gr_ll->SetLineWidth(1);
   
  TGraph *gr_vo = new TGraph();
  gr_vo->SetNameTitle(Form("rvol_pot_%i",kbpot),Form("RealVolumePotential_%i; R [fm]; V [MeV]",kbpot));    
  gr_vo->SetLineColor(2);  
  gr_vo->SetLineWidth(1);
 
  TGraph *gr_su = new TGraph();
  gr_su->SetNameTitle(Form("rsur_pot_%i",kbpot),Form("RealSurfacePotential_%i; R [fm]; V [MeV]",kbpot)); 
  gr_su->SetLineColor(3);  
  gr_su->SetLineWidth(1);
  
  TGraph *gr_ls = new TGraph();
  gr_ls->SetNameTitle(Form("rspo_pot_%i",kbpot),Form("RealSpinOrbitPotential_%i; R [fm]; V [MeV]",kbpot)); 
  gr_ls->SetLineColor(4);  
  gr_ls->SetLineWidth(1);
  
  TGraph *gimag = new TGraph();
  gimag->SetNameTitle(Form("imag_pot_%i",kbpot),Form("ImagPotential_%i; R [fm]; V [MeV]",kbpot));
  gimag->SetLineColor(1);
  gimag->SetLineStyle(2);
  gimag->SetLineWidth(1);
  
  
  TGraph *gi_vo = new TGraph();
  gi_vo->SetNameTitle(Form("ivol_pot_%i",kbpot),Form("ImagVolumePotential_%i; R [fm]; V [MeV]",kbpot));   
  gi_vo->SetLineColor(2);      
  gi_vo->SetLineStyle(2);
  gi_vo->SetLineWidth(1);
     
  TGraph *gi_su = new TGraph();
  gi_su->SetNameTitle(Form("isur_pot_%i",kbpot),Form("ImagSurfacePotential_%i; R [fm]; V [MeV]",kbpot)); 
  gi_su->SetLineColor(3);      
  gi_su->SetLineStyle(2);
  gi_su->SetLineWidth(1);
  
  TGraph *gi_ls = new TGraph();
  gi_ls->SetNameTitle(Form("ispo_pot_%i",kbpot),Form("ImagSpinOrbitPotential_%i; R [fm]; V [MeV]",kbpot));   
  gi_ls->SetLineColor(4);    
  gi_ls->SetLineStyle(2);
  gi_ls->SetLineWidth(1);
    
  
  Double_t x, xmax, Aval, LS = 0.5*(JO*(JO+1)-L*(L+1)-0.75);
  Double_t Vr_co, Vr_ll, Vr_vo, Vr_su, Vr_ls, Vi_vo, Vi_su, Vi_ls;

  if(kbpot==2)
    Aval = pow(RECOA,1.0/3.0);
  else 
    Aval = pow(BEAMA,1.0/3.0);    
    
   xmax = 4*Aval;  
    
  for(int i=0; i<100; i++){
  
    x = xmax*(double)(i+1)/100.0;
    // REAL
    Vr_co = Coulomb(x,rpars[0]*Aval); // coulomb
    Vr_ll = L*(L+1)*0.4379/(2*x*x); // ang. mom. barrier
    Vr_vo = WS_Volume(x,rpars[1],rpars[2]*Aval,rpars[3]);  // WS volume
    Vr_su = WS_Volume(x,rpars[4],rpars[5]*Aval,rpars[6]);  // WS surface
    Vr_ls = WS_Surface(x,rpars[7],rpars[8]*Aval,rpars[9])*LS/x; // modified WS surface
    // IMAG
    Vi_vo = WS_Volume(x,ipars[1],ipars[2]*Aval,ipars[3]);  // WS volume
    Vi_su = WS_Volume(x,ipars[4],ipars[5]*Aval,ipars[6]);  // WS surface
    Vi_ls = WS_Surface(x,ipars[7],ipars[8]*Aval,ipars[9])*LS/x; // modified WS surface
            
    gr_co->SetPoint(i,x,Vr_co);
    gr_ll->SetPoint(i,x,Vr_ll);
    gr_vo->SetPoint(i,x,Vr_vo);
    gr_su->SetPoint(i,x,Vr_su);
    gr_ls->SetPoint(i,x,Vr_ls);
    greal->SetPoint(i,x,Vr_co + Vr_ll + Vr_vo + Vr_su + Vr_ls);

    gi_vo->SetPoint(i,x,Vi_vo);
    gi_su->SetPoint(i,x,Vi_su);
    gi_ls->SetPoint(i,x,Vi_ls);
    gimag->SetPoint(i,x,Vi_vo + Vi_su + Vi_ls);
    
  }
 
  TList *list = new TList();
  
  
  if(CheckState(greal))list->Add(greal);
  if(CheckState(gr_co))list->Add(gr_co);
  if(CheckState(gr_vo))list->Add(gr_vo);
  if(CheckState(gr_su))list->Add(gr_su);
  if(CheckState(gr_ls))list->Add(gr_ls);  
  if(CheckState(gr_ll))list->Add(gr_ll);
  if(CheckState(gimag))list->Add(gimag);
  if(CheckState(gi_vo))list->Add(gi_vo);
  if(CheckState(gi_su))list->Add(gi_su);
  if(CheckState(gi_ls))list->Add(gi_ls);
  
  for(int i=0; i<100; i++){
  
    x = xmax*(double)(i+1)/100.0;
    // REAL
    Vr_co = Coulomb(x,rpars[0]*Aval); // coulomb
    Vr_ll = L*(L+1)*0.4379/(2*x*x); // ang. mom. barrier
    Vr_vo = WS_Volume(x,rpars[1],rpars[2]*Aval,rpars[3]);  // WS volume
    Vr_su = WS_Surface(x,rpars[4],rpars[5]*Aval,rpars[6]);  // WS surface
    Vr_ls = WS_Surface(x,rpars[7],rpars[8]*Aval,rpars[9])*LS/x; // modified WS surface
    // IMAG
    Vi_vo = WS_Volume(x,ipars[1],ipars[2]*Aval,ipars[3]);  // WS volume
    Vi_su = WS_Surface(x,ipars[4],ipars[5]*Aval,ipars[6]);  // WS surface
    Vi_ls = WS_Surface(x,ipars[7],ipars[8]*Aval,ipars[9])*LS/x; // modified WS surface
            
    gr_co->SetPoint(i,x,Vr_co);
    gr_ll->SetPoint(i,x,Vr_ll);
    gr_vo->SetPoint(i,x,Vr_vo);
    gr_su->SetPoint(i,x,Vr_su);
    gr_ls->SetPoint(i,x,Vr_ls);
    greal->SetPoint(i,x,Vr_co + Vr_ll + Vr_vo + Vr_su + Vr_ls);

    gi_vo->SetPoint(i,x,Vi_vo);
    gi_su->SetPoint(i,x,Vi_su);
    gi_ls->SetPoint(i,x,Vi_ls);
    gimag->SetPoint(i,x,Vi_vo + Vi_su + Vi_ls);
    
  }  
   
  return list;
}

Double_t TFrescoAnalysis::Coulomb(Double_t x, Double_t R){
  Double_t val;
  if(x<R)
    val = 38*1.44/(2*R)*(3-x*x/(R*R));
  else  
    val = 38*1.44/x;
    
  if(TMath::IsNaN(val) || fabs(val)>1e6)
    return 0;   
  else
    return val;      
}

Double_t TFrescoAnalysis::WS_Volume(Double_t x, Double_t V, Double_t R, Double_t A){
  Double_t val =  -V/(1+exp((x-R)/A));
  
  if(TMath::IsNaN(val) || fabs(val)>1e6)
    return 0;
  else
    return val; 
}

Double_t TFrescoAnalysis::WS_Surface(Double_t x, Double_t V, Double_t R, Double_t A){
  Double_t fx = exp((x-R)/A);  
  Double_t val = -V/A*fx/pow(1+fx,2.0);
  
  if(TMath::IsNaN(val) || fabs(val)>1e6)
    return 0;
  else
    return val;     
}

void TFrescoAnalysis::Print(Option_t *opt){

  printf("\n\n\t\t ___TFrescoAnalysis___\n");
  if(!set_info){
    printf("\n\n Warning :  Set the reaction using SetInfo function\n\n");
    return;
  }
  
  printf("\n Reaction :  %iSr(%c,%c)%iSr @ %.3f MeV/u [beame = %.1f MeV]",
    BEAMA,REAC.at(0),REAC.at(1),RECOA,EPERU,ELAB);
  
  if(TYPE==0)
    printf("\n\n   %iSr   :  Exc = %.3f MeV,  I=%.1f x [ L=%i , j=%.0f/2 ]-> J=%.1f, OM = '%s'[p]+'%s'[d]: -\n",
            RECOA,REXC,BSPIN,L,2*JO,JF,OM2.c_str(),OM1.c_str()); 
            
  if(strcmp(opt,"all")==0)
    for(int i=0; i<HOLDERNAME.size(); i++)
      printf("\n%i.\t%20s\t=\t%s",i,HOLDERNAME.at(i).c_str(),HOLDERVAL.at(i).c_str());
              
  if(strcmp(opt,"fit")==0){
    for(int i=0; i<NFITPARS; i++){
      UInt_t indx = GetHolderIndex(Form("<FITPARS%i>",i));
      printf("\n%i.\t%20s\t=\t%s",indx,HOLDERNAME.at(indx).c_str(),HOLDERVAL.at(indx).c_str());
    }
    for(int i=0; i<NFITDATA; i++){
      UInt_t indx = GetHolderIndex(Form("<FITDATA%i>",i));
      printf("\n%i.\t%20s\t=\t%s",indx,HOLDERNAME.at(indx).c_str(),HOLDERVAL.at(indx).c_str());
    }
  }
  
  std::string optstr;
  optstr.assign(opt);
  if(optstr.length()){
    for(int i=0; i<(int)HOLDERVAL.size(); i++){
      if(HOLDERNAME.at(i).find(opt)!=npos){ // find any matching string
        if(HOLDERNAME.at(i).find("KP")!=npos){
          UInt_t kind=1, kbpot, kbtype, parnum; // print real potential name
          GetParInfo(HOLDERNAME.at(i),kind,kbpot,kbtype,parnum);
          std::string pname = GetParName(kbpot,kbtype,parnum,false);
          printf("\n%i.\t%10s%10s\t=\t%s",i,HOLDERNAME.at(i).c_str(),pname.c_str(),HOLDERVAL.at(i).c_str());        
        } else    
          printf("\n%i.\t%20s\t=\t%s",i,HOLDERNAME.at(i).c_str(),HOLDERVAL.at(i).c_str());
      }
    }  
  }
  
  if(NORM)
    printf("\n\n\t  Normalization         : %.3e +/- %.3e",NORM,NORM_ERR);
  
  if(GAM)
    printf("\n\n\t  TIGRESS Efficiency    : %.3e +/- %.3e @ %.1f keV",GAMEFF,GAMEFF_ERR,GAM);
  
  if(SF)
    printf("\n\n\t  Spectroscopic Factor  : %.3e +/- %.3e",SF,SF_ERR);    
  
  if(CS)
    printf("\n\n\t  Total Cross Section   : %.3e +/- %.3e",CS,CS_ERR);   
    
  if(CHI2)
    printf("\n\n\t  ChiSq/N               : %.3e",CHI2);
        
  printf("\n\n");  
}
 
void TFrescoAnalysis::Clear(Option_t *opt){

  set_info = false;
  TYPE = -1;
    
  BEAMA=0;
  BMASS=0.0;
  BSPIN=0.0;
  EPERU=0.0;
  ELAB=0.0;
  QVAL=0.0;
  
  TARGA=0;
  TMASS=0.0;
  
  RECOA=0;
  RMASS=0.0;
  REXC=0.0;
  RSN=0.0;
  
  OM.clear();
  OM1.clear();
  OM2.clear();

  REAC.clear();
  
  HOLDERNAME.clear();
  HOLDERVAL.clear();
  NFITDATA = 0;
  NFITPARS = 0;  
  
  JF=0.0;
  L=0;
  JO=0.0;    
  NO=0.0;
  NBE=0.0;  
  
  EXCHI=0.0;   
  NORM=0.0;
  NORM_ERR=0.0;
  GAM=0.0;
  GAMEFF=0.0;
  GAMEFF_ERR=0.0;     

  SA=0.0;
  SA_ERR=0.0;
  SF=0.0;
  SF_ERR=0.0;  
  CS=0.0;
  CS_ERR=0.0;

  CHI2=0.0;  

}    

void TFrescoAnalysis::GetRid(const char *name, Bool_t delete_all){
	TObject *obj;

	obj = gROOT->FindObjectAny(name);
	if(obj) gROOT->Remove(obj);

  if(delete_all)
    return;
}

Bool_t TFrescoAnalysis::CheckState(TGraph *g){

  double ytot=0.0;
  
  for(int i=0; i<g->GetN(); i++)
    ytot = ytot + fabs(g->GetY()[i]);
         
  if(TMath::IsNaN(ytot) || ytot<1.0)
    return false;
  
   return true; 
}
