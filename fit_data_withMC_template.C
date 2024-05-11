#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TCut.h>
#include <TStyle.h>
#include <TF1.h>

TH1F *h_prompt_MC;
TH1F *h_nonprompt_MC;

Double_t ftotal(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(xx);
	Double_t br = par[1]*(1-par[0])*h_nonprompt_MC->GetBinContent(bin)/h_nonprompt_MC->Integral();
	Double_t sr = par[1]*par[0]*h_prompt_MC->GetBinContent(bin)/h_prompt_MC->Integral();
	return sr+ br;


}

Double_t funNonPrompt(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(xx);
	Double_t br = par[1]*(1-par[0])*h_nonprompt_MC->GetBinContent(bin)/h_nonprompt_MC->Integral();
	return br;


}

Double_t funPrompt(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(xx);
	Double_t sr = par[1]*par[0]*h_prompt_MC->GetBinContent(bin)/h_prompt_MC->Integral();
	return sr;


}

void fit_data_withMC_template(TString promptF = "output_prompt.root", TString nonpromptF = "output_nonprompt.root", TString dataF = "output_data.root"){

TFile *promptMC_DCA = new TFile(promptF);
TFile *nonpromptMC_DCA = new TFile(nonpromptF);
TFile *data_DCA = new TFile(dataF);


h_prompt_MC = (TH1F*) promptMC_DCA->Get("h_dcas")->Clone("h_prompt_MC");
h_nonprompt_MC = (TH1F*) nonpromptMC_DCA->Get("h_dcas")->Clone("h_nonprompt_MC");
TH1F *h_data_DCA = (TH1F*) data_DCA->Get("h_dcas")->Clone("h_data_DCA");


TH1F *h_nonprompt_d = (TH1F*)h_nonprompt_MC->Clone("h_nonprompt_d");
TH1F *h_prompt_d = (TH1F*)h_prompt_MC->Clone("h_prompt_d");

TF1 *ftot = new TF1("ftot",ftotal,0,40,2);
double total_yield=h_data_DCA->Integral();
cout<<total_yield<<endl;
ftot->SetParameter(0,0.7);
ftot->SetParameter(1,total_yield);
ftot->SetParLimits(0,0,1);
ftot->SetParLimits(1,0,2*total_yield);


h_data_DCA->Fit("ftot","b ");

cout<<"Prompt fraction = "<<ftot->GetParameter(0)<<" ± "<<ftot->GetParError(0)<<endl;
cout<<"Normalized chi2 = "<<(ftot->GetChisquare()/ftot->GetNDF())<<endl;

}


