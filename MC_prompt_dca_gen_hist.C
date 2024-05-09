#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TCut.h>
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"

#include <cstdlib> 
#include <map>
#include <string>
#include <fstream>
#include "TString.h"
#include "TObjString.h"

#include "/home/chand140/Continue_Event_shape_engineering/Event_shape_engineering/Analysis/CMSSW_10_3_3_patch1/src/prompt_D0/include/q2_bin_distribution.h" //latest q2 distribution

using namespace std;

#define  MAX_XB    1500
#define N_DCABINS 15

void divideBin(TH1* h)
{
  h->Sumw2();
  for(int i=1;i<=h->GetNbinsX();i++)
  {
    Float_t contentvalue=h->GetBinContent(i);
    Float_t content_err=h->GetBinError(i);
    contentvalue/=h->GetBinWidth(i);
    content_err/=h->GetBinWidth(i);
    h->SetBinContent(i,contentvalue);
    h->SetBinError(i,content_err);
  }
}

/*Double_t GetBDTCut(Double_t pt, Int_t cent) {
  Double_t bdtcut = -1.;

  //Cuts after changing |y|<1 and DtrkpT > 1~GeV/c below (new production)  
  Double_t bdtcut_0_10[10] = {0.29,0.36,0.37,0.33,0.25,0.16,0.15,0.10,0.00,0.00};
  Double_t bdtcut_10_30[10] = {0.29,0.37,0.36,0.29,0.20,0.14,0.11,0.06,-0.04,-0.06};
  Double_t bdtcut_30_50[10] = {0.28,0.34,0.32,0.23,0.16,0.09,0.05,0.00,-0.05,-0.05};

  Double_t ptBDTbin[11] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};

  for (int i=0; i<10; i++) {
    if (cent >= 0 && cent < 20 && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_0_10[i];
    if (cent >= 20 && cent < 60 && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_10_30[i];
    if (cent >= 60 && cent < 100 && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_30_50[i];
  }
  return bdtcut;
}*/

//New function returns BDT value for a given centrality

Double_t GetBDTCut(Double_t pt, string cent_name){ //(Double_t pt, Int_t cent) {
  Double_t bdtcut = -1.;

  //Cuts after changing |y|<1 and DtrkpT > 1~GeV/c below (new production)  
  Double_t bdtcut_0_10[10] = {0.29,0.36,0.37,0.33,0.25,0.16,0.15,0.10,0.00,0.00};
  Double_t bdtcut_10_30[10] = {0.29,0.37,0.36,0.29,0.20,0.14,0.11,0.06,-0.04,-0.06};
  Double_t bdtcut_30_50[10] = {0.28,0.34,0.32,0.23,0.16,0.09,0.05,0.00,-0.05,-0.05};

  Double_t ptBDTbin[11] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};

  for (int i=0; i<10; i++) {
    if (cent_name=="0_10" && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_0_10[i];
    if (cent_name=="10_20" && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_10_30[i];
    if (cent_name=="20_30" && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_10_30[i];
    if (cent_name=="30_40" && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_30_50[i];
    if (cent_name=="40_50" && ptBDTbin[i]<=pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_30_50[i];
  }
  return bdtcut;
}

int MC_prompt_dca_gen_hist(TString sc_val="", TString cent_name_input="")
{
  //N_CENTBINS = 1; //for 0-10% centrality only
  //TString out_file = "DCA_dist";
  //TFile *outf = new TFile(Form(out_file ),"recreate");

  // Define the histogram to plot on
  //
  /*double edge[9]={0,0.0016,0.0032,0.0048,0.0064,0.008,0.0096,0.016,0.16};
    TH1D *MC_prompt = new TH1D("MC_prompt","MC_prompt",8,edge);
    TH1D *MC_nonprompt = new TH1D("MC_nonprompt","MC_nonprompt",8,edge);*/
  float sc = sc_val.Atof();
  string cent_name = cent_name_input.Data();

  TVector3 pgen1, dir;
  double Dgenalpha;
  double Dgendl;

  //Set Variables for branch address
  int   Dgen[MAX_XB];
  int   parent_id[MAX_XB];
  int   candSize;
  int   cen_val;
  float mass_val[MAX_XB];
  float pT_val[MAX_XB];
  float dca_val[MAX_XB];
  float y_val[MAX_XB];

  float Dgenpt[MAX_XB];
  float Dgeneta[MAX_XB];
  float Dgenphi[MAX_XB];
  float DgenprodvtxX[MAX_XB];
  float DgenprodvtxY[MAX_XB];
  float DgenprodvtxZ[MAX_XB];
  float DgendecayvtxX[MAX_XB];
  float DgendecayvtxY[MAX_XB];
  float DgendecayvtxZ[MAX_XB];
  float PVx;
  float PVy;
  float PVz;

  float Dl3D[MAX_XB];
  float Dl3DErr[MAX_XB];
  float alpha[MAX_XB];
  float Dtrk1Chi2n[MAX_XB];
  float Dtrk2Chi2n[MAX_XB];

  double BDT_weight[MAX_XB];

  //TString filename = "/scratch/bell/chand140/ESE_output_Dec_26/MC_output/MC_prompt_with_BDT.root";
  string filename = "/scratch/bell/chand140/ESE_output_Mar_11_2024/MC_output/MC_prompt_with_BDT_cen"+cent_name+".root";
  if (!TFile::Open(filename.c_str()))   { cout << " fail to open file" << endl; return 0;}	 
  TFile* f = TFile::Open(filename.c_str());
  cout << "Opened file " << filename <<endl;
  TTree* t1 = (TTree*) f->Get("mvaTree");

  //Set Branch Address
  t1->SetBranchAddress("candSize",&candSize);
  t1->SetBranchAddress("centrality",&cen_val);
  t1->SetBranchAddress("pT",&pT_val);
  t1->SetBranchAddress("y",&y_val);
  t1->SetBranchAddress("mass",&mass_val);
  t1->SetBranchAddress("dca",&dca_val);
  t1->SetBranchAddress("Dgen",&Dgen);
  t1->SetBranchAddress("Dtrk1Chi2n",&Dtrk1Chi2n);
  t1->SetBranchAddress("Dtrk2Chi2n",&Dtrk2Chi2n);
  t1->SetBranchAddress("DgenBAncestorpdgId",&parent_id);

  t1->SetBranchAddress("Dgenpt",&Dgenpt);
  t1->SetBranchAddress("Dgeneta",&Dgeneta);
  t1->SetBranchAddress("Dgenphi",&Dgenphi);
  t1->SetBranchAddress("DgenprodvtxX",&DgenprodvtxX);
  t1->SetBranchAddress("DgenprodvtxY",&DgenprodvtxY);
  t1->SetBranchAddress("DgenprodvtxZ",&DgenprodvtxZ);
  t1->SetBranchAddress("DgendecayvtxX",&DgendecayvtxX);
  t1->SetBranchAddress("DgendecayvtxY",&DgendecayvtxY);
  t1->SetBranchAddress("DgendecayvtxZ",&DgendecayvtxZ);
  t1->SetBranchAddress("PVx",&PVx);
  t1->SetBranchAddress("PVy",&PVy);
  t1->SetBranchAddress("PVz",&PVz);
  t1->SetBranchAddress("3DDecayLength",&Dl3D);
  t1->SetBranchAddress("3DDecayLengthSignificance",&Dl3DErr);
  t1->SetBranchAddress("3DPointingAngle",&alpha);
  t1->SetBranchAddress("BDT_weight",&BDT_weight);

  string output_name = "MC_DCA_template/MC_prompt_cen"+cent_name+"_dca_hist_with_res_"+to_string(sc)+".root"; //MC_prompt_cen"+cent_name+"_with_res_"+to_string(sc)+".root";
  TFile *outf = new TFile(output_name.c_str(),"RECREATE");

  TH1D *hist_MC_prompt_dca[N_PTBINS];
  //TH1D *hist_MC_nonprompt_dca[N_CENTBINS][N_PTBINS];
  //TH1D *hist_data_dca[N_PTBINS];
  Double_t dca_edges[N_DCABINS+1]={0,0.00109,0.0024,0.0039,0.0059,0.0085,0.01179,0.016,0.0214,0.028,0.0366,0.0475,0.0615,0.079,0.10145,0.135};

  //auto inf_data = TFile::Open("dca_hist_output.root");

  //dca_edges[N_DCABINS+1]={0,0.00109,0.0024,0.0039,0.0059,0.0085,0.012,0.145};
  //for (int i_cen=0; i_cen<1; i_cen++) {
    
  for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
    hist_MC_prompt_dca[i_pt] = new TH1D("hist_MC_prompt_dca_cen_"+cent_name+"_pt_"+pt_name[i_pt],"hist_MC_prompt_dca_cen_"+cent_name+"_pt_"+pt_name[i_pt],N_DCABINS,dca_edges);
  }

  //}

  /*TString nt_name;
    TNtuple *nt[N_CENTBINS][N_PTBINS];

    for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
    nt_name = "ntuple_cen_"+cen_name[i_cen]+"_pt_"+pt_name[i_pt];
    nt[i_cen][i_pt] = new TNtuple(nt_name,nt_name,"mass:dca");
    }
    }*/

  int n_entries = t1->GetEntries();
  cout << "Reading input MC file \n" << endl;
  for (int i_entry = 0; i_entry < n_entries; i_entry++)
  {
    //if (i_entry > 100) break;
    t1->GetEntry(i_entry);
    if(i_entry%40000==0) cout <<i_entry<<" / "<<n_entries<< "  "<<100*i_entry/n_entries<<"%"<<endl;

    //for (int i_cen=0; i_cen<1; i_cen++)
    //{
      //if (cen_val>=cen_edges[i_cen] && cen_val<cen_edges[i_cen+1])  //centrality cuts taken from the input file here
      //{
        for (int i_cand=0; i_cand<candSize; i_cand++)
        {
          if(fabs(y_val[i_cand])>=1.0 || Dtrk1Chi2n[i_cand]>=0.18 || Dtrk2Chi2n[i_cand]>=0.18 || abs(parent_id[i_cand])>400 || alpha[i_cand]>=0.2) continue;

          for(int i_pt=0; i_pt<N_PTBINS; i_pt++) 
          {        
            if(pT_val[i_cand]<pt_edges[i_pt] || pT_val[i_cand]>=pt_edges[i_pt+1] || Dgen[i_cand]!=23333 || BDT_weight[i_cand]<GetBDTCut(pT_val[i_cand],cent_name)) continue; //pT/cen/Dgen/BDT cuts

            pgen1.SetPtEtaPhi(Dgenpt[i_cand],Dgeneta[i_cand],Dgenphi[i_cand]);
            dir.SetXYZ(DgendecayvtxX[i_cand]-DgenprodvtxX[i_cand],DgendecayvtxY[i_cand]-DgenprodvtxY[i_cand],DgendecayvtxZ[i_cand]-DgenprodvtxZ[i_cand]);
            Dgenalpha = dir.Angle(pgen1);
            Dgendl = dir.Mag();

            // to check the smearing code
            /*if(pT_val[i_cand]>=20 && pT_val[i_cand]<30 &&  cen_val>=40 && cen_val<60) 
              {
              cout << " Dgenalpha: "<< Dgenalpha << " Dgendl: " << Dgendl << " DCA: " << dca_val[i_cand] << " scaled DCA: " << Dgendl*sin(Dgenalpha)+ sc*(dca_val[i_cand]-Dgendl*sin(Dgenalpha)) << endl;
              }*/

            hist_MC_prompt_dca[i_pt]->Fill(Dgendl*sin(Dgenalpha)+ sc*(dca_val[i_cand]-Dgendl*sin(Dgenalpha)));
          }		
        }
      //}
    //}
  } //for loop i_entry

  f->Close();
  cout << "Closed file "<<filename<<endl;

  //divideBin(MC_prompt);
  //divideBin(MC_nonprompt);
  outf->cd();
  //for (int i_cen=0; i_cen<1; i_cen++) {
  for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
    divideBin(hist_MC_prompt_dca[i_pt]);
    hist_MC_prompt_dca[i_pt]->Write();
  }
  //}
  outf->Write(0,TObject::kOverwrite);
  cout<<"--- Writing finished"<<endl;
  outf->Close();
  cout <<"\nSuccessful!!!\n";
  return 0;
}

int main(int argc, char *argv[])
{
    if(argc==3)
    {
      MC_prompt_dca_gen_hist(argv[1],argv[2]);
    }
    else
    {
      std::cout << "Error in arguments" << std::endl;
      return 1;
    }
    return 0;
}