#define HWWAnalysis_cxx
// The class definition in HWWAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("HWWAnalysis.C")
// root> T->Process("HWWAnalysis.C","some options")
// root> T->Process("HWWAnalysis.C+")
//


#include "HWWAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <map> 
#include <math.h>
#include <stdio.h>
#include "TROOT.h"
#include "TObjString.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>

// Define som global variables

TString option; 
TString channel; 

vector<TString> channels;  
vector<TString> categories;

int events_tot = 0;  
int events_ee = 0; 
int events_uu = 0;
int events_eu;

Int_t DSID = 0; 
Int_t prev_DSID = 0; 

TLorentzVector l1; 
TLorentzVector l2; 
TLorentzVector dileptons; 
TLorentzVector MeT;

Float_t wgt; 

void HWWAnalysis::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    events_tot = 0;  
    events_ee = 0; 
    events_uu = 0;
    events_eu = 0;

    DSID = 0; 
    prev_DSID = 0; 

    gROOT->SetBatch(kTRUE);

    option = GetOption(); // Print TString from Process(obj, str)

    lumi = 10.06;

    categories = {"Diboson", "Higgs", "Wjetsincl", "Zjetsincl", "Zjets", "Wjets", "singleTop", "topX", "ttbar"};
    
    CreateHistograms();

}

void HWWAnalysis::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    option = GetOption();

}

Bool_t HWWAnalysis::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // When processing keyed objects with PROOF, the object is already loaded
    // and is available via the fObject pointer.
    //
    // This function should contain the \"body\" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.

    Int_t nbin = 2; // used for adding to cutflow histogram

    fReader.SetLocalEntry(entry);

    // Count events 
    events_tot++;  	
    
    DSID = *channelNumber;
    
    type = GetType(DSID);
    
    isData = false;
    if (type == "data") isData = true;

    // Calculate event weight 
    wgt = 1.0; 
    if( !isData )
    {
        wgt = (*mcWeight)*(*scaleFactor_PILEUP)*(*scaleFactor_ELE)*(*scaleFactor_MUON)*(*scaleFactor_LepTRIGGER)*((*XSection*lumi*1000.)/(*SumWeights)); 
    }

    //------------------------------//
    // Event selection & Kinematics //
    //------------------------------//

    // Require (exactly) 2 leptons
    /*if(*lep_n < 2){ return kTRUE; }

    // Identify the leptons 
    channel = "eu";
    //if(fabs(lep_type[0])==11){ channel = "ee";  } // Electrons  
    //if(fabs(lep_type[0])==13){ channel = "uu";  } // Muons 

    // Require opposite charge
    if(lep_charge[0] == lep_charge[1]){ return kTRUE; } 

    // Require opposite flavour (electron + muon)
    if(lep_type[0] == lep_type[1]){ return kTRUE; } 

    // Set Lorentz vectors 
    l1.SetPtEtaPhiE(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    l2.SetPtEtaPhiE(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);

    for (unsigned int i = 0; i < *lep_n; i++)
    {
        // Cut on z_0 impact parameter 
        if(fabs(lep_z0[i])*TMath::Sin(l1.Theta()) > 0.5){ return kTRUE; }
        
        // Lepton isolation cut 
        if( lep_etcone20[i]/lep_pt[i] > 0.1 && lep_ptcone30[i]/lep_pt[i] > 0.1 ){ return kTRUE; } 
        
        // Require tight id
        if ( !lep_isTightID->at(i)) {return kTRUE; }
        
        // general leptons cuts
        if (lep_pt[i] < 15000) {return kTRUE;}
        if (lep_type[i] == 11)
        {
            if (fabs(lep_eta[i]) > 2.47) {return kTRUE;}
            if( TMath::Abs(lep_trackd0pvunbiased[i])/lep_tracksigd0pvunbiased[i] > 5) {return kTRUE; }
        }
        if (lep_type[i] == 13)
        {
            if (fabs(lep_eta[i]) > 2.5) {return kTRUE;}
            if( TMath::Abs(lep_trackd0pvunbiased[i])/lep_tracksigd0pvunbiased[i] > 3) {return kTRUE; }
        }
    }

    // Sort Lorentz vector (leading lepton = l1, subleading lepton = l2)  
    if(lep_pt[1]>lep_pt[0])
    {
        l1.SetPtEtaPhiE(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);
        l2.SetPtEtaPhiE(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    }
    
    if (lep_pt[0] < 22000) {return kTRUE; }
    // Variables are stored in the TTree with unit MeV, so we need to divide by 1000 
    // to get GeV, which is a more practical and commonly used unit. 
    
    if (*met_et < 30000.) {return kTRUE; }
    
    dileptons = l1 + l2; 
    
    if (dileptons.M() < 10000 || dileptons.M() > 55000){return kTRUE; }
    if (dileptons.Pt() < 30000) {return kTRUE; }
    
    MeT.SetPtEtaPhiE(*met_et, 0, *met_phi, *met_et);
    
    float dPhi_ll = TMath::Abs(lep_phi[0] - lep_phi[1]);
    dPhi_ll = dPhi_ll < TMath::Pi() ? dPhi_ll : 2*TMath::Pi() - dPhi_ll;
    if (dPhi_ll > 1.8) {return kTRUE; }
    
    float dPhillmet = TMath::Abs( dileptons.Phi() - MeT.Phi() );
    dPhillmet = dPhillmet < TMath::Pi() ? dPhillmet : 2*TMath::Pi() - dPhillmet;
    if (dPhillmet < TMath::Pi()/2) {return kTRUE; }
    
    float mt    = sqrt(2*dileptons.Pt()*MeT.Et()*(1-cos(dileptons.DeltaPhi(MeT))))/1000.;
    
    //Preselection of good jets
    int goodjet_n = 0;
    int goodbjet_n = 0;

    int goodjet_index[*jet_n];
    int jet_index = 0;

    int goodbjet_index[*jet_n];
    int bjet_index = 0;
    for(unsigned int i = 0; i < *jet_n; i++)
    {
        if(jet_pt[i] > 20000. && TMath::Abs(jet_eta[i]) < 2.5)
        {
            // JVT cleaning
            bool jvt_pass=true;
            if (jet_pt[i] < 60000. && TMath::Abs(jet_eta[i]) < 2.4 && jet_jvt[i] < 0.59) jvt_pass=false;
            if (jvt_pass)
            {
                // cut on 85% WP
                if ( jet_MV2c10[i] > 0.1758475  && TMath::Abs(jet_eta[i]) < 2.5 )
                {
                    goodbjet_n++;
                    goodbjet_index[bjet_index] = i;
                    bjet_index++;
                }

                if (jet_pt[i] > 30000.)
                {
                    goodjet_n++;
                    goodjet_index[jet_index] = i;
                    jet_index++;
                }
            }
        }
    }
    
    if (goodjet_n > 1) {return kTRUE; }
    if (goodbjet_n > 0) {return kTRUE; }*/
    
    //Preselection cut for electron/muon trigger 
    if(*trigE || *trigM)
    {

        // Preselection of good leptons
        int goodlep_index[2];
        int goodlep_n = 0;
        int lep_index =0;

        for(unsigned int i=0; i<*lep_n; i++)
        {

            TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(lep_pt[i]/1000., lep_eta[i], lep_phi[i], lep_E[i]/1000.);

            // Lepton is Tight
            if( lep_isTightID->at(i) )
            {

                // standard lepton isolation requirement => strict isolation
                if( lep_pt[i] >15000. && ( (lep_ptcone30[i]/lep_pt[i]) < 0.1) && ( (lep_etcone20[i] / lep_pt[i]) < 0.1 ) )
                {
                    if ( lep_type[i]==11 && TMath::Abs(lep_eta[i]) < 2.47 && ( TMath::Abs(lep_eta[i]) < 1.37 || TMath::Abs(lep_eta[i]) > 1.52 ) ) 
                    {
                        if( TMath::Abs(lep_trackd0pvunbiased[i])/lep_tracksigd0pvunbiased[i] < 5 && TMath::Abs(lep_z0[i]*TMath::Sin(leptemp.Theta())) < 0.5) 
                        {
                            goodlep_n = goodlep_n + 1;
                            goodlep_index[lep_index] = i;
                            lep_index++;
                        }
                    }
                    // muon selection
                    if ( lep_type[i] ==13 && TMath::Abs(lep_eta[i]) < 2.5 ) 
                    {
                        if( TMath::Abs(lep_trackd0pvunbiased[i])/lep_tracksigd0pvunbiased[i] < 3 && TMath::Abs(lep_z0[i]*TMath::Sin(leptemp.Theta())) < 0.5) 
                        {
                            goodlep_n = goodlep_n + 1;
                            goodlep_index[lep_index] = i;
                            lep_index++;
                        }
                    }
                }
            }// tight
        }


        //Exactly two good leptons, leading lepton with pT > 22 GeV and the subleading lepton with pT > 15 GeV
        if(goodlep_n==2)
        {

            int goodlep1_index = goodlep_index[0];
            int goodlep2_index = goodlep_index[1];

            if(lep_pt[goodlep1_index] > 22000)
            {

                //two different-flavour opposite-sign leptons
                if ( lep_charge[goodlep1_index] * lep_charge[goodlep2_index]  < 0 ) 
                {
                    if ( lep_type[goodlep1_index] != lep_type[goodlep2_index] )
                    {

                        // TLorentzVector definitions
                        TLorentzVector Lepton_1  = TLorentzVector();
                        TLorentzVector Lepton_2  = TLorentzVector();
                        TLorentzVector      MeT  = TLorentzVector();

                        Lepton_1.SetPtEtaPhiE(lep_pt[goodlep1_index], lep_eta[goodlep1_index], lep_phi[goodlep1_index],lep_E[goodlep1_index]);
                        Lepton_2.SetPtEtaPhiE(lep_pt[goodlep2_index], lep_eta[goodlep2_index], lep_phi[goodlep2_index],lep_E[goodlep2_index]);
                        MeT.SetPtEtaPhiE(*met_et, 0, *met_phi , *met_et);

                        TLorentzVector     Lepton_12 = TLorentzVector();
                        Lepton_12 = Lepton_1 + Lepton_2;

                        float mLL       = Lepton_12.Mag()/1000.;
                        float ptLL      = Lepton_12.Pt()/1000.;

                        float dPhi_LL  = TMath::Abs(lep_phi[goodlep1_index] - lep_phi[goodlep2_index] );
                        dPhi_LL        = dPhi_LL < TMath::Pi() ? dPhi_LL : 2*TMath::Pi() - dPhi_LL;

                        Float_t MET = *met_et/1000.;

                        float dPhiLLmet = TMath::Abs( Lepton_12.Phi() - MeT.Phi() );
                        dPhiLLmet    = dPhiLLmet < TMath::Pi() ? dPhiLLmet : 2*TMath::Pi() - dPhiLLmet;

                        float mt    = sqrt(2*Lepton_12.Pt()*MeT.Et()*(1-cos(Lepton_12.DeltaPhi(MeT))))/1000.;


                        //Preselection of good jets
                        int goodjet_n = 0;
                        int goodbjet_n = 0;

                        int goodjet_index[*jet_n];
                        int jet_index = 0;

                        int goodbjet_index[*jet_n];
                        int bjet_index = 0;

                        for(unsigned int i=0; i<*jet_n; i++)
                        {
                            if(jet_pt[i] > 20000. && TMath::Abs(jet_eta[i]) < 2.5)
                            {
                                // JVT cleaning
                                bool jvt_pass=true;
                                if (jet_pt[i] < 60000. && TMath::Abs(jet_eta[i]) < 2.4 && jet_jvt[i] < 0.59) jvt_pass=false;
                                if (jvt_pass)
                                {

                                    // cut on 85% WP
                                    if ( jet_MV2c10[i] > 0.1758475  && TMath::Abs(jet_eta[i]) < 2.5 )
                                    {
                                        goodbjet_n++;
                                        goodbjet_index[bjet_index] = i;
                                        bjet_index++;
                                    }

                                    if (jet_pt[i]>30000.)
                                    {
                                        goodjet_n++;
                                        goodjet_index[jet_index] = i;
                                        jet_index++;
                                    }

                                }
                            }
                        }

                        //  remove low mass meson resonances and DY events; ggF regions, at least 1 jet
                        if ( mLL > 10 && goodjet_n <= 1 && MET > 20)
                        {
                            if ( dPhiLLmet > TMath::Pi()/2 )
                            {

                                if ( ptLL > 30 )
                                {

                                    if ( mLL < 55 )
                                    {

                                        if ( dPhi_LL < 1.8 ) 
                                        {      

                                            if ( goodbjet_n ==0 ) 
                                            {

                                                events_eu++;
                                                histograms["h_mll"][type]->Fill(dileptons.M()/1000., wgt);
                                                histograms["h_mt"][type]->Fill(mt, wgt);
                                                histograms["h_met"][type]->Fill(*met_et/1000., wgt);
                                                
                                                /*if (DSID != prev_DSID)
                                                {
                                                    cout << type << " ; " << DSID << endl;
                                                    prev_DSID = DSID;
                                                }*/

                                            }
                                        }
                                    }
                                } // selection      
                            } // jet cut
                        }
                    }
                }
            }
        }
    }
    
    // Fill histograms
    
    /*histograms["h_mll"][type]->Fill(dileptons.M()/1000., wgt);
    histograms["h_mt"][type]->Fill(mt, wgt);
    histograms["h_met"][type]->Fill(*met_et/1000., wgt);*/

    return kTRUE;
}

void HWWAnalysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void HWWAnalysis::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    cout << "Total number of processed events: " << events_tot << endl; 
    cout << "Number of events in eu channel: " << events_eu << endl; 

}

void HWWAnalysis::CreateHistograms()
{
    for (const TString &type : categories)
    {
        histograms["h_mll"][type] = new TH1D("h_mll", "Invariant mass", 30, 10, 55);
        histograms["h_mt"][type] = new TH1D("h_mt", "Transverse mass", 15, 50, 200);
        histograms["h_met"][type] = new TH1D("h_met", "Missing transverse mass", 30, 0, 150);
    }
    histograms["h_mll"]["data"] = new TH1D("h_mll_d", "Invariant mass", 30, 10, 55);
    histograms["h_mt"]["data"] = new TH1D("h_mt_d", "Transverse mass", 15, 50, 200);
    histograms["h_met"]["data"] = new TH1D("h_met_d", "Missing transverse energu", 30, 0, 150);
}

// A whole load of switch statements
// Basically takes a sample id as input and return the corresponding category at basically no cost
// other than my sanity writing down all the cases
TString HWWAnalysis::GetType(Int_t dsid)
{
    switch(dsid)
    {
        case 363356:
            return "Diboson";
        case 363358:
            return "Diboson";
        case 363359:
            return "Diboson";
        case 363360:
            return "Diboson";
        case 363489:
            return "Diboson";
        case 363490:
            return "Diboson";
        case 363491:
            return "Diboson";
        case 363492:
            return "Diboson";
        case 363493:
            return "Diboson";
        case 364250:
            return "Diboson";
         
        case 341081:
            return "Higgs";
        case 343981:
            return "Higgs";
        case 341122:
            return "Higgs";
        case 341155:
            return "Higgs";
        case 341947:
            return "Higgs";
        case 341964:
            return "Higgs";
        case 344235:
            return "Higgs";
        case 345060:
            return "Higgs";
        case 345323:
            return "Higgs";
        case 345324:
            return "Higgs";
        case 345325:
            return "Higgs";
        case 345327:
            return "Higgs";
        case 345336:
            return "Higgs";
        case 345337:
            return "Higgs";
        case 345445:
            return "Higgs";
        case 345041:
            return "Higgs";
        case 345318:
            return "Higgs";
        case 345319:
            return "Higgs";
                
        case 361100:
            return "Wjetsincl";
        case 361101:
            return "Wjetsincl";
        case 361102:
            return "Wjetsincl";
        case 361103:
            return "Wjetsincl";
        case 361104:
            return "Wjetsincl";
        case 361105:
            return "Wjetsincl";
            
        case 361106:
            return "Zjetsincl";
        case 361107:
            return "Zjetsincl";
        case 361108:
            return "Zjetsincl";

        case 364100:
            return "Zjets";
        case 364101:
            return "Zjets";
        case 364102:
            return "Zjets";
        case 364103:
            return "Zjets";
        case 364104:
            return "Zjets";
        case 364105:
            return "Zjets";
        case 364106:
            return "Zjets";
        case 364107:
            return "Zjets";
        case 364108:
            return "Zjets";
        case 364109:
            return "Zjets";
        case 364110:
            return "Zjets";
        case 364111:
            return "Zjets";
        case 364112:
            return "Zjets";
        case 364113:
            return "Zjets";
        case 364114:
            return "Zjets";
        case 364115:
            return "Zjets";
        case 364116:
            return "Zjets";
        case 364117:
            return "Zjets";
        case 364118:
            return "Zjets";
        case 364119:
            return "Zjets";
        case 364120:
            return "Zjets";
        case 364121:
            return "Zjets";
        case 364122:
            return "Zjets";
        case 364123:
            return "Zjets";
        case 364124:
            return "Zjets";
        case 364125:
            return "Zjets";
        case 364126:
            return "Zjets";
        case 364127:
            return "Zjets";
        case 364128:
            return "Zjets";
        case 364129:
            return "Zjets";
        case 364130:
            return "Zjets";
        case 364131:
            return "Zjets";
        case 364132:
            return "Zjets";
        case 364133:
            return "Zjets";
        case 364134:
            return "Zjets";
        case 364135:
            return "Zjets";
        case 364136:
            return "Zjets";
        case 364137:
            return "Zjets";
        case 364138:
            return "Zjets";
        case 364139:
            return "Zjets";
        case 364140:
            return "Zjets";
        case 364141:
            return "Zjets";  
            
        case 364156:
            return "Wjets";
        case 364157:
            return "Wjets";
        case 364158:
            return "Wjets";
        case 364159:
            return "Wjets";
        case 364160:
            return "Wjets";
        case 364161:
            return "Wjets";
        case 364162:
            return "Wjets";
        case 364163:
            return "Wjets";
        case 364164:
            return "Wjets";
        case 364165:
            return "Wjets";
        case 364166:
            return "Wjets";
        case 364167:
            return "Wjets";
        case 364168:
            return "Wjets";
        case 364169:
            return "Wjets";
        case 364170:
            return "Wjets";
        case 364171:
            return "Wjets";
        case 364172:
            return "Wjets";
        case 364173:
            return "Wjets";
        case 364174:
            return "Wjets";
        case 364175:
            return "Wjets";
        case 364176:
            return "Wjets";
        case 364177:
            return "Wjets";
        case 364178:
            return "Wjets";
        case 364179:
            return "Wjets";
        case 364180:
            return "Wjets";
        case 364181:
            return "Wjets";
        case 364182:
            return "Wjets";
        case 364183:
            return "Wjets";
        case 364184:
            return "Wjets";
        case 364185:
            return "Wjets";
        case 364186:
            return "Wjets";
        case 364187:
            return "Wjets";
        case 364188:
            return "Wjets";
        case 364189:
            return "Wjets";
        case 364190:
            return "Wjets";
        case 364191:
            return "Wjets";
        case 364192:
            return "Wjets";
        case 364193:
            return "Wjets";
        case 364194:
            return "Wjets";
        case 364195:
            return "Wjets";
        case 364196:
            return "Wjets";
        case 364197:
            return "Wjets";
            
        case 410011:
            return "singleTop";
        case 410012:
            return "singleTop";
        case 410013:
            return "singleTop";
        case 410014:
            return "singleTop";
        case 410025:
            return "singleTop";
        case 410026:
            return "singleTop";
                
        case 410155:
            return "topX";
        case 410218:
            return "topX";
        case 410219:
            return "topX";
            
        case 410000:
            return "ttbar";
                
        default:
            return "data";
    }
}

