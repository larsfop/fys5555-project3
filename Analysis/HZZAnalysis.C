#define HZZAnalysis_cxx
// The class definition in HZZAnalysis.h has been generated automatically
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
// root> T->Process("HZZAnalysis.C")
// root> T->Process("HZZAnalysis.C","some options")
// root> T->Process("HZZAnalysis.C+")
//


#include "HZZAnalysis.h"
#include "lester_mt2_bisect.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <map> 
#include <math.h>
#include <stdio.h>
#include "TROOT.h"
#include "TObjString.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>

void HZZAnalysis::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).
    
    asymm_mt2_lester_bisect::disableCopyrightMessage();

    events_tot = 0;  
    events_ee = 0; 
    events_uu = 0;
    events_eu = 0;
    events_selected = 0;
    n_goodjets = 0;
    
    Nsig = 0;
    Nsig_tot = 0;

    DSID = 0; 
    prev_DSID = 0; 

    gROOT->SetBatch(kTRUE);

    option = GetOption(); // Print TString from Process(obj, str)
    
    signal = option;

    lumi = 10.06;

    categories = {"Diboson", "Higgs", "Wjetsincl", "Zjetsincl", "Zjets", "Wjets", "singleTop", "topX", "ttbar"};
    
    CreateHistograms();
    
    SetupTree();
    
}

void HZZAnalysis::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    option = GetOption();

}

Bool_t HZZAnalysis::Process(Long64_t entry)
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
    n_el = 0;
    n_mu = 0;

    fReader.SetLocalEntry(entry);

    // Count events 
    events_tot++;  	
    
    DSID = *channelNumber;
    
    type = GetType(DSID);
    
    isData = false;
    if (type == "data") isData = true;
    
    if (type == signal) Nsig_tot++;

    // Calculate event weight 
    scalef = (*scaleFactor_PILEUP)*(*scaleFactor_ELE)*(*scaleFactor_MUON)*(*scaleFactor_LepTRIGGER);
    
    wgt = 1.0; 
    if( !isData )
    {
        wgt = (*mcWeight)*scalef*((*XSection*lumi*1000.)/(*SumWeights)); 
    }

    //------------------------------//
    // Event selection & Kinematics //
    //------------------------------//

if(*trigE || *trigM)
{
    
    event_id = *eventNumber;

    // Preselection of good leptons
    int goodlep_index[*lep_n];
    int goodlep_n = 0;
    int lep_index =0;

    lep_masses.resize(0);
    for(unsigned int i=0; i<*lep_n; i++)
    {
        TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(lep_pt[i]/1000., lep_eta[i], lep_phi[i], lep_E[i]/1000.);

        // loosely isolated and very soft 
        if( lep_pt[i] > 5000. && TMath::Abs(lep_eta[i]) < 2.5 && ( (lep_ptcone30[i]/lep_pt[i]) < 0.3) && ( (lep_etcone20[i] / lep_pt[i]) < 0.3 ) ) 
        {
            // electron
            if ( lep_type[i] == 11 && lep_pt[i] > 7000. && TMath::Abs(lep_eta[i]) <2.47 ) 
            {
                if( TMath::Abs(lep_trackd0pvunbiased[i])/lep_tracksigd0pvunbiased[i] < 5 && TMath::Abs(lep_z0[i]*TMath::Sin(leptemp.Theta())) < 0.5) 
                {
                    goodlep_n = goodlep_n + 1;
                    goodlep_index[lep_index] = i;
                    lep_index++;
                    n_el++;
                    lep_masses.push_back(0.511);
                }
            }
            //muon
            if ( lep_type[i] == 13) 
            {
                if( TMath::Abs(lep_trackd0pvunbiased[i])/lep_tracksigd0pvunbiased[i] < 3 && TMath::Abs(lep_z0[i]*TMath::Sin(leptemp.Theta())) < 0.5) 
                {
                    goodlep_n = goodlep_n + 1;
                    goodlep_index[lep_index] = i;
                    lep_index++;
                    n_mu++;
                    lep_masses.push_back(105.66);
                }
            }
        }
    }
    


    //Exactly four good leptons
    if(goodlep_n == 4 )
    {

        goodlep1_index = goodlep_index[0];
        goodlep2_index = goodlep_index[1];
        goodlep3_index = goodlep_index[2];
        goodlep4_index = goodlep_index[3];
        
        lep1_mass = lep_masses[0];
        lep2_mass = lep_masses[1];
        lep3_mass = lep_masses[2];
        lep4_mass = lep_masses[3];
        
        lep1_pt = lep_pt[goodlep1_index];
        lep2_pt = lep_pt[goodlep2_index];
        lep3_pt = lep_pt[goodlep3_index];
        lep4_pt = lep_pt[goodlep4_index];
        
        lep1_eta = lep_eta[goodlep1_index];
        lep2_eta = lep_eta[goodlep2_index];
        lep3_eta = lep_eta[goodlep3_index];
        lep4_eta = lep_eta[goodlep4_index];

        lep1_phi = lep_phi[goodlep1_index];
        lep2_phi = lep_phi[goodlep2_index];
        lep3_phi = lep_phi[goodlep3_index];
        lep4_phi = lep_phi[goodlep4_index];
        
        ev_met_pt = *met_et;
        ev_met_phi = *met_phi;
        
        //first lepton pT > 25 GeV, second > 15 GeV and third > 10 GeV		      
        if (lep_pt[goodlep1_index] > 25000. && lep_pt[goodlep2_index] > 15000. && lep_pt[goodlep3_index] > 10000. ) 
        {
            
            if (*met_et < 40000)
            {


                // TLorentzVector definitions
                Lepton_1.SetPtEtaPhiE(lep_pt[goodlep1_index], lep_eta[goodlep1_index], lep_phi[goodlep1_index],lep_E[goodlep1_index]);
                Lepton_2.SetPtEtaPhiE(lep_pt[goodlep2_index], lep_eta[goodlep2_index], lep_phi[goodlep2_index],lep_E[goodlep2_index]);
                Lepton_3.SetPtEtaPhiE(lep_pt[goodlep3_index], lep_eta[goodlep3_index], lep_phi[goodlep3_index],lep_E[goodlep3_index]);
                Lepton_4.SetPtEtaPhiE(lep_pt[goodlep4_index], lep_eta[goodlep4_index], lep_phi[goodlep4_index],lep_E[goodlep4_index]);

                MeT.SetPtEtaPhiE(*met_et, 0, *met_phi , *met_et);

                // minimisation of difference from the Z mass
                float delta_Z1=0; 
                float delta_Z2=0; 
                InvMassZ1=0; 
                float InvMassZ2=0;
                float delta_Z1_1=0; float delta_Z1_2=0; float delta_Z1_3=0;
                float delta_Z2_1=0; float delta_Z2_2=0; float delta_Z2_3=0;
                float InvMassZ1_1=0; float InvMassZ1_2=0; float InvMassZ1_3=0;
                float InvMassZ2_1=0; float InvMassZ2_2=0; float InvMassZ2_3=0;
                float sum_ZZ1=0; float sum_ZZ2=0; float sum_ZZ3=0;

                // final values
                float InvMassZ1_min=0; float InvMassZ2_min=0; float sum_ZZ_fin=0;

                lep1_ch = lep_charge[goodlep1_index];
                lep2_ch = lep_charge[goodlep2_index];
                lep3_ch = lep_charge[goodlep3_index];
                lep4_ch = lep_charge[goodlep4_index];

                float sum_charges = lep1_ch + lep2_ch + lep3_ch + lep4_ch;

                // step-by-step
                // opposite charge leptons
                if ( sum_charges == 0  ) 
                {

                    lep1_type = lep_type[goodlep1_index];
                    lep2_type = lep_type[goodlep2_index];
                    lep3_type = lep_type[goodlep3_index];
                    lep4_type = lep_type[goodlep4_index];

                    int sum_types  = lep1_type + lep2_type + lep3_type + lep4_type;

                    // type e=11, mu=13
                    // begin case e+e-e+e- or mu+mu-mu+mu-
                    if ( sum_types == 44 || sum_types == 52  )
                    {
                        if ( lep_type[goodlep1_index] == lep_type[goodlep2_index] && ( (lep_charge[goodlep1_index] * lep_charge[goodlep2_index]) < 0 )  )
                        {
                            InvMassZ1_1=(Lepton_1+Lepton_2).Mag()/1000.;
                            InvMassZ2_1=(Lepton_3+Lepton_4).Mag()/1000.;
                            delta_Z1_1 =  TMath::Abs(InvMassZ1_1 - 91.18); 
                            delta_Z2_1 =  TMath::Abs(InvMassZ2_1 - 91.18);
                        }
                        if ( lep_type[goodlep1_index] == lep_type[goodlep3_index]  && ( (lep_charge[goodlep1_index] * lep_charge[goodlep3_index]) < 0 ) )
                        {
                            InvMassZ1_2=(Lepton_1+Lepton_3).Mag()/1000.;
                            InvMassZ2_2=(Lepton_2+Lepton_4).Mag()/1000.;
                            delta_Z1_2 =  TMath::Abs(InvMassZ1_2 - 91.18); 
                            delta_Z2_2 =  TMath::Abs(InvMassZ2_2 - 91.18);
                        }
                        if ( lep_type[goodlep1_index] == lep_type[goodlep4_index]  && ( (lep_charge[goodlep1_index] * lep_charge[goodlep4_index]) < 0 ) )
                        {
                            InvMassZ1_3=(Lepton_1+Lepton_4).Mag()/1000.;
                            InvMassZ2_3=(Lepton_2+Lepton_3).Mag()/1000.;
                            delta_Z1_3 =  TMath::Abs(InvMassZ1_3 - 91.18); 
                            delta_Z2_3 =  TMath::Abs(InvMassZ2_3 - 91.18);
                        }

                        if(delta_Z1_1 < delta_Z2_1) { InvMassZ1_min = InvMassZ1_1; InvMassZ2_min = InvMassZ2_1;}
                        if(delta_Z2_1 < delta_Z1_1) { InvMassZ1_min = InvMassZ2_1; InvMassZ2_min = InvMassZ1_1;}

                        if(delta_Z1_2 < delta_Z2_2) { InvMassZ1_min = InvMassZ1_2; InvMassZ2_min = InvMassZ2_2;}
                        if(delta_Z2_2 < delta_Z1_2) { InvMassZ1_min = InvMassZ2_2; InvMassZ2_min = InvMassZ1_2;}

                        if(delta_Z1_3 < delta_Z2_3) { InvMassZ1_min = InvMassZ1_3; InvMassZ2_min = InvMassZ2_3;}
                        if(delta_Z2_3 < delta_Z1_3) { InvMassZ1_min = InvMassZ2_3; InvMassZ2_min = InvMassZ1_3;}

                    } // cases of eeee or mumumumu

                    ////////////////////////////////////
                    // case eemumu 
                    if ( sum_types == 48 )
                    {

                        if ( lep_type[goodlep1_index] == lep_type[goodlep2_index]  && ( (lep_charge[goodlep1_index] * lep_charge[goodlep2_index]) < 0 ) )
                        {
                            InvMassZ1=(Lepton_1+Lepton_2).Mag()/1000.;
                            InvMassZ2=(Lepton_3+Lepton_4).Mag()/1000.;
                            delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
                            delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
                        }
                        if ( lep_type[goodlep1_index] == lep_type[goodlep3_index]  && ( (lep_charge[goodlep1_index] * lep_charge[goodlep3_index]) < 0 ) )
                        {
                            InvMassZ1=(Lepton_1+Lepton_3).Mag()/1000.;
                            InvMassZ2=(Lepton_2+Lepton_4).Mag()/1000.;
                            delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
                            delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
                        }
                        if ( lep_type[goodlep1_index] == lep_type[goodlep4_index]  && ( (lep_charge[goodlep1_index] * lep_charge[goodlep4_index]) < 0 ) )
                        {
                            InvMassZ1=(Lepton_1+Lepton_4).Mag()/1000.;
                            InvMassZ2=(Lepton_2+Lepton_3).Mag()/1000.;
                            delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
                            delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
                        }

                        if(delta_Z1 < delta_Z2) { InvMassZ1_min = InvMassZ1; InvMassZ2_min = InvMassZ2;}
                        if(delta_Z2 < delta_Z1) { InvMassZ1_min = InvMassZ2; InvMassZ2_min = InvMassZ1;}
                    } // eemumu overe

                    if ( (sum_types == 44 || sum_types == 52 || sum_types == 48) )
                    {

                        TLorentzVector FourLepSystem = TLorentzVector();
                        FourLepSystem = Lepton_1 + Lepton_2 + Lepton_3 + Lepton_4;
                        float FourLepSystem_M = FourLepSystem.Mag()/1000.;
                        float FourLepSystem_pt = FourLepSystem.Pt()/1000.;
                        float FourLepSystem_y = FourLepSystem.Rapidity();

                        ev_mllll = FourLepSystem.M();
                        ev_mll = (Lepton_1 + Lepton_2).M();

                        //ev_mt2 = (Lepton_1 + Lepton_2 + MeT).Mt2();
                        ev_mt2 = asymm_mt2_lester_bisect::get_mT2(lep1_mass, lep1_pt*cos(lep1_phi), lep1_pt*sin(lep1_phi),
                                                                   lep2_mass, lep2_pt*cos(lep2_phi), lep2_pt*sin(lep2_phi),
                                                                   *met_et*cos(*met_phi), *met_et*sin(*met_phi),
                                                                   0, 0,
                                                                   0);


                        //Preselection of good jets
                        goodjet_n = 0;
                        goodjet_index.resize(0); // Forces vector to have zero entries

                        if (*jet_n > 0)
                        {
                            for(unsigned int i=0; i<*jet_n; i++)
                            {
                                if(jet_pt[i] > 30000. && TMath::Abs(jet_eta[i]) < 4.4)
                                {
                                    goodjet_n++;
                                    goodjet_index.push_back(i);
                                }
                            }
                        }
                        
                        n_goodjets += goodjet_n;
                        vector<Double_t> njet_pt (4,-999.0);
                        vector<Double_t> njet_eta (4,-999.0);
                        vector<Double_t> njet_phi (4,-999.0);
                        for (Int_t i = 0; i < goodjet_n; i++)
                        {
                            njet_pt.push_back(jet_pt[goodjet_index[i]]);
                            njet_eta.push_back(jet_eta[goodjet_index[i]]);
                            njet_phi.push_back(jet_phi[goodjet_index[i]]);
                        }
                        
                        if (*jet_n < 2)
                        {
                            jet1_pt = njet_pt[0];
                            jet1_eta = njet_eta[0];
                            jet1_phi = njet_phi[0];
                        }
                        
                        if (*jet_n < 3)
                        {
                            jet2_pt = njet_pt[1];
                            jet2_eta = njet_eta[1];
                            jet2_phi = njet_phi[1];
                        }
                        
                        if (*jet_n < 4)
                        {
                            jet3_pt = njet_pt[2];
                            jet3_eta = njet_eta[2];
                            jet3_phi = njet_phi[2];
                        }
                        
                        if (*jet_n < 5)
                        {
                            jet4_pt = njet_pt[3];
                            jet4_eta = njet_eta[3];
                            jet4_phi = njet_phi[3];
                        }
                                
                        events_selected++;
                        if (type == signal) Nsig++;

                        //Start to fill histograms : definitions of x-axis variables
                        histograms["h_mllll"][type]->Fill(ev_mllll/1000., wgt);

                        if (n_el == 4) {histograms["h_meeee"][type]->Fill(ev_mllll/1000., wgt); }
                        if (n_mu == 4) {histograms["h_muuuu"][type]->Fill(ev_mllll/1000., wgt); }
                        if (n_el == 2) {histograms["h_meeuu"][type]->Fill(ev_mllll/1000., wgt); }

                        // write down a tree filled with features for ML
                        if (!isData)
                        {
                            // signal = 1; Background = 0
                            if (type == signal) label = 1;
                            else label = 0;

                            dataframe->Fill();
                        }

                    }
                }
            }
        }
    }
}
    
    return kTRUE;
}

void HZZAnalysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void HZZAnalysis::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    cout << "Total number of processed events: " << events_tot << endl; 
    cout << "Number of selected events: " << events_selected << endl; 
    cout << "Number of jets: " << n_goodjets << endl;

}

void HZZAnalysis::CreateHistograms()
{
    for (const TString &type : categories)
    {
        histograms["h_mllll"][type] = new TH1D("h_mllll", "4L invariant mass", 24, 80, 170);
        histograms["h_meeee"][type] = new TH1D("h_meeee", "4 electron invariant mass", 24, 80, 170);
        histograms["h_muuuu"][type] = new TH1D("h_muuuu", "4 muon invariant mass", 24, 80, 170);
        histograms["h_meeuu"][type] = new TH1D("h_meeuu", "2 electron and 2 muon invariant mass", 24, 80, 170);
    }
    histograms["h_mllll"]["data"] = new TH1D("h_mllll_d", "4L invariant mass", 24, 80, 170);
    histograms["h_meeee"]["data"] = new TH1D("h_meeee_d", "4 electron invariant mass", 24, 80, 170);
    histograms["h_muuuu"]["data"] = new TH1D("h_muuuu_d", "4 muon invariant mass", 24, 80, 170);
    histograms["h_meeuu"]["data"] = new TH1D("h_meeuu_d", "2 electron and 2 muon invariant mass", 24, 80, 170);
}

void HZZAnalysis::SetupTree()
{
    
    // dataframe->Branch("name", var, "type");
    dataframe->Branch("eventID", &event_id, "eventID/I");
    dataframe->Branch("n_el", &n_el, "n_el/I");
    dataframe->Branch("n_mu", &n_mu, "n_mu/I");
    dataframe->Branch("n_jet", &goodjet_n, "n_jet/I");
    dataframe->Branch("weight", &wgt, "weight/D");
    dataframe->Branch("lep1_pt", &lep1_pt, "lep1_pt/D");
    dataframe->Branch("lep2_pt", &lep2_pt, "lep2_pt/D");
    dataframe->Branch("lep3_pt", &lep3_pt, "lep3_pt/D");
    dataframe->Branch("lep4_pt", &lep4_pt, "lep4_pt/D");
    dataframe->Branch("lep1_eta", &lep1_eta, "lep1_eta/D");
    dataframe->Branch("lep2_eta", &lep2_eta, "lep2_eta/D");
    dataframe->Branch("lep3_eta", &lep3_eta, "lep3_eta/D");
    dataframe->Branch("lep4_eta", &lep4_eta, "lep4_eta/D");
    dataframe->Branch("lep1_phi", &lep1_phi, "lep1_phi/D");
    dataframe->Branch("lep2_phi", &lep2_phi, "lep2_phi/D");
    dataframe->Branch("lep3_phi", &lep3_phi, "lep3_phi/D");
    dataframe->Branch("lep4_phi", &lep4_phi, "lep4_phi/D");
    dataframe->Branch("lep1_ch", &lep1_ch, "lep1_ch/I");
    dataframe->Branch("lep2_ch", &lep2_ch, "lep2_ch/I");
    dataframe->Branch("lep3_ch", &lep3_ch, "lep3_ch/I");
    dataframe->Branch("lep4_ch", &lep4_ch, "lep4_ch/I");
    dataframe->Branch("lep1_type", &lep1_type, "lep1_type/I");
    dataframe->Branch("lep2_type", &lep2_type, "lep2_type/I");
    dataframe->Branch("lep3_type", &lep3_type, "lep3_type/I");
    dataframe->Branch("lep4_type", &lep4_type, "lep4_type/I");
    dataframe->Branch("lep1_mass", &lep1_mass, "lep1_mass/D");
    dataframe->Branch("lep2_mass", &lep2_mass, "lep2_mass/D");
    dataframe->Branch("lep3_mass", &lep3_mass, "lep3_mass/D");
    dataframe->Branch("lep4_mass", &lep4_mass, "lep4_mass/D");
    dataframe->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
    dataframe->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
    dataframe->Branch("jet3_pt", &jet3_pt, "jet3_pt/D");
    dataframe->Branch("jet4_pt", &jet4_pt, "jet4_pt/D");
    dataframe->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
    dataframe->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
    dataframe->Branch("jet3_eta", &jet3_eta, "jet3_eta/D");
    dataframe->Branch("jet4_eta", &jet4_eta, "jet4_eta/D");
    dataframe->Branch("jet1_phi", &jet1_phi, "jet1_phi/D");
    dataframe->Branch("jet2_phi", &jet2_phi, "jet2_phi/D");
    dataframe->Branch("jet3_phi", &jet3_phi, "jet3_phi/D");
    dataframe->Branch("jet4_phi", &jet4_phi, "jet4_phi/D");
    dataframe->Branch("ev_mll", &ev_mll, "ev_mll/D");
    dataframe->Branch("ev_mllll", &ev_mllll, "ev_mllll/D");
    dataframe->Branch("met_pt", &ev_met_pt, "met_pt/D");
    dataframe->Branch("met_phi", &ev_met_phi, "met_phi/D");
    dataframe->Branch("ev_mt2", &ev_mt2, "ev_mt2/D");
    dataframe->Branch("label", &label, "label/I");
    
}

// A whole load of switch statements
// Basically takes a sample id as input and return the corresponding category at basically no cost
// other than my sanity writing down all the cases
TString HZZAnalysis::GetType(Int_t dsid)
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

