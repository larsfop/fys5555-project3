#define HyyAnalysis_cxx
// The class definition in HyyAnalysis.h has been generated automatically
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
// root> T->Process("HyyAnalysis.C")
// root> T->Process("HyyAnalysis.C","some options")
// root> T->Process("HyyAnalysis.C+")
//


#include "HyyAnalysis.h"
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

void HyyAnalysis::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    events_tot = 0;  
    events_yy = 0;
        
    Nsig = 0;
    Nsig_tot = 0;

    DSID = 0; 

    gROOT->SetBatch(kTRUE);

    option = GetOption(); // Print TString from Process(obj, str)
    
    signal = option;

    lumi = 10.06;

    categories = {"Diboson", "Higgs", "Wjetsincl", "Zjetsincl", "Zjets", "Wjets", "singleTop", "topX", "ttbar"};
    
    CreateHistograms();

}

void HyyAnalysis::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    option = GetOption();

}

Bool_t HyyAnalysis::Process(Long64_t entry)
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
    
    if (type == signal) Nsig_tot++;

    // Calculate event weight 
    wgt = 1.0; 
    if( !isData )
    {
        wgt = (*mcWeight)*(*scaleFactor_PILEUP)*(*scaleFactor_PHOTON)*(*scaleFactor_PhotonTRIGGER)*((*XSection)*lumi*1000.)/(*SumWeights); 
        //wgt = (*mcWeight)*(*scaleFactor_PILEUP)*(*scaleFactor_PHOTON)*(*scaleFactor_PhotonTRIGGER);
    }

    //------------------------------//
    // Event selection & Kinematics //
    //------------------------------//
    
    if (*trigP)
    {
        
        int goodphoton_index[2];
        int goodphoton_n = 0;
        int photon_index =0;
        
        // Selection if good photons
        for (unsigned int i = 0; i < *photon_n; i++)
        {
            if (photon_isTightID->at(i))
            {
                if (photon_pt[i] > 25000 && TMath::Abs(photon_eta[i]) < 2.37 && (TMath::Abs(photon_eta[i]) < 1.37 || TMath::Abs(photon_eta[i]) > 1.52))
                {
                    goodphoton_n = goodphoton_n + 1;
                    goodphoton_index[photon_index] = i;
                    photon_index++;
                }
            }
        }
        
        if (goodphoton_n == 2)
        {
            int goodphoton1_index = goodphoton_index[0];
            int goodphoton2_index = goodphoton_index[1];

            if (photon_pt[goodphoton1_index] > 35000)
            {
                // isolated photons
                if( ( (photon_ptcone30[goodphoton1_index]/photon_pt[goodphoton1_index]) < 0.065) && ( (photon_etcone20[goodphoton1_index] / photon_pt[goodphoton1_index]) < 0.065 ) )
                {
                    if( ( (photon_ptcone30[goodphoton2_index]/photon_pt[goodphoton2_index]) < 0.065) && ( (photon_etcone20[goodphoton2_index] / photon_pt[goodphoton2_index]) < 0.065 ) )
                    {


                        // create 2 vectors   
                        TLorentzVector Photon_1  = TLorentzVector();
                        TLorentzVector Photon_2  = TLorentzVector();

                        Photon_1.SetPtEtaPhiE(photon_pt[goodphoton1_index], photon_eta[goodphoton1_index], photon_phi[goodphoton1_index],photon_E[goodphoton1_index]);
                        Photon_2.SetPtEtaPhiE(photon_pt[goodphoton2_index], photon_eta[goodphoton2_index], photon_phi[goodphoton2_index],photon_E[goodphoton2_index]);

                        // calculate dPhi(photon-photon)
                        float dPhi_yy = TMath::Abs(photon_phi[goodphoton1_index] - photon_phi[goodphoton2_index] );
                        dPhi_yy       = dPhi_yy < TMath::Pi() ? dPhi_yy : 2*TMath::Pi() - dPhi_yy;

                        // diphoton mass
                        float m_yy  = sqrt( 2 * Photon_1.Pt()/1000. * Photon_2.Pt()/1000. * (cosh( Photon_1.Eta() - Photon_2.Eta()) - cos(dPhi_yy)));
                        // kinematics
                        float Photon_1_kin = Photon_1.Pt()/1000. / m_yy;
                        float Photon_2_kin = Photon_2.Pt()/1000. / m_yy;

                        // kinematical selection 
                        if ( Photon_1_kin > 0.35 && Photon_2_kin > 0.25 ) 
                        { 

                            // mass-window cut
                            if(m_yy > 105 && m_yy < 160 ) 
                            {

                                events_yy++;
                                if (type == signal) Nsig++;
                                if (type != type_old)
                                histograms["h_myy"][type]->Fill(m_yy, wgt);

                            }
                        }
                    }
                }
            }
        }
    }
    
    // Fill histograms
    
    //histograms["h_myy"][type]->Fill(dileptons.M()/1000., wgt);

    return kTRUE;
}

void HyyAnalysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void HyyAnalysis::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
    
    // Signal plus background fit
    fit = new TF1("fit", "([0]+[1]*x+[2]*x^2+[3]*x^3)+[4]*exp(-0.5*((x-[5])/[6])^2)", 105, 160);
    fit->FixParameter(5, 125.0);
    fit->FixParameter(4, 119.1);
    fit->FixParameter(6, 2.39);
    
    // Fit to data
    histograms["h_myy"]["data"]->Fit("fit", "0", "", 105, 160);
    TF1 *fitresults = histograms["h_myy"]["data"]->GetFunction("fit");
    
    // Background fit
    bkg = new TF1("bkg", "([0]+[1]*x+[2]*x^2+[3]*x^3)", 105, 160);
    for (int i = 0; i < 4; i++)
    {
        bkg->SetParameter(i, fit->GetParameter(i));
    }

    cout << "Total number of processed events: " << events_tot << endl; 
    cout << "Number of events in yy channel: " << events_yy << endl; 

}

void HyyAnalysis::CreateHistograms()
{
    for (const TString &type : categories)
    {
        histograms["h_myy"][type] = new TH1D("h_myy", "Invariant mass", 24, 105, 160);
    }
    histograms["h_myy"]["data"] = new TH1D("h_myy_d", "Invariant mass", 24, 105, 160);
}

TF1* HyyAnalysis::GetFit(TString fit_name)
{
    if (fit_name == "bkg+sig")
    {
        return fit;
    }
    else if (fit_name == "bkg")
    {
        return bkg;
    }
}

// A whole load of switch statements
// Basically takes a sample id as input and return the corresponding category at basically no cost
// other than my sanity writing down all the cases
TString HyyAnalysis::GetType(Int_t dsid)
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

