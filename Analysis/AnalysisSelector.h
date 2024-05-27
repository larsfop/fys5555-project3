#ifndef RunSelector_h
#define RunSelector_h

#include <map>
#include <memory>
#include <iostream>
#include <string>
#include <any>

#include <TSelector.h>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TProof.h"

//void RunSelector(std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> x);

class AnalysisSelector
{
public:
    AnalysisSelector(TChain *chain, TString Analysis);
    TSelector *Selector();
    void Process(TString option);
    
private:
    TSelector *m_selector;
    TChain *m_chain;
};

#endif