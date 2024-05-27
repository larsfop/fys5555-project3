#define AnalysisSelector_cxx
//////////////////////////////////////////////////////////////////////////////////////////
// Creates a TChain to be used by the Analysis.C class
#include <stdio.h>

#include "AnalysisSelector.h"

using namespace std;

AnalysisSelector::AnalysisSelector(TChain *chain, TString Analysis)
{
    m_chain = chain;
    if (Analysis == "GamGam")
    {
        m_selector = TSelector::GetSelector("Analysis/HyyAnalysis.C+");
    }
    else if (Analysis == "2lep")
    {
        m_selector = TSelector::GetSelector("Analysis/HWWAnalysis.C+");
    }
    else if (Analysis == "4lep")
    {
        m_selector = TSelector::GetSelector("Analysis/HZZAnalysis.C+");
    }
}

TSelector *AnalysisSelector::Selector()
{
    return m_selector;
}

void AnalysisSelector::Process(TString option)
{
    printf("-------------------------------------------\n");
    printf("Processing MC and Data\n");
    printf("Number of events to process: %lld\n", m_chain->GetEntries());
    printf("-------------------------------------------\n");

    m_chain->Process(m_selector, option);
}

