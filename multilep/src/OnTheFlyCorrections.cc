#include "heavyNeutrino/multilep/interface/OnTheFlyCorrections.h"

OnTheFlyCorrections::OnTheFlyCorrections(std::string path, std::string gt, bool isdata){
        std::vector<std::string> runs;
        if(isdata) runs = {"BCDV4_DATA", "EFV4_DATA", "GV4_DATA", "HV4_DATA"};
        else       runs = {"V4_MC"};

       for(std::string run : runs){
	 	  jetUncertainties[run] = new JetCorrectionUncertainty(path+gt+run+"_Uncertainty_AK4PFchs.txt");
          std::vector<JetCorrectorParameters> jcParam;
          jcParam.push_back(JetCorrectorParameters(path+gt+run+"_L1FastJet_AK4PFchs.txt"));
          jcParam.push_back(JetCorrectorParameters(path+gt+run+"_L2Relative_AK4PFchs.txt"));
          jcParam.push_back(JetCorrectorParameters(path+gt+run+"_L3Absolute_AK4PFchs.txt"));
          if(isdata) jcParam.push_back(JetCorrectorParameters(path+gt+run+"_L2L3Residual_AK4PFchs.txt"));
	  	  jetCorrectors[run] = new FactorizedJetCorrector(jcParam);
        }
	fIsData = isdata;
}

std::string getRunName(unsigned long runNumber){
  if(runNumber < 1)           return "V4_MC";
  else if(runNumber < 276812) return "BCDV4_DATA";
  else if(runNumber < 278809) return "EFV4_DATA";
  else if(runNumber < 280385) return "GV4_DATA";
  else                        return "HV4_DATA";

}


float OnTheFlyCorrections::getJECUncertainty(float pt, float eta, int runnumber){
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	jetUncertainties[fIsData ? getRunName(runnumber) : "V4_MC"]->setJetPt(pt);
	jetUncertainties[fIsData ? getRunName(runnumber) : "V4_MC"]->setJetEta(eta);
	return jetUncertainties[fIsData ? getRunName(runnumber) : "V4_MC"]->getUncertainty(true);
}

float OnTheFlyCorrections::getJetCorrection(float pt, float corr, float eta, float rho, float area, std::string level = "L3Absolute"){
	// sets the pT back to raw and returns the raw-pT correction factor
	return getJetCorrectionRawPt(pt/corr, eta, rho, area, level);
}

float OnTheFlyCorrections::getJetCorrectionRawPt(float rawpt, float eta, float rho, float area, std::string level = "L3Absolute", int runnumber){
	// slighly redundant considering we have what we have below, but I think that's what frederic was thinking about
        FactorizedJetCorrector* jetCorrector = jetCorrectors[fIsData ? getRunName(runnumber) : "V4_MC"];
	jetCorrector->setJetEta(eta);
	jetCorrector->setRho(rho);
	jetCorrector->setJetA(area);
	jetCorrector->setJetPt(rawpt); // new raw-pT here...! this is called with fTR->JPt[jetindex]/fTR->JEcorr[jetindex]; in the useranalysisbase.
	std::vector< float > corrections = jetCorrector->getSubCorrections();

	if (level == "L1FastJet")    return corrections.front();
	if (level == "L2Relative")   return corrections[1];
	if (level == "L2L3Residual") return corrections.back();
	return corrections[2];
}



std::pair<float,float> OnTheFlyCorrections::getCorrections(float rawpt, float raweta, float rawnomupt, float phi, float emf, float rho, float area, int runnumber) {
  std::pair<float, float> corr = std::pair<float, float>(0, 0);
  FactorizedJetCorrector* jetCorrector = jetCorrectors[fIsData ? getRunName(runnumber) : "V4_MC"];
  jetCorrector->setJetEta(raweta);
  jetCorrector->setRho(rho);
  jetCorrector->setJetA(area);
  jetCorrector->setJetPt(rawpt);
  
  std::vector< float > corrections = jetCorrector->getSubCorrections();
  

  float l1corrpt   = rawnomupt*corrections.front(); // l1fastjet corrections were pushed pack first
  float fullcorrpt = rawnomupt*corrections.back();  // full corrections are the last in the vector


  // the corretions for the MET are the difference between l1fastjet and the full corrections on the jet!
  if(emf > 0.9 or fullcorrpt < 15.) return corr;       // skip jets with EMF > 0.9
  
  corr.first  = getPx(l1corrpt - fullcorrpt, phi); // fill the px of the correction in the pair.first
  corr.second = getPy(l1corrpt - fullcorrpt, phi); // and the py in the pair.second
  
  return corr;
}
