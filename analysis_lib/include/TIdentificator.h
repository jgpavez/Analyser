#ifndef TIDENTIFICATOR_H
#define TIDENTIFICATOR_H

#include "TClasTool.h"


class TIdentificator {
protected:
    const Double_t kEbeam;    // The energy of incoming electron beam
    const Double_t kMpi;      // The mass of the pion
    const Double_t kGOOD;     // The key for the exceptions (should be improved to avoid it at all !!!)

    TClasTool *fCT;
    TEVNTClass *fEVNT;
    TGSIMClass *fGSIM;
    TCCPBClass *fCCPB;
    TECPBClass *fECPB;
    TSCPBClass *fSCPB;
    TDCPBClass *fDCPB;

public:
    TIdentificator(TClasTool *);
    ~TIdentificator();

    // HEADER bank
    Float_t NEvent();                           //inline

    // EVNT bank
    Double_t Betta(Int_t);                      //inline
    Double_t Id(Int_t, Bool_t = 0);
    Double_t Charge(Int_t);                     //inline
    Double_t Px(Int_t, Bool_t = 0);
    Double_t Py(Int_t, Bool_t = 0);
    Double_t Pz(Int_t, Bool_t = 0);
    Double_t X(Int_t);                          //inline
    Double_t Y(Int_t);                          //inline
    Double_t Z(Int_t);                          //inline
    Int_t StatCC(Int_t);                        //inline - CCPB key
    Int_t StatSC(Int_t);                        //inline - SCPB key
    Int_t StatDC(Int_t);                        //inline - DCPB key
    Int_t StatEC(Int_t);                        //inline - ECPB key
    Double_t Status(Int_t);                     //inline - Event general status

    // CCPB bank
    Double_t Nphe(Int_t);
    Double_t Chi2CC(Int_t);
    Double_t CCStatus(Int_t);

    // DCPB bank
    Double_t DCStatus(Int_t);

    // ECPB bank
    Double_t Etot(Int_t k);
    Double_t Ein(Int_t k);
    Double_t Eout(Int_t k);
    Double_t ECStatus(Int_t);

    // SCPB bank
    Double_t PathSC(Int_t k);
    Double_t TimeSC(Int_t k);
    Double_t EdepSC(Int_t k);
    Double_t SCStatus(Int_t);

    // Derived observables
    Double_t Momentum(Int_t, Bool_t = 0);
    Double_t Moment(Int_t, Bool_t = 0);  // Deprecated
    Double_t Mass2(Int_t);
    Double_t ThetaLab(Int_t, Bool_t = 0);
    Double_t PhiLab(Int_t, Bool_t = 0);
    Double_t ThetaVirtLab(Bool_t = 0); // Check if it is correct !!!
    Double_t PhiVirtLab(Bool_t = 0); // Check if it is correct !!!
    Double_t ThetaPQ(Int_t, Bool_t = 0);
    Double_t PhiPQ(Int_t, Bool_t = 0);
    Double_t CosThetaPQ(Int_t, Bool_t = 0);
    Double_t PTrans2PQ(Int_t, Bool_t = 0);
    Double_t PLong2PQ(Int_t, Bool_t = 0);
    Int_t Sector(Int_t, Bool_t = 0); // Check if it is correct !!! Add k==1

    // Kinematic variables
    Double_t Q2(Bool_t = 0);
    Double_t W(Bool_t = 0);
    Double_t Nu(Bool_t = 0);
    Double_t ZhPi(Int_t, Bool_t = 0, Double_t = 0.139);

    // Correction functions
    Double_t TimeCorr4(Double_t, Int_t);

    // Particle Identification cuts
    TString GetCategorization(Int_t);
    TString* GetCategorization();
    void PrintCategorization();
    void PrintCategorization(TString*);

    // Fiducial Cut
    Double_t FidTheta(Int_t, Bool_t = 0);
    Double_t FidThetaMin();
    Double_t FidFunc(Int_t, Int_t);
    Double_t FidPhi(Int_t, Bool_t = 0);
    Double_t FidPhiMin();
    Double_t FidPhiMax();
    Bool_t FidCheckCut();
    Int_t FidSector(Int_t, Bool_t = 0);

private:
    TString* fPartIds;
};

#include "TIdentificator.icc"

#endif
