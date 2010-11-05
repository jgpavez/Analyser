#ifndef TIDENTIFICATOR_H
#define TIDENTIFICATOR_H

#include "TClasTool.h"


class TIdentificator {
public:
    explicit TIdentificator(TClasTool* CT = 0);
    ~TIdentificator();

    // HEADER bank
    Float_t NEvent();                             // inline

    // EVNT bank
    Double_t Betta(Int_t k);                      // inline
    Double_t Id(Int_t k, Bool_t kind = 0);        // inline
    Double_t Charge(Int_t k);                     // inline
    Double_t Px(Int_t k, Bool_t kind = 0);        // inline
    Double_t Py(Int_t k, Bool_t kind = 0);        // inline
    Double_t Pz(Int_t k, Bool_t kind = 0);        // inline
    Double_t X(Int_t k);                          // inline
    Double_t Y(Int_t k);                          // inline
    Double_t Z(Int_t k);                          // inline
    Int_t StatCC(Int_t k);                        // inline
    Int_t StatSC(Int_t k);                        // inline
    Int_t StatDC(Int_t k);                        // inline
    Int_t StatEC(Int_t k);                        // inline
    Double_t Status(Int_t k);                     // inline

    // CCPB
    Double_t Nphe(Int_t k);                       // inline
    Double_t Chi2CC(Int_t k);                     // inline
    Double_t CCStatus(Int_t k);                   // inline

    // DCPB
    Double_t DCStatus(Int_t k);                   // inline

    // ECPB
    Double_t Etot(Int_t k);                       // inline
    Double_t Ein(Int_t k);                        // inline
    Double_t Eout(Int_t k);                       // inline
    Double_t ECStatus(Int_t k);                   // inline

    // SCPB
    Double_t PathSC(Int_t k);                     // inline
    Double_t TimeSC(Int_t k);                     // inline
    Double_t EdepSC(Int_t k);                     // inline
    Double_t SCStatus(Int_t k);                   // inline

    // Derived observables
    Double_t Momentum(Int_t k, Bool_t kind = 0);
    Double_t Moment(Int_t k, Bool_t kind = 0);    // Deprecated
    Double_t Mass2(Int_t k);
    Double_t ThetaLab(Int_t k, Bool_t kind = 0);
    Double_t PhiLab(Int_t k, Bool_t kind = 0);
    Double_t ThetaVirtLab(Bool_t kind = 0);
    Double_t PhiVirtLab(Bool_t kind = 0);
    Double_t ThetaPQ(Int_t k, Bool_t kind = 0);
    Double_t PhiPQ(Int_t k, Bool_t kind = 0);
    Double_t CosThetaPQ(Int_t k, Bool_t kind = 0);
    Double_t PTrans2PQ(Int_t k, Bool_t kind = 0);
    Double_t PLong2PQ(Int_t k, Bool_t kind = 0);
    Int_t Sector(Int_t k, Bool_t kind = 0);

    // Kinematic variables
    Double_t Q2(Bool_t kind = 0);
    Double_t W(Bool_t kind = 0);
    Double_t Nu(Bool_t kind = 0);
    Double_t ZhPi(Int_t k, Double_t Mass, Bool_t kind = 0);

    // Correction functions
    Double_t TimeCorr4(Double_t mass, Int_t k);

    // Particle Identification cuts
    TString GetCategorization(Int_t k);
    TString* GetCategorization();
    void PrintCategorization();
    void PrintCategorization(TString* partIds);

    // Fiducial Cut
    Double_t FidTheta(Int_t k, Bool_t kind = 0);
    Double_t FidThetaMin();
    Double_t FidFunc(Int_t side, Int_t param);
    Double_t FidPhi(Int_t k, Bool_t kind = 0);
    Double_t FidPhiMin();
    Double_t FidPhiMax();
    Bool_t FidCheckCut();
    Int_t FidSector(Int_t k, Bool_t kind = 0);
    Int_t ElecVertTarg(Int_t k, Bool_t = 0);
    Bool_t PionVertTarg(Int_t k, Bool_t = 0);

private:
    const Double_t kEbeam;    // The energy of incoming electron beam
    const Double_t kMpi;      // The mass of the pion
    const Double_t kGOOD;     // The key for the exceptions (should be improved to avoid it at all !!!)

    TClasTool *fCT;           // Pointer to the main ClasTool object
    TEVNTClass *fEVNT;        // Pointer to the EVNT object
    TGSIMClass *fGSIM;        // Pointer to the GSIM object
    TCCPBClass *fCCPB;        // Pointer to the CCPB object
    TECPBClass *fECPB;        // Pointer to the ECPB object
    TSCPBClass *fSCPB;        // Pointer to the SCPB object
    TDCPBClass *fDCPB;        // Pointer to the DCPB object

    TString* fPartIds;        // Array with the categories of the particles belonging to an event.
};

#include "TIdentificator.icc"

#endif
