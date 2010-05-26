
TString TIdentificator::GetCategorization(Int_t k)
{
    Int_t number_dc = fCT->GetNRows("DCPB");
    Int_t number_cc = fCT->GetNRows("CCPB");
    Int_t number_sc = fCT->GetNRows("SCPB");
    Int_t number_ec = fCT->GetNRows("ECPB");

    TString partId;

    partId = "not recognized";

    if (number_dc != 0) {
        if (k == 0 &&
                    Status(0) > 0 && Status(0) < 100 &&
                    Charge(0) == -1 &&
                    number_cc != 0 && number_ec != 0 && number_sc != 0 &&
                    StatCC(0) > 0 && StatSC(0) > 0 &&
                    StatDC(0) > 0 && StatEC(0) > 0 &&
                    DCStatus(0) > 0 && SCStatus(0) == 33 &&
                    Nphe(0) > 25 &&
                    Etot(0) / 0.27 + 0.4 > Moment(0) &&
                    Etot(0) / 0.27 - 0.4 < Moment(0) &&
                    Ein(0) + Eout(0) > 0.8 * 0.27 * Moment(0) &&
                    Ein(0) + Eout(0) < 1.2 * 0.27 * Moment(0) &&
                    FidCheckCut() == 1)
            partId = "electron";


        if (k > 0) {
            if (Charge(k)==0 && Betta(k)>0.95 && ECStatus(k)>0)
                partId = "photon";

            //positive particles
            if (Charge(k) == 1 &&
                        Status(k) > 0 && Status(k) < 100 &&
                        StatDC(k) > 0 && DCStatus(k) > 0) {
                if (Moment(k)>=2.7) {
                    if (number_cc != 0 && StatCC(k) > 0 &&
                                Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3)
                        partId = "high energy pion +";
                }

                if (Moment(k) < 2.7) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Moment(k) < 1 &&
                                    TimeCorr4(0.139,k) >= -1.46 &&
                                    TimeCorr4(0.139,k) <= 0.15) ||
                                (Moment(k) >=1 &&
                                    TimeCorr4(0.139,k) >= -1.38 &&
                                    TimeCorr4(0.139,k) <= 0.53)))
                        partId = "low energy pion +";
                }

                if (Moment(k) < 2.) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Moment(k) >= 1. &&
                                    TimeCorr4(0.938,k) >= -0.69 &&
                                    TimeCorr4(0.938,k) <= 1.38) ||
                                (Moment(k) < 1. &&
                                    TimeCorr4(0.938,k) >= -3.78 &&
                                    TimeCorr4(0.938,k) <= 6.75)))
                        partId = "low energy proton";
                }

                if (Charge(k) == 1 && number_cc != 0 &&
                            number_ec != 0 && number_sc != 0 &&
                            StatCC(k) > 0 && StatSC(k) > 0 &&
                            StatDC(k) > 0 && StatEC(k) > 0 &&
                            DCStatus(k) > 0 && Nphe(k) > 25 &&
                            Etot(k) / 0.27 + 0.4 > Moment(k) &&
                            Etot(k) / 0.27 - 0.4 < Moment(k) &&
                            Ein(k) + Eout(k) > 0.8 * 0.27 * Moment(k) &&
                            Ein(k) + Eout(k) < 1.2 * 0.27 * Moment(k))
                    partId = "positron";
            }
        }
    }

    return partId;
}
