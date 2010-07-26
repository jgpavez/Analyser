#include <string.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "TIdentificator.h"
#include "massConst.h"

using namespace std;

const Double_t kFidThetaMax = 54;

//For theta_min calculation for electron
const Double_t kThetaPar0[6] = {15        , 15        ,  15       , 15        ,  13       ,  13};
const Double_t kThetaPar1[6] = {-0.425145 , -1.02217  , -0.7837   , -1.47798  ,   3.47361 ,   3.5714};
const Double_t kThetaPar2[6] = {-0.666294 , -0.616567 , -0.673602 , -0.647113 ,  -0.34459 ,  -0.398458};
const Double_t kThetaPar3[6] = { 5.73077  ,  5.51799  ,  8.05224  ,  7.74737  ,   8.45226 ,   9.54265};
const Double_t kThetaPar4[6] = {10.4976   , 14.0557   , 15.2178   , 16.7291   , -63.4556  , -22.649};
const Double_t kThetaPar5[6] = {-1.13254  , -1.16189  , -2.08386  , -1.79939  ,  -3.3791  ,  -1.89746};

//For parameter 0 of the phi_min calculation for electron
const Double_t kFidPar0Low0[6] = {25     ,25       ,25       ,24.6345 ,23.4731 ,24.8599};
const Double_t kFidPar0Low1[6] = {-12    ,-12      ,-12      ,-12     ,-12     ,-12};
const Double_t kFidPar0Low2[6] = {0.5605 ,0.714261 ,0.616788 ,0.62982 ,1.84236 ,1.00513};
const Double_t kFidPar0Low3[6] = {4.4    ,4.4      ,4.4      ,4.4     ,4.4     ,4.4};

//For parameter 1 of the phi_min calculation for electron
const Double_t kFidPar1Low0[6] = {2.1945    ,4        ,3.3352  ,2.22769   ,1.63143   ,3.19807};
const Double_t kFidPar1Low1[6] = {1.51417   ,1.56882  ,2       ,2         ,1.90179   ,0.173168};
const Double_t kFidPar1Low2[6] = {-0.354081 ,-2       ,-2      ,-0.760895 ,-0.213751 ,-0.1};
const Double_t kFidPar1Low3[6] = {0.5       ,0.5      ,1.01681 ,1.31808   ,0.786844  ,1.6};

//For parameter 0 of the phi_max calculation for electron
const Double_t kFidPar0High0[6] = {25       ,25       ,25       ,25       ,23.7067 ,25};
const Double_t kFidPar0High1[6] = {-8       ,-10.3277 ,-12      ,-11.3361 ,-12     ,-11.4641};
const Double_t kFidPar0High2[6] = {0.479446 ,0.380908 ,0.675835 ,0.636018 ,2.92146 ,0.55553};
const Double_t kFidPar0High3[6] = {4.8      ,4.79964  ,4.4      ,4.4815   ,4.4     ,4.41327};

//For parameter 1 of the phi_max calculation for electron
const Double_t kFidPar1High0[6] = {3.57349  ,3.02279      ,2.02102    ,3.1948    ,3.0934    ,2.48828};
const Double_t kFidPar1High1[6] = {2        ,0.966175     ,2          ,0.192701  ,0.821726  ,2};
const Double_t kFidPar1High2[6] = {-2       ,-2           ,-1.70021   ,-1.27578  ,-0.233492 ,-2};
const Double_t kFidPar1High3[6] = {0.5      ,0.527823     ,0.68655    ,1.6       ,1.6       ,0.70261};


//ClassImp(TIdentificator);
	


TIdentificator::TIdentificator(TClasTool *CT)
    : kEbeam(5.014), kMpi(0.139570), kGOOD(-1000.)
{
    this->fCT = CT;
}



TIdentificator::~TIdentificator()
{
    fCT = 0;
}



Double_t TIdentificator::Moment(Int_t k, Bool_t kind)
{
  
    if (kind == 0) {
        return sqrt(pow(Px(k),2) + pow(Py(k),2) + pow(Pz(k),2));
    } else {                            // Fix this in case kind != 1
        return sqrt(pow(Px(k,1),2) + pow(Py(k,1),2) + pow(Pz(k,1),2));
    }
}



Double_t TIdentificator::Mass2(Int_t k)
{
      return pow(Moment(k), 2)* (pow(Betta(k), -2) - 1);
}



Double_t TIdentificator::ThetaLab(Int_t k, Bool_t kind)
{
    Double_t theta;

    if (kind == 0) {
        theta = acos(Pz(k) / Moment(k));
    } else {                            // Fix this in case kind != 1
        theta = acos(Pz(k,1) / Moment(k,1));
    }

    return theta;
}



Double_t TIdentificator::PhiLab(Int_t k, Bool_t kind)
{
    Double_t phi;

    if (kind == 0) {
        if (Py(k) >= 0 && Px(k) >= 0)     phi =   fabs(57.3 * atan(Px(k) / Py(k)));
        else if (Py(k) >= 0 && Px(k) < 0) phi = - fabs(57.3 * atan(Px(k) / Py(k)));
        else if (Py(k) < 0 && Px(k) >= 0) phi =  180 - fabs(57.3 * atan(Px(k) / Py(k)));
        else                              phi = -180 + fabs(57.3 * atan(Px(k) / Py(k)));
    } else {                            // Fix this in case kind != 1
        if (Py(k,1) >= 0 && Px(k,1) >= 0)     phi =   fabs(57.3 * atan(Px(k,1) / Py(k,1)));
        else if (Py(k,1) >= 0 && Px(k,1) < 0) phi = - fabs(57.3 * atan(Px(k,1) / Py(k,1)));
        else if (Py(k,1) < 0 && Px(k,1) >= 0) phi =  180 - fabs(57.3 * atan(Px(k,1) / Py(k,1)));
        else                                  phi = -180 + fabs(57.3 * atan(Px(k,1) / Py(k,1)));
    }

    return phi;
}



Double_t TIdentificator::ThetaVirtLab(Bool_t kind) // Check if it is correct !!!
{
    Double_t theta_virt;

    if (kind == 0) {
        theta_virt = acos((kEbeam - Moment(0) * cos(ThetaLab(0))) /
                           (sqrt(Q2() + Nu() * Nu())));
    } else {                            // Fix this in case kind != 1
        theta_virt = acos((kEbeam - Moment(0,1) * cos(ThetaLab(0,1))) /
                           (sqrt(Q2(1) + Nu(1) * Nu(1))));
    }

    return theta_virt;
}



Double_t TIdentificator::PhiVirtLab(Bool_t kind) // Check if it is correct !!!
{
    Double_t phi_virt;

    if (PhiLab(0,kind) > 0) phi_virt = - 180 + PhiLab(0,kind);
    else phi_virt = 180 + PhiLab(0,kind);

    return phi_virt;
}



Double_t TIdentificator::ThetaPQ(Int_t k, Bool_t kind)
{
    if (k == 0){
        cout <<"message HERE"<<endl;
    }
    if (kind == 0) {
        return acos(Pz(k) / Moment(k));
    } else {                            // Fix this in case kind != 1
        return acos(Pz(k,1) / Moment(k,1));
    }
}



Double_t TIdentificator::PhiPQ(Int_t k, Bool_t kind)
{
    Double_t phi_pq;
    if (k == 0){
        cout <<"message HERE"<<endl;
    }
    if (kind == 0) {
        TVector3 vpi(Px(k), Py(k), Pz(k));
        TVector3 vvirt(-Px(0), -Py(0), kEbeam - Pz(0));
        Double_t phi_z = TMath::Pi() - vvirt.Phi();
        vvirt.RotateZ(phi_z);
        vpi.RotateZ(phi_z);
            TVector3 vhelp(0., 0., 1.);
        Double_t phi_y = vvirt.Angle(vhelp);
        vvirt.RotateY(phi_y);
        vpi.RotateY(phi_y);
        phi_pq = vpi.Phi() * 180 / (TMath::Pi());
    } else {                            // Fix this in case kind != 1
        TVector3 vpi(Px(k,1), Py(k,1), Pz(k,1));
        TVector3 vvirt(-Px(0,1), -Py(0,1), kEbeam - Pz(0,1));
        Double_t phi_z = TMath::Pi() - vvirt.Phi();
        vvirt.RotateZ(phi_z);
        vpi.RotateZ(phi_z);
        TVector3 vhelp(0., 0., 1.);
        Double_t phi_y = vvirt.Angle(vhelp);
        vvirt.RotateY(phi_y);
        vpi.RotateY(phi_y);
        phi_pq = vpi.Phi() * 180 / (TMath::Pi());
    }

    return phi_pq;
}



Double_t TIdentificator::CosThetaPQ(Int_t k, Bool_t kind)
{
    if (k == 0){
        cout <<"message HERE"<<endl;
    }
    if (kind == 0)
        return (Pz(k) * (kEbeam - Pz(0)) - Px(k) * Px(0) - Py(k) * Py(0)) /
                            (sqrt(Nu() * Nu() + Q2()) * Moment(k));
    else                                // Fix this in case kind != 1
        return (Pz(k,1) * (kEbeam - Pz(0,1)) - Px(k,1) * Px(0,1) -
                Py(k,1) * Py(0,1)) /
                            (sqrt(Nu(1) * Nu(1) + Q2(1)) * Moment(k,1));
}



Double_t TIdentificator::PTrans2PQ(Int_t k, Bool_t kind)//proyeccion -> //cambiar TVector3
{   
    if (k == 0){
        cout <<"message HERE"<<endl;
    }
    if (kind == 0)
        return pow(Moment(k), 2) * (1 - pow(CosThetaPQ(k), 2));
    else                                // Fix this in case kind != 1
        return pow(Moment(k,1), 2) * (1 - pow(CosThetaPQ(k,1),2));
}
    


Double_t TIdentificator::PLong2PQ(Int_t k, Bool_t kind)
{
    if (kind == 0)
        return pow(Moment(k), 2) * pow(CosThetaPQ(k), 2);
    else                                // Fix this in case kind != 1
        return pow(Moment(k,1), 2) * pow(CosThetaPQ(k,1), 2);
}



Double_t TIdentificator::Sector(Int_t k, Bool_t kind) // Check if it is correct !!! Add k==1
{
    const Double_t r2d = 57.2957795;

    Double_t pex = Px(k);
    Double_t pey = Py(k);
    Double_t phi = TMath::ATan2(pey,pex) * r2d;

    if (phi >= -30) phi = phi+30;
    else phi = phi + 390;
    return int(phi / 60.);
}



Double_t TIdentificator::Q2(Bool_t kind) //agregar excepcion si se le da kind=0:
{
    if (kind == 0) {
        return 4. * kEbeam * Moment(0) * pow(sin(ThetaPQ(0)/2),2);
                            
    } else {                            // Fix this in case kind != 1
        return 4. * kEbeam * Moment(0,1) * pow(sin(ThetaPQ(0,1)/2),2);
    }
}



Double_t TIdentificator::W(Bool_t kind)
{
    if (kind == 0) {
        return sqrt(0.938 * 0.938 + 2 * 0.938 * Nu() - Q2());
    } else {                            // Fix this in case kind != 1
        return sqrt(0.938 * 0.938 + 2 * 0.938 * Nu(1) - Q2(1));
    }
}



Double_t TIdentificator::Nu(Bool_t kind)
{
    if (kind == 0) {
        return kEbeam - Moment(0);
    } else {                            // Fix this in case kind != 1
        return kEbeam - Moment(0,1);
    }
}




Double_t TIdentificator::ZhPi(Int_t k, Bool_t kind, Double_t Mass)
{
    if (kind == 0)
        return sqrt(pow(Mass,2) + pow(Moment(k),2)) / Nu(fCT);
    else                                // Fix this in case kind != 1
        return sqrt(pow(Mass,2) + pow(Moment(k,1),2)) / Nu(1);
}



Double_t TIdentificator::TimeCorr4(Double_t mass, Int_t k)
{
    if (k == 0){ 
        cout <<"message HERE"<<endl;
        
    }
    return (PathSC(0)/30.) - TimeSC(0) + TimeSC(k) - 0.08 -
                (PathSC(k) / 30.) * sqrt(pow(mass/Moment(k),2) + 1);
}



Double_t TIdentificator::FidTheta(Int_t k, Bool_t kind)
{
    Double_t fid_theta_val;

    if (kind == 0)
    {
        TVector3 v3p(Px(k), Py(k), Pz(k));
        fid_theta_val = v3p.Theta() * 180 / TMath::Pi();
    } else {
        TVector3 v3p(Px(k,1), Py(k,1), Pz(k,1));
        fid_theta_val = v3p.Theta() * 180 / TMath::Pi();
    }

    return fid_theta_val;
}



Double_t TIdentificator::FidThetaMin()
{
    Int_t sector = FidSector(0);

    Double_t theta_min_val = kThetaPar0[sector] +
            kThetaPar1[sector] / pow(Moment(0),2) +
            kThetaPar2[sector] * Moment(0) +
            kThetaPar3[sector] / Moment(0) +
            kThetaPar4[sector] *exp(kThetaPar5[sector] * Moment(0));

    return theta_min_val;
}



Double_t TIdentificator::FidFunc(Int_t side, Int_t param)
{
    Int_t sector = FidSector(0);
    Double_t fid_func_val = 0.0; // dummy value to avoid that uninitialized warning

    if (side == 0 && param==0)
        fid_func_val = kFidPar0Low0[sector] +
                    kFidPar0Low1[sector] * exp(kFidPar0Low2[sector] *
                            (Moment(0) - kFidPar0Low3[sector]));
    else if (side == 1 && param==0)
        fid_func_val = kFidPar0High0[sector] +
                    kFidPar0High1[sector] * exp(kFidPar0High2[sector] *
                            (Moment(0) - kFidPar0High3[sector]));
    else if (side == 0 && param==1)
        fid_func_val=kFidPar1Low0[sector] +
                    kFidPar1Low1[sector] * Moment(0) *
                    exp(kFidPar1Low2[sector] * pow((Moment(0) -
                                kFidPar1Low3[sector]),2));
    else if (side == 1 && param==1)
        fid_func_val = kFidPar1High0[sector] +
                    kFidPar1High1[sector] * Moment(0) *
                    exp(kFidPar1High2[sector] * pow((Moment(0) -
                                kFidPar1High3[sector]),2));

    return fid_func_val;
}												



Double_t TIdentificator::FidPhi(Int_t k, Bool_t kind)
{
    Double_t fid_phi_val;

    if (kind == 0) {
        TVector3 v3p(Px(k), Py(k), Pz(k));
        fid_phi_val = v3p.Phi() * 180 / TMath::Pi();
    } else {
        TVector3 v3p(Px(k,1), Py(k,1), Pz(k,1));
        fid_phi_val = v3p.Phi() * 180 / TMath::Pi();
    }

    if (fid_phi_val < -30)
        fid_phi_val += 360;
    else if (fid_phi_val > 330)
        fid_phi_val -= 360;

    return fid_phi_val;
}



Double_t TIdentificator::FidPhiMin()
{
    Int_t sector = FidSector(0);
    Double_t fid_phi_min_val;

    if (FidTheta(0) <= FidThetaMin()) {
        fid_phi_min_val = 60. * sector;
        return fid_phi_min_val;
    } else {
        fid_phi_min_val = 60. * sector - FidFunc(0,0) *
                (1 - 1 / (1 + (FidTheta(0) - FidThetaMin()) / FidFunc(0,1)));
        return fid_phi_min_val;
    }
}



Double_t TIdentificator::FidPhiMax()
{
    Int_t sector = FidSector(0);
    Double_t fid_phi_max_val;

    if (FidTheta(0) <= FidThetaMin()){
        fid_phi_max_val = 60. * sector;
        return fid_phi_max_val;
    } else {
        fid_phi_max_val = 60. * sector + FidFunc(1,0) *
                (1 - 1 / (1 + (FidTheta(0) - FidThetaMin()) / FidFunc(1,1)));
        return fid_phi_max_val;
    }
}



Bool_t TIdentificator::FidCheckCut()
{
    if (FidTheta(0) > FidThetaMin() &&
                FidPhi(0) > FidPhiMin() &&
                FidPhi(0) < FidPhiMax())
        return 1;                               // Fiducial Cut passed
    else
        return 0;                               // Fiducial Cut not passed
}



Int_t TIdentificator::FidSector(Int_t k, Bool_t kind)
{
    Int_t sector;

    if (kind == 0) {
        if (FidPhi(k) != 330) {
            sector = int((FidPhi(k) + 90) / 60) - 1;
            return sector;
        } else {
            return 5;
        }
    }
    else {
        if (FidPhi(k,1) != 330) {
            sector = int((FidPhi(k,1) + 90) / 60) - 1;
            return sector;
        } else {
            return 5;
        }
    }
}
//--------------- Added functions from old analyser----------------


Int_t TIdentificator::ElecVertTarg(Int_t k, Bool_t kind){
  Int_t p_vertex_cut_elec=0;
  Double_t ele_liq_lim[6][2];
  Double_t ele_sol_low[6];
  ele_liq_lim[0][0]=-32.5;
  ele_liq_lim[0][1]=-28;
  ele_liq_lim[1][0]=-32.5;
  ele_liq_lim[1][1]=-27.5;
  ele_liq_lim[2][0]=-32;
  ele_liq_lim[2][1]=-27.25;
  ele_liq_lim[3][0]=-32;
  ele_liq_lim[3][1]=-27.75;
  ele_liq_lim[4][0]=-32.5;
  ele_liq_lim[4][1]=-28.35;
  ele_liq_lim[5][0]=-33.5;
  ele_liq_lim[5][1]=-28.75;
  Double_t ele_sol_high=-20;
  ele_sol_low[0]=-26.5;
  ele_sol_low[1]=-26.;
  ele_sol_low[2]=-25.65;
  ele_sol_low[3]=-25.85;
  ele_sol_low[4]=-26.65;
  ele_sol_low[5]=-27.15;
  Int_t n_sector=Sector(0);

  if(Z(0)>=ele_liq_lim[n_sector][0] && Z(0)<=ele_liq_lim[n_sector][1]) p_vertex_cut_elec=1;
  if(Z(0)>=ele_sol_low[n_sector] && Z(0)<=ele_sol_high) p_vertex_cut_elec=2;

  return p_vertex_cut_elec;
}
Bool_t TIdentificator::PionVertTarg(Int_t k, Bool_t kind) {
  Bool_t vertex_cut_pion=0;
  Double_t pion_liq_low,pion_liq_high;
  Int_t n_ele_sector=Sector(0);
  Int_t n_pion_sector=Sector(k);
  if(n_pion_sector==5 || (n_ele_sector==3 && n_pion_sector==4) || (n_pion_sector==0 && n_ele_sector!=1 && n_ele_sector!=4)) pion_liq_low=-36.;
  else pion_liq_low=-35.;
  if(n_ele_sector==3 && n_pion_sector==2) pion_liq_high=-24.;
  else if((n_ele_sector==5 && n_pion_sector!=2 && n_pion_sector!=3)|| (n_pion_sector==5 && n_ele_sector!=2)|| (n_ele_sector==0 && n_pion_sector==0)|| (n_ele_sector==1 && n_pion_sector==1)|| (n_pion_sector==4 && (n_ele_sector==3 || n_ele_sector==4))) pion_liq_high=-26.;
  else pion_liq_high=-25.;

  if(ElecVertTarg(k)==1 && Z(k)>=pion_liq_low && Z(k)<=pion_liq_high) vertex_cut_pion=1;
  if(ElecVertTarg(k)==2 && Z(k)>=-30 && Z(k)<=-18) vertex_cut_pion=1;

  return vertex_cut_pion;
}


#include <../src/Categorize.C>
