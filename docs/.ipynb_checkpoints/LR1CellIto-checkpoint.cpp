//
//  LR1CellIto.cpp
//
//  Implementation of the LR1 Model, with UCLA Ito
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#include <math.h>
#include <fstream>
#include "LR1CellIto.h"

#ifndef LR1CellIto_cpp
#define LR1CellIto_cpp

#define nai 18.0
#define nao 140.0
#define ki 145.0
#define ko 5.4
#define cao 1.8

#define gna 23.0

#define ibarca 0.09

#ifndef IBARCAFAC
#define IBARCAFAC 1.0
#endif

#ifndef dshift
#define dshift 0.0
#endif

#ifndef fshift
#define fshift 0.0
#endif

#ifndef ZSHIFT
#define ZSHIFT 0.0
#endif

#ifndef YSHIFT
#define YSHIFT 0.0
#endif

#ifndef IKFAC
#define IKFAC 1.0
#endif

#ifndef IKIFAC
#define IKIFAC 1.0
#endif

#ifndef TAUXFAC
#define TAUXFAC 1.0
#endif

#define gkr 0.282
#define prnak 0.01833
#define gk1 0.6047
#define gpk 0.0183
#define ibbar 0.03921

#define gtof 0.2

#ifdef slowito
#define gtos 0.073
#else
#define gtos 0.04
#endif

#define R 8314.0
#define frdy 96485.0
#define temp 310.0
#define zna 1.0
#define zk 1.0
#define zca 2.0
#define frt 0.037435883507803
#define rtf 26.712338705498265

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

#ifndef stimulus
#define stimulus -80.0
#endif

#ifndef stimduration
#define stimduration 0.5
#endif

template <int ncells>
 LR1CellIto<ncells>::LR1CellIto()
{
    for (int i = 0; i < ncells; i++) {
        v[i] = -88.654973;
        cai[i] = 0.0002;
        m[i] = 0.001;
        h[i] = 0.94;
        j[i] = 0.94;
        xr[i] = 0.01;
        d[i] = 0.01;
        f[i] = 0.95;
        xtos[i] = 0.0;
        ytos[i] = 1.0;
        xtof[i] = 0.0;
        ytof[i] = 1.0;
        diffcurrent[i] = 0.0;
        itofac[i] = 1.0;
        itoslowfac[i] = 0.0;
        yshift[i] = YSHIFT;
        zshift[i] = ZSHIFT;
        inafac[i] = 1.0;
        icalfac[i] = IBARCAFAC;
        ikfac[i] = IKFAC;
        ikifac[i] = IKIFAC;
        tauXfac[i] = TAUXFAC;
        tauyslowfac[i] = 1.0;
        typedestal[i] = 0.0;
    }
}

template <int ncells>
  bool LR1CellIto<ncells>::iterate(const int id, double dt, double st, double dv_max) {
    double dm, dh, dj, dd, df, dxtos, dytos, dxtof, dytof, dxr, dcai;
    double ina, ical, ito, ikr, ik1, ipk, ib, dv;
    
    ina = comp_ina(id, dt, dm, dh, dj);
    ical = comp_ical(id, dt, dd, df);
    ito = comp_ito(id, dt, dxtos, dytos, dxtof, dytof);
    ikr = comp_ikr(id, dt, dxr);
    ik1 = comp_ik1(id);
    ipk = comp_ipk(id);
    ib = comp_ib(id);
    dv = (diffcurrent[id] - (ikr + ik1 + ito + ina + ical + ib + ipk + st))*dt;
    
    //Comment this out to remove adaptive time-stepping
    if (dv_max > 0 && dv*dv > dv_max*dv_max) {return false;}
      
    comp_calcdyn(id, ical, dcai);
    
    v[id] += dv;
    m[id] += dm*dt;
    h[id] += dh*dt;
    j[id] += dj*dt;
    d[id] += dd*dt;
    f[id] += df*dt;
    xtos[id] += dxtos*dt;
    ytos[id] += dytos*dt;
    xtof[id] += dxtof*dt;
    ytof[id] += dytof*dt;
    xr[id] += dxr*dt;
    cai[id] += dcai*dt;
    
    return true;
}

template <int ncells>
void LR1CellIto<ncells>::stepdt (const int id, double dt, double st) {
    if (id > -1 && id < ncells) {
        bool success = iterate(id, dt, st, DV_MAX);
        if (!success) {
            for (int i = 0; i < ADAPTIVE; i++) {
                iterate(id, dt/ADAPTIVE, st, -1);
            }
        }
    }
}

template <int ncells>
  double LR1CellIto<ncells>::comp_ina (int id, double dt, double& dm, double& dh, double& dj) //Fast Sodium Current
{
    double ena, mtau, htau, jtau, mss, hss, jss, ina;
    double am, bm, ah, bh, aj, bj;
    
    ena = (rtf/zna)*log(nao/nai);
    
    am = 0.32*(v[id]+47.13)/(1-exp(-0.1*(v[id]+47.13)));
    bm = 0.08*exp(-v[id]/11.0);
    
    if (v[id] < -40.0) {
        ah = 0.135*exp((80+v[id])/-6.8);
        bh = 3.56*exp(0.079*v[id])+310000*exp(0.35*v[id]);
        aj = (-127140*exp(0.2444*v[id])-0.00003474*exp(-0.04391*v[id]))*((v[id]+37.78)/(1.0+exp(0.311*(v[id]+79.23))));
        bj = (0.1212*exp(-0.01052*v[id]))/(1.0+exp(-0.1378*(v[id]+40.14)));
    }
    else {
        ah = 0.0;
        bh = 1.0/(0.13*(1.0+exp((v[id]+10.66)/-11.1)));
        aj = 0.0;
        bj = (0.3*exp(-0.0000002535*v[id]))/(1.0+exp(-0.1*(v[id]+32.0)));
    }
    
    mtau = 1.0/(am+bm);
    htau = 1.0/(ah+bh);
    jtau = 1.0/(aj+bj);
    
    mss = am*mtau;
    hss = ah*htau;
    jss = aj*jtau;
    
    ina = inafac[id]*gna*m[id]*m[id]*m[id]*h[id]*j[id]*(v[id]-ena);
    
    dm = (mss-(mss-m[id])*exp(-dt/mtau) - m[id])/dt;
    dh = (hss-(hss-h[id])*exp(-dt/htau) - h[id])/dt;
    dj = (jss-(jss-j[id])*exp(-dt/jtau) - j[id])/dt;
    
    return ina;
}

#ifndef TAUDFAC
#define TAUDFAC 1.0
#endif

#ifndef TAUFFAC
#define TAUFFAC 1.0
#endif

template <int ncells>
  double LR1CellIto<ncells>::comp_ical (int id, double dt, double& dd, double& df) // L-type Calcium Current
{
    double ad, bd, af, bf, taud, tauf, dss, fss, ical;
    
    ad = 0.095*exp(-0.01*(v[id] + dshift-5.0))/(1.0+exp(-0.072*(v[id] + dshift-5.0)));
    bd = 0.07*exp(-0.017*(v[id] + dshift+44.0))/(1.0+exp(0.05*(v[id] + dshift+44.0)));
    af = 0.012*exp(-0.008*(v[id]-fshift+28.0))/(1.0+exp(0.15*(v[id]-fshift+28.0)));
    bf = 0.0065*exp(-0.02*(v[id]-fshift+30.0))/(1.0+exp(-0.2*(v[id]-fshift+30.0)));
    taud = 1.0/(ad+bd);
    tauf = 1.0/(af+bf);
    dss = ad*taud;
    fss = af*tauf;
    taud = TAUDFAC*taud;
    tauf = TAUFFAC*tauf;
    
    ical = ibarca*d[id]*f[id]*(v[id] - (7.7 - 13.0287*log(cai[id])));
    
    dd = (dss - (dss - d[id])*exp(-dt/taud) - d[id])/dt;
    df = (fss - (fss - f[id])*exp(-dt/tauf) - f[id])/dt;
    
    return icalfac[id]*ical;
}

template <int ncells>
double LR1CellIto<ncells>::comp_ito (int id, double dt, double& dxtos, double& dytos, double& dxtof, double& dytof)
{
#ifndef slowito
    double ek, rt1, rt2, rt3, xtos_inf, ytos_inf, rs_inf, txs, tys, xitos;
    double xtof_inf, ytof_inf, rt4, rt5, txf, tyf, xitof;
    
    ek = (rtf/zk)*log(ko/ki);
    rt1 = -(v[id] + 3.0 + zshift[id])/15.0;
    rt2 = (v[id] + 33.5 - yshift[id])/10.0;
    rt3 = (v[id] + 60.0)/10.0;
    xtos_inf = 1.0/(1.0 + exp(rt1));
    ytos_inf = 1.0/(1.0 + exp(rt2));
    rs_inf = 1.0/(1.0 + exp(rt2));
    txs = 9.0/(1.0 + exp(-rt1)) + 0.5;
    tys = 3000.0/(1.0 + exp(rt3)) + 30.0;
    tys = tauyslowfac[id]*tys;
    if (tys < typedestal[id]) {
        tys = typedestal[id];  
    }
    //tys = tauyslowfac[id]*tys;
    xitos = gtos*xtos[id]*(ytos[id] + 0.5*rs_inf)*(v[id] - ek);
    dxtos = (xtos_inf - (xtos_inf - xtos[id])*exp(-dt/txs) - xtos[id])/dt;
    dytos = (ytos_inf - (ytos_inf - ytos[id])*exp(-dt/tys) - ytos[id])/dt;
    
    xtof_inf = xtos_inf;
    ytof_inf = ytos_inf;
    rt4 = -((v[id] + zshift[id])/30.0)*((v[id] + zshift[id])/30.0);
    rt5 = (v[id] + 33.5 - yshift[id])/10.0;
    txf = 3.5*exp(rt4) + 1.5;
    tyf = 20.0/(1.0 + exp(rt5)) + 20.0;
    xitof = gtof*xtof[id]*ytof[id]*(v[id] - ek);
    dxtof = (xtof_inf - (xtof_inf - xtof[id])*exp(-dt/txf) - xtof[id])/dt;
    dytof = (ytof_inf - (ytof_inf - ytof[id])*exp(-dt/tyf) - ytof[id])/dt;
    
    return itoslowfac[id]*xitos + itofac[id]*xitof;

#else
    double ek, zssdv, tauzdv, yssdv, tauydv, ito;
    ek = (rtf/zk)*log(ko/ki);
    
    zssdv = 1.0/(1.0 + exp((20.0 - (v[id] - zshift[id]))/6.0));
    tauzdv = 9.5*exp(-((v[id] - zshift[id]) + 40.0)*((v[id] - zshift[id]) + 40.0)/1800.0) + 0.8;
                 
    yssdv = 1.0/(1.0 + exp(((v[id] - yshift[id]) + 28.0)/5.0));
    tauydv = 1000.0*exp(-((v[id] - yshift[id]) + 67.0)*((v[id] - yshift[id]) + 67.0)/1000.0) + 8.0;
    tauydv = tauyslowfac[id]*tauydv;
    
    //ito = itoslowfac[id]*gitodv*xtos[id]*ytos[id]*(v[id]-ek);
    ito = itoslowfac[id]*gtof*xtos[id]*ytos[id]*(v[id]-ek);
    
    dxtos = (zssdv-(zssdv-xtos[id])*exp(-dt/tauzdv) - xtos[id])/dt; 
    dytos = (yssdv-(yssdv-ytos[id])*exp(-dt/tauydv) - ytos[id])/dt;
    
    dxtof = 0;
    dytof = 0;
    
    return ito;  
    
#endif
}

template <int ncells>
double LR1CellIto<ncells>::comp_ikr (int id, double dt, double& dxr)
{
    double ekr, r, ax, bx, tauxr, xrss,  ikr;
    
    ekr = (rtf/zk)*log((ko + prnak*nao)/(ki + prnak*nai));
    
    r = (v[id] > -100.0) ? 2.837*(exp(0.04*(v[id] + 77.0)) - 1.0)/((v[id] + 77.0)*exp(0.04*(v[id] + 35.0))) : 1.0;
    
    ax = 0.0005*exp(0.083*(v[id] + 50))/(1.0 + exp(0.057*(v[id] + 50)));
    bx = 0.0013*exp(-0.06*(v[id]+20))/(1.0 + exp(-0.04*(v[id]+20)));
    tauxr = 1.0/(ax+bx);
    xrss = ax*tauxr; //This is the proper order, do this before tauxr is multiplied by the tauX factor
    tauxr = tauXfac[id]*tauxr;
    
    ikr = ikfac[id]*gkr*sqrt(ko/5.4)*xr[id]*r*(v[id] - ekr);
    
    dxr = (xrss - (xrss - xr[id])*exp(-dt/tauxr) - xr[id])/dt;
    
    return ikr;
}

template <int ncells>
double LR1CellIto<ncells>::comp_ik1 (int id)
{
    double ek, ak1, bk1, xk1ss, ik1;
    
    ek = (rtf/zk)*log(ko/ki);
    
    ak1 = 1.02/(1.0+exp(0.2385*(v[id]-ek-59.215)));
    bk1 = (0.49124*exp(0.08032*(v[id]-ek+5.476))+exp(0.06175*(v[id]-ek-594.31)))/(1.0+exp(-0.5143*(v[id]-ek+4.753)));
    xk1ss = ak1/(ak1 + bk1);
    ik1 = ikifac[id]*gk1*sqrt(ko/5.4)*xk1ss*(v[id] - ek);
    
    return ik1;
}

template <int ncells>
double LR1CellIto<ncells>::comp_ipk (int id)
{
    double ek, ipk, kp;
    
    ek = (rtf/zk)*log(ko/ki);
    kp = 1.0/(1.0+exp((7.488-v[id])/5.98));
    ipk = gpk*kp*(v[id]-ek);
    return ipk;
}

template <int ncells>
double LR1CellIto<ncells>::comp_ib (int id)
{
    return ibbar*(v[id] + 59.87);
}

template <int ncells>
void LR1CellIto<ncells>::comp_calcdyn (int id, double ical, double& dcai) // MAKE SURE THIS IS COMPUTED AFTER ALL OTHER CURRENTS
{
    dcai = -0.0001*ical + 0.07*(0.0001 - cai[id]);
}

template <int ncells>
void LR1CellIto<ncells>::setcell (int id, LR1CellIto<1>* newcell)
{
    v[id] = newcell->v[0];
    cai[id] = newcell->cai[0];
    m[id] = newcell->m[0];
    h[id] = newcell->h[0];
    j[id] = newcell->j[0];
    xr[id] = newcell->xr[0];
    d[id] = newcell->d[0];
    f[id] = newcell->f[0];
    xtos[id] = newcell->xtos[0];
    ytos[id] = newcell->ytos[0];
    xtof[id] = newcell->xtof[0];
    ytof[id] = newcell->ytof[0];
    diffcurrent[id] = newcell->diffcurrent[0];
    itofac[id] = newcell->itofac[0];
    itoslowfac[id] = newcell->itoslowfac[0];
    tauXfac[id] = newcell->tauXfac[0];
    icalfac[id] = newcell->icalfac[0];
    ikfac[id] = newcell->ikfac[0];
    ikifac[id] = newcell->ikifac[0];
    yshift[id] = newcell->yshift[0];
    zshift[id] = newcell->zshift[0];
    inafac[id] = newcell->inafac[0]; 
    tauyslowfac[id] = newcell->tauyslowfac[0]; 
    typedestal[id] = newcell->typedestal[0]; 
}

template <int ncells>
void LR1CellIto<ncells>::getcell (int id, LR1CellIto<1>* newcell)
{
    newcell->v[0] = v[id];
    newcell->cai[0] = cai[id];
    newcell->m[0] = m[id];
    newcell->h[0] = h[id];
    newcell->j[0] = j[id];
    newcell->xr[0] = xr[id];
    newcell->d[0] = d[id];
    newcell->f[0] = f[id];
    newcell->xtos[0] = xtos[id];
    newcell->ytos[0] = ytos[id];
    newcell->xtof[0] = xtof[id];
    newcell->ytof[0] = ytof[id];
    newcell->diffcurrent[0] = diffcurrent[id];
    newcell->itofac[0] = itofac[id];
    newcell->itoslowfac[0] = itoslowfac[id];
    newcell->tauXfac[0] = tauXfac[id];
    newcell->icalfac[0] = icalfac[id];
    newcell->ikfac[0] = ikfac[id];
    newcell->ikifac[0] = ikifac[id];
    newcell->yshift[0] = yshift[id];
    newcell->zshift[0] = zshift[id];
    newcell->inafac[0] = inafac[id];
    newcell->tauyslowfac[0] = tauyslowfac[id];
    newcell->typedestal[0] = typedestal[id];
}

template <int ncells>
void LR1CellIto<ncells>::saveconditions(FILE* file, int id, bool header, double t) {
    if (header) {
        fprintf(file,"t\tv\tm\th\tj\td\tf\txtos\tytos\txtof\tytof\txr\tcai\n");
    }
    fprintf(file,"%g\t",t);
    fprintf(file,"%.12f\t",v[id]);
    fprintf(file,"%.12f\t",m[id]);
    fprintf(file,"%.12f\t",h[id]);
    fprintf(file,"%.12f\t",j[id]);
    fprintf(file,"%.12f\t",d[id]);
    fprintf(file,"%.12f\t",f[id]);
    fprintf(file,"%.12f\t",xtos[id]);
    fprintf(file,"%.12f\t",ytos[id]);
    fprintf(file,"%.12f\t",xtof[id]);
    fprintf(file,"%.12f\t",ytof[id]);
    fprintf(file,"%.12f\t",xr[id]);
    fprintf(file,"%.12f\n",cai[id]);
}
#undef nai
#undef nao
#undef ki
#undef ko
#undef cao
#undef gna
#undef ibarca
#undef gkr
#undef prnak
#undef gk1
#undef gpk
#undef ibbar
#undef gitodv
#undef R
#undef frdy
#undef temp
#undef zna
#undef zk
#undef zca
#undef frt
#undef rtf

#endif // LR1CellIto_cpp

