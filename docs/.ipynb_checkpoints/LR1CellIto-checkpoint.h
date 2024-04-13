//
//  LR1CellIto.h
//
//  Implementation of the LR1 model, with UCLA Ito
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef LR1CellIto_h
#define LR1CellIto_h

template <int ncells>
class LR1CellIto
{
public:
    double v[ncells]; //Membrane voltage
    double diffcurrent[ncells];
    double cai[ncells];
    
    /* Fast Sodium Current (time dependant) */
    double m[ncells]; // Na activation
    double h[ncells]; // Na inactivation
    double j[ncells]; // Na inactivation
    
    /* L-type Calcium Window Current ([Ca] dependent) */
    double d[ncells];
    double f[ncells];
    
    /* Ito Transient Outward Current (Dumaine et al. Circ Res 1999;85:803-809) */
    double xtos[ncells];
    double ytos[ncells];
    double xtof[ncells];
    double ytof[ncells];
    
    //double zdv[ncells]; // Ito activation
    //double ydv[ncells]; // Ito inactivation
    double itofac[ncells];
    double itoslowfac[ncells];
    double zshift[ncells];
    double yshift[ncells];
    double inafac[ncells];
    double icalfac[ncells];
    double ikfac[ncells];
    double ikifac[ncells];
    double tauXfac[ncells];
    
    double xr[ncells];
    double tauyslowfac[ncells];
    double typedestal[ncells];
    
    LR1CellIto();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    
    void stepdt(const int id, double dt, double st);
    
    double comp_ina (int id, double dt, double& dm, double& dh, double& dj); //Fast Sodium Current
    
    double comp_ical (int id, double dt, double& dd, double& df); // L-type Calcium Current
    
    double comp_ito (int id, double dt, double& dxtos, double& dytos, double& dxtof, double& dytof);
    
    double comp_ikr (int id, double dt, double& dxr);
    
    double comp_ik1 (int id);
    
    double comp_ipk (int id);
    
    double comp_ib (int id);
    
    void comp_calcdyn (int id, double ical, double& dcai); // MAKE SURE THIS IS COMPUTED AFTER ALL OTHER CURRENTS
    
    void setcell (int id, LR1CellIto<1>* newcell);
    
    void getcell (int id, LR1CellIto<1>* newcell);
    
    void saveconditions(FILE* file, int id, bool header, double t=0.0);
};

#endif // LR1CellIto_h

