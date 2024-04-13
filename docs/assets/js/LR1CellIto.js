class LR1CellIto {
    nai = 18.0
    nao = 140.0
    ki = 145.0
    ko = 5.4
    cao = 1.8
    gna = 23.0
    ibarca = 0.09
    gkr = 0.282
    prnak = 0.01833
    gk1 = 0.6047
    gpk = 0.0183
    ibbar = 0.03921
    gtof = 0.2
    gtos = 0.04
    R = 8314.0
    frdy = 96485.0
    temp = 310.0
    zna = 1.0
    zk = 1.0
    zca = 2.0
    frt = 0.037435883507803
    rtf = 26.712338705498265
    dv_max = 0.1
    adaptive = 10
    stimulus = -80.0
    stimduration = 0.5
    threshold = -75.0

    dm = 0.0
    dh = 0.0
    dj = 0.0
    dd = 0.0
    df = 0.0 
    dxtos = 0.0
    dytos = 0.0
    dxtof = 0.0
    dytof = 0.0
    dxr = 0.0
    dcai = 0.0
    ina = 0.0
    ical = 0.0
    ito = 0.0
    ikr = 0.0
    ik1 = 0.0
    ipk = 0.0
    ib = 0.0
    dv = 0.0

    v = -88.654973;
    cai = 0.0002;
    m = 0.001;
    h = 0.94;
    j = 0.94;
    xr = 0.01;
    d = 0.01;
    f = 0.95;
    xtos = 0.0;
    ytos = 1.0;
    xtof = 0.0;
    ytof = 1.0;
    diffcurrent = 0.0;
    itofac = 1.0;
    itoslowfac = 0.0;
    yshift = 0.0;
    zshift = 0.0;
    dshift = 0.0;
    fshift = 0.0;
    inafac = 1.0;
    icalfac = 1.0;
    ikfac = 1.0;
    ikifac = 1.0;
    tauXfac = 1.0;
    tauyslowfac = 1.0;
    typedestal = 0.0;

    dt = 0.01;
    
    constructor() {
        
    }
    iterate(st, DV_MAX) {
        this.comp_ina();
        this.comp_ical();
        this.comp_ito();
        this.comp_ikr();
        this.comp_ik1();
        this.comp_ipk();
        this.comp_ib();
        this.dv = (this.diffcurrent - (this.ikr + this.ik1 + this.ito + this.ina + this.ical + this.ib + this.ipk + st))*this.dt;

        if (DV_MAX > 0 & this.dv*this.dv > DV_MAX*DV_MAX) {
            return false;
        }
        else {
            this.comp_calcdyn();
            
            this.v = this.v + this.dv;
            this.m = this.m + this.dm*this.dt;
            this.h = this.h + this.dh*this.dt;
            this.j = this.j + this.dj*this.dt;
            this.d = this.d + this.dd*this.dt;
            this.f = this.f + this.df*this.dt;
            this.xtos = this.xtos + this.dxtos*this.dt;
            this.ytos = this.ytos + this.dytos*this.dt;
            this.xtof = this.xtof + this.dxtof*this.dt;
            this.ytof = this.ytof + this.dytof*this.dt;
            this.xr = this.xr + this.dxr*this.dt;
            this.cai = this.cai + this.dcai*this.dt;
        }
        return true;
    };
    stepdt(st) {
        var success = this.iterate(st, this.dv_max);
        if (!success) {
            this.dt = this.dt/this.adaptive;
            for (let i = 0; i < this.adaptive; i++) {
                this.iterate(st, -1);
            }
            this.dt = this.dt*this.adaptive;
        }
    };
    comp_ina() {
        const ena = (this.rtf/this.zna)*Math.log(this.nao/this.nai);
        const am = 0.32*(this.v+47.13)/(1.0-Math.exp(-0.1*(this.v+47.13)));
        const bm = 0.08*Math.exp(-this.v/11.0);
        let ah, bh, aj, bj = 0.0;
    
        if (this.v < -40.0) {
            ah = 0.135*Math.exp((80+this.v)/-6.8);
            bh = 3.56*Math.exp(0.079*this.v)+310000*Math.exp(0.35*this.v);
            aj = (-127140*Math.exp(0.2444*this.v)-0.00003474*Math.exp(-0.04391*this.v))*((this.v+37.78)/(1.0+Math.exp(0.311*(this.v+79.23))));
            bj = (0.1212*Math.exp(-0.01052*this.v))/(1.0+Math.exp(-0.1378*(this.v+40.14)));
        }
        else {
            ah = 0.0;
            bh = 1.0/(0.13*(1.0+Math.exp((this.v+10.66)/-11.1)));
            aj = 0.0;
            bj = (0.3*Math.exp(-0.0000002535*this.v))/(1.0+Math.exp(-0.1*(this.v+32.0)));
        }   

        const mtau = 1.0/(am+bm);
        const htau = 1.0/(ah+bh);
        const jtau = 1.0/(aj+bj);
    
        const mss = am*mtau;
        const hss = ah*htau;
        const jss = aj*jtau;

        this.ina = this.inafac*this.gna*this.m*this.m*this.m*this.h*this.j*(this.v-ena);

        this.dm = (mss-(mss-this.m)*Math.exp(-this.dt/mtau) - this.m)/this.dt;
        this.dh = (hss-(hss-this.h)*Math.exp(-this.dt/htau) - this.h)/this.dt;
        this.dj = (jss-(jss-this.j)*Math.exp(-this.dt/jtau) - this.j)/this.dt;
    };
    comp_ical() {
        let ad, bd, af, bf, taud, tauf, dss, fss, ical = 0;
    
        ad = 0.095*Math.exp(-0.01*(this.v + this.dshift-5.0))/(1.0+Math.exp(-0.072*(this.v + this.dshift-5.0)));
        bd = 0.07*Math.exp(-0.017*(this.v + this.dshift+44.0))/(1.0+Math.exp(0.05*(this.v + this.dshift+44.0)));
        af = 0.012*Math.exp(-0.008*(this.v - this.fshift+28.0))/(1.0+Math.exp(0.15*(this.v-this.fshift+28.0)));
        bf = 0.0065*Math.exp(-0.02*(this.v - this.fshift+30.0))/(1.0+Math.exp(-0.2*(this.v - this.fshift+30.0)));
        taud = 1.0/(ad+bd);
        tauf = 1.0/(af+bf);
        dss = ad*taud;
        fss = af*tauf;
        //taud = TAUDFAC*taud;
        //tauf = TAUFFAC*tauf;
    
        ical = this.ibarca*this.d*this.f*(this.v - (7.7 - 13.0287*Math.log(this.cai)));
    
        this.dd = (dss - (dss - this.d)*Math.exp(-this.dt/taud) - this.d)/this.dt;
        this.df = (fss - (fss - this.f)*Math.exp(-this.dt/tauf) - this.f)/this.dt;
        this.ical = this.icalfac*ical;
    }
    comp_ito() {
        let ek, rt1, rt2, rt3, xtos_inf, ytos_inf, rs_inf, txs, tys, xitos = 0;
        let xtof_inf, ytof_inf, rt4, rt5, txf, tyf, xitof = 0;
        ek = (this.rtf/this.zk)*Math.log(this.ko/this.ki);
        rt1 = -(this.v + 3.0 + this.zshift)/15.0;
        rt2 = (this.v + 33.5 - this.yshift)/10.0;
        rt3 = (this.v + 60.0)/10.0;
        xtos_inf = 1.0/(1.0 + Math.exp(rt1));
        ytos_inf = 1.0/(1.0 + Math.exp(rt2));
        rs_inf = 1.0/(1.0 + Math.exp(rt2));
        txs = 9.0/(1.0 + Math.exp(-rt1)) + 0.5;
        tys = 3000.0/(1.0 + Math.exp(rt3)) + 30.0;
        tys = this.tauyslowfac*tys;
        if (tys < this.typedestal) {
            tys = this.typedestal;  
        }
        //tys = tauyslowfac[id]*tys;
        xitos = this.gtos*this.xtos*(this.ytos + 0.5*rs_inf)*(this.v - ek);
        this.dxtos = (xtos_inf - (xtos_inf - this.xtos)*Math.exp(-this.dt/txs) - this.xtos)/this.dt;
        this.dytos = (ytos_inf - (ytos_inf - this.ytos)*Math.exp(-this.dt/tys) - this.ytos)/this.dt;
    
        xtof_inf = xtos_inf;
        ytof_inf = ytos_inf;
        rt4 = -((this.v + this.zshift)/30.0)*((this.v + this.zshift)/30.0);
        rt5 = (this.v + 33.5 - this.yshift)/10.0;
        txf = 3.5*Math.exp(rt4) + 1.5;
        tyf = 20.0/(1.0 + Math.exp(rt5)) + 20.0;
        xitof = this.gtof*this.xtof*this.ytof*(this.v - ek);
        this.dxtof = (xtof_inf - (xtof_inf - this.xtof)*Math.exp(-this.dt/txf) - this.xtof)/this.dt;
        this.dytof = (ytof_inf - (ytof_inf - this.ytof)*Math.exp(-this.dt/tyf) - this.ytof)/this.dt;
        this.ito = this.itoslowfac*xitos + this.itofac*xitof;
    }
    comp_ikr() {
        let ekr, r, ax, bx, tauxr, xrss, ikr = 0.0;
    
        ekr = (this.rtf/this.zk)*Math.log((this.ko + this.prnak*this.nao)/(this.ki + this.prnak*this.nai));
    
        r = (this.v > -100.0) ? 2.837*(Math.exp(0.04*(this.v + 77.0)) - 1.0)/((this.v + 77.0)*Math.exp(0.04*(this.v + 35.0))) : 1.0;
    
        ax = 0.0005*Math.exp(0.083*(this.v + 50))/(1.0 + Math.exp(0.057*(this.v + 50)));
        bx = 0.0013*Math.exp(-0.06*(this.v + 20))/(1.0 + Math.exp(-0.04*(this.v + 20)));
        tauxr = 1.0/(ax+bx);
        xrss = ax*tauxr; //This is the proper order, do this before tauxr is multiplied by the tauX factor
        tauxr = this.tauXfac*tauxr;
    
        this.ikr = this.ikfac*this.gkr*Math.sqrt(this.ko/5.4)*this.xr*r*(this.v - ekr);
    
        this.dxr = (xrss - (xrss - this.xr)*Math.exp(-this.dt/tauxr) - this.xr)/this.dt;
    }
    comp_ik1() {
        let ek, ak1, bk1, xk1ss, ik1 = 0;
    
        ek = (this.rtf/this.zk)*Math.log(this.ko/this.ki);
    
        ak1 = 1.02/(1.0+Math.exp(0.2385*(this.v - ek - 59.215)));
        bk1 = (0.49124*Math.exp(0.08032*(this.v - ek+5.476))+Math.exp(0.06175*(this.v - ek - 594.31)))/(1.0+Math.exp(-0.5143*(this.v - ek+4.753)));
        xk1ss = ak1/(ak1 + bk1);
        this.ik1 = this.ikifac*this.gk1*Math.sqrt(this.ko/5.4)*xk1ss*(this.v - ek);
    }
    comp_ipk() {
        let ek, ipk, kp = 0;
    
        ek = (this.rtf/this.zk)*Math.log(this.ko/this.ki);
        kp = 1.0/(1.0+Math.exp((7.488-this.v)/5.98));
        this.ipk = this.gpk*kp*(this.v - ek);
    }
    comp_ib() {
        this.ib = this.ibbar*(this.v + 59.87);
    }
    comp_calcdyn() {
        this.dcai = -0.0001*this.ical + 0.07*(0.0001 - this.cai);
    }
    setcell() {

    }
    getcell() {

    }
    saveconditions() {

    }
}