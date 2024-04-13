class LR1CellIto {
    nai = 18.0
    nao = 140.0
    ki = 145.0
    ko = 5.4
    cao = 1.8
    gna = 23.0
    ibarca = 0.09
    IBARCAFAC = 1.0
    dshift = 0.0
    fshift = 0.0
    ZSHIFT = 0.0
    YSHIFT = 0.0
    IKFAC = 1.0
    IKIFAC = 1.0
    TAUXFAC = 1.0
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

    dm = 0.0
    dh = 0.0
    dj = 0.0
    dd = 0.0
    df = 0.0 dxtos = 0.0
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
    yshift = YSHIFT;
    zshift = ZSHIFT;
    inafac = 1.0;
    icalfac = IBARCAFAC;
    ikfac = IKFAC;
    ikifac = IKIFAC;
    tauXfac = TAUXFAC;
    tauyslowfac = 1.0;
    typedestal = 0.0;

    dt = 0.01;
    
    constructor() {
        
    }
    iterate(st) {
        comp_ina();
        comp_ical();
        comp_ito();
        comp_ikr();
        comp_ik1();
        comp_ipk();
        comp_ib();
        this.dv = (this.diffcurrent - (this.ikr + this.ik1 + this.ito + this.ina + this.ical + this.ib + this.ipk + st))*this.dt;

        if (this.dv_max > 0 & this.dv*this.dv > this.dv_max*this.dv_max) {
            return false;
        }
        else {
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
    }
    stepdt(st) {
        success = iterate(st);
        if (!success) {
            this.dt = this.dt/this.adaptive;
            for (int i = 0; i < this.adaptive; i++) {
                iterate(st);
            }
            this.dt = this.dt*this.adaptive;
        }
    }
    comp_ina() {
        const ena = (this.rtf/this.zna)*Math.log(this.nao/this.nai);
        const am = 0.32*(this.v+47.13)/(1.0-Math.exp(-0.1*(this.v+47.13)));
        const bm = 0.08*Math.exp(-this.v/11.0);
    
        if (this.v < -40.0) {
            const ah = 0.135*Math.exp((80+this.v)/-6.8);
            const bh = 3.56*Math.exp(0.079*this.v[id])+310000*Math.exp(0.35*this.v);
            const aj = (-127140*Math.exp(0.2444*this.v)-0.00003474*Math.exp(-0.04391*this.v))*((this.v+37.78)/(1.0+Math.exp(0.311*(this.v+79.23))));
            const bj = (0.1212*Math.exp(-0.01052*this.v))/(1.0+Math.exp(-0.1378*(this.v+40.14)));
        }
        else {
            const ah = 0.0;
            const bh = 1.0/(0.13*(1.0+Math.exp((this.v+10.66)/-11.1)));
            const aj = 0.0;
            const bj = (0.3*Math.exp(-0.0000002535*this.v))/(1.0+Math.exp(-0.1*(this.v+32.0)));
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
    }
    comp_ical() {

    }
    comp_ito() {

    }
    comp_ikr() {

    }
    comp_ik1() {

    }
    comp_ipk() {

    }
    comp_ib() {

    }
    comp_calcdyn() {

    }
    setcell() {

    }
    getcell() {

    }
    saveconditions() {

    }
}