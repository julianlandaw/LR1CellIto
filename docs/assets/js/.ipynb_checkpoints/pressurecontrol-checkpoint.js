var PCslider = document.getElementById("PCslider");
var PChtml = document.getElementById("PChtml");
PChtml.innerHTML = "Pressure Control (cm H<sub>2</sub>O)";
var PCnum = document.getElementById("PC");
PCnum.value = PCslider.value/10;

var PEEPslider = document.getElementById("PEEPslider");
var PEEPhtml = document.getElementById("PEEPhtml");
PEEPhtml.innerHTML = "PEEP (cm H<sub>2</sub>O)";
var PEEPnum = document.getElementById("PEEP");
PEEPnum.value = PEEPslider.value/10;

var RRslider = document.getElementById("RRslider");
var RRhtml = document.getElementById("RRhtml");
RRhtml.innerHTML = "RR (breaths/min)";
var RRnum = document.getElementById("RR");
RRnum.value = RRslider.value/10;

var IFslider = document.getElementById("IFslider");
var IFhtml = document.getElementById("IFhtml");
IFhtml.innerHTML = "Inspiratory Fraction";
var IFnum = document.getElementById("IF");
IFnum.value = IFslider.value/100;

var Rslider = document.getElementById("Rslider");
var Rhtml = document.getElementById("Rhtml");
Rhtml.innerHTML = "Resistance R (cm H<sub>2</sub>O*s/mL)";
var Rnum = document.getElementById("R");
Rnum.value = Rslider.value/10000;

var Cslider = document.getElementById("Cslider");
var Chtml = document.getElementById("Chtml");
Chtml.innerHTML = "Compliance C (mL/cm H<sub>2</sub>O)";
var Cnum = document.getElementById("C");
Cnum.value = Cslider.value/10;

var tfinalslider = document.getElementById("tfinalslider");
var tfinalhtml = document.getElementById("tfinalhtml");
tfinalhtml.innerHTML = "Time";
var tfinalnum = document.getElementById("tfinal");
tfinalnum.value = tfinalslider.value/10;

function dfsolve() {
    // Exactly computes 

    function expn(x) {
        return math.exp(x); //1 + x + x*x/2 + x*x*x/6 + x*x*x*x/24 + x*x*x*x*x/120 + x*x*x*x*x*x/720 + x*x*x*x*x*x*x/5040;
    }

    function iterate1(t, C1, C2, C3, params) { //C1 Pressure, C2 is Volume, C3 is Flow
        const pc = params.PC;
        const peep = params.PEEP;
        const RR = params.RR;
        const ifraction = params.IF;
        const R = params.R;
        const C = params.C;
        const breathduration = ifraction*60.0/RR;
        const timenextbreath = params.timenextbreath;
        const dt = params.dt;
        if (t > timenextbreath - dt/2 & t < timenextbreath + breathduration) {
            return pc;
        }
        else {
            return peep;
        }
    }

    function iterate2(t, C1, C2, C3, params) {
        const pc = params.PC;
        const peep = params.PEEP;
        const RR = params.RR;
        const ifraction = params.IF;
        const R = params.R;
        const C = params.C;
        const breathduration = ifraction*60.0/RR;
        const timenextbreath = params.timenextbreath;
        const dt = params.dt;
        return C*C1 + (C2 - C*C1)*expn(-dt/(R*C));
        //return C2 + dt*C3;    //C2*expn(-dt/(R*C)) + (C1/R)*(expn(dt/(R*C)) - 1.0)
    }

    function iterate3(t, C1, C2, C3, params) {
        const pc = params.PC;
        const peep = params.PEEP;
        const RR = params.RR;
        const ifraction = params.IF;
        const R = params.R;
        const C = params.C;
        const breathduration = ifraction*60.0/RR;
        const timenextbreath = params.timenextbreath;
        const dt = params.dt;
        return (C1/R - C2/(R*C))*expn(-dt/(R*C)); //[R] = [cm H2O*s/mL], [C] = [mL/cm H2O]
    }
    
    let params = { PC: [], PEEP: [], RR: [], IF: [], R: [], C: [], timenextbreath: [], dt: []};
    
    params.PC = PCslider.value/10.0;
    params.PEEP = PEEPslider.value/10.0;
    params.RR = RRslider.value/10.0;
    params.IF = IFslider.value/100.0;
    params.R = Rslider.value/10000.0;
    params.C = Cslider.value/10.0;
    params.timenextbreath = 0;
    params.dt = 0.1;
    params.tfinal = tfinalslider.value/10.0;

    let t0 = 0;
    let y01 = params.PEEP;
    let y02 = params.C*params.PEEP;
    let y03 = 0;

    let N = Math.ceil(params.tfinal/params.dt);
    let breathN = Math.ceil(params.IF*60.0/params.RR);

    let ts = new Array(N + 1);
    let ys1 = new Array(N + 1);
    let ys2 = new Array(N + 1);
    let ys3 = new Array(N + 1);

    for (let i = 0; i < N+1; i++) {
       ts[i] = i*params.dt;
    }

    ys1[0] = y01;
    ys2[0] = y02;
    ys3[0] = y03;

    let counter = 0;
    while (counter < N) {
        if (ts[counter] > params.timenextbreath + params.IF*60.0/params.RR) {
            params.timenextbreath = params.timenextbreath + 60.0/params.RR;
        }
        ys1[counter + 1] = iterate1(ts[counter], ys1[counter], ys2[counter], ys3[counter], params);
        ys2[counter + 1] = iterate2(ts[counter], ys1[counter], ys2[counter], ys3[counter], params);
        ys3[counter + 1] = iterate3(ts[counter], ys1[counter], ys2[counter], ys3[counter], params);
        counter = counter + 1;
    }

    let trace_u1 = { x: [], y: [], name: 'Pressure'};
    let trace_u2 = { x: [], y: [], name: 'Volume'};
    let trace_u3 = { x: [], y: [], name: 'Flow'};
    
    for (let i = 0; i < N + 1; i++) {
        trace_u1.x.push(ts[i]);  
        trace_u1.y.push(ys1[i]);
        trace_u2.x.push(ts[i]);
        trace_u2.y.push(ys2[i]);
        trace_u3.x.push(ts[i]);
        trace_u3.y.push(ys3[i]);
    }

    const trace1 = {
        x: trace_u1.x,
        y: trace_u1.y,
        type: 'scatter',
        mode: 'lines',
        name: 'Pressure (cm H<sub>2</sub>O)'
    };

    const trace2 = {
        x: trace_u2.x,
        y: trace_u2.y,
        type: 'scatter',
        mode: 'lines',
        xaxis: 'x2',
        yaxis: 'y2',
        name: 'Volume (mL)'
    };

    const trace3 = {
        x: trace_u3.x,
        y: trace_u3.y,
        type: 'scatter',
        mode: 'lines',
        xaxis: 'x3',
        yaxis: 'y3',
        name: 'Flow (mL/s)'
    };

    const layout = {
        grid: {
            rows: 3,
            columns: 1,
            pattern: 'independent',
            roworder: 'top to bottom'
        },
        title: 'Pressure Control'
    };

    Plotly.newPlot('myDiv1', [trace1, trace2, trace3], layout);
}

dfsolve();

function reset() {
    PCslider.value = 120;
    PCnum.value = 12;
    PEEPslider.value = 50;
    PEEPnum.value = 5;
    RRslider.value = 120;
    RRnum.value = 12;
    IFslider.value = 25;
    IFnum.value = 0.25;
    Rslider.value = 30;
    Rnum.value = 0.003;
    Cslider.value = 500;
    Cnum.value = 50;
    tfinalslider.value = 600;
    tfinalnum.value = 60;
    dfsolve();
}

//function is called when slider value changes
PCslider.addEventListener("change", function() {
    PCnum.value = PCslider.value/10;
    dfsolve();    
});

PCnum.addEventListener("change", function() {
    PCslider.value = PCnum.value*10;
    dfsolve();    
});

PEEPslider.addEventListener("change", function() {
    PEEPnum.value = PEEPslider.value/10;
    dfsolve();    
});

PEEPnum.addEventListener("change", function() {
    PEEPslider.value = PEEPnum.value*10;
    dfsolve();    
});

RRslider.addEventListener("change", function() {
    RRnum.value = RRslider.value/10;
    dfsolve();    
});

RRnum.addEventListener("change", function() {
    RRslider.value = RRnum.value*10;
    dfsolve();    
});

IFslider.addEventListener("change", function() {
    IFnum.value = IFslider.value/100;
    dfsolve();    
});

IFnum.addEventListener("change", function() {
    IFslider.value = IFnum.value*100;
    dfsolve();    
});

Rslider.addEventListener("change", function() {
    Rnum.value = Rslider.value/10000;
    dfsolve();    
});

Rnum.addEventListener("change", function() {
    Rslider.value = Rnum.value*10000;
    dfsolve();    
});

Cslider.addEventListener("change", function() {
    Cnum.value = Cslider.value/10;
    dfsolve();    
});

Cnum.addEventListener("change", function() {
    Cslider.value = Cnum.value*10;
    dfsolve();    
});

tfinalslider.addEventListener("change", function() {
    tfinalnum.value = tfinalslider.value/10;
    dfsolve();
});

tfinalnum.addEventListener("change", function() {
    tfinalslider.value = tfinalnum.value*10;
    dfsolve();
});