var Vd1slider = document.getElementById("Vd1slider");
var Vd1html = document.getElementById("Vd1html");
Vd1html.innerHTML = "<i>V</i><sub><i>d</i>, central</sub> (mL/kg)";
var Vd1num = document.getElementById("Vd1");
Vd1num.value = Vd1slider.value/10;

var Clslider = document.getElementById("Clslider");
var Clhtml = document.getElementById("Clhtml");
Clhtml.innerHTML = "<i>C</i><sub><i>l</i></sub> (mL/kg/min)";
var Clnum = document.getElementById("Cl");
Clnum.value = Clslider.value/10;

var bslider = document.getElementById("bslider");
var bhtml = document.getElementById("bhtml");
bhtml.innerHTML = "Bolus Amount (mg/kg)";
var bnum = document.getElementById("b");
bnum.value = bslider.value/10;

var infusionslider = document.getElementById("infusionslider");
var infusionhtml = document.getElementById("infusionhtml");
infusionhtml.innerHTML = "Infusion (mg/kg/min)";
var infusionnum = document.getElementById("infusion");
infusionnum.value = infusionslider.value/10;

var tbolusslider = document.getElementById("tbolusslider");
var tbolushtml = document.getElementById("tbolushtml");
tbolushtml.innerHTML = "Bolus Time (min)";
var tbolusnum = document.getElementById("tbolus");
tbolusnum.value = tbolusslider.value/10;

var tperiodslider = document.getElementById("tperiodslider");
var tperiodhtml = document.getElementById("tperiodhtml");
tperiodhtml.innerHTML = "Bolus Frequency (min)";
var tperiodnum = document.getElementById("tperiod");
tperiodnum.value = tperiodslider.value/10;

var tfinalslider = document.getElementById("tfinalslider");
var tfinalhtml = document.getElementById("tfinalhtml");
tfinalhtml.innerHTML = "<i>t</i><sub>final</sub> (min)";
var tfinalnum = document.getElementById("tfinal");
tfinalnum.value = tfinalslider.value/10;

function dfsolve() {
    // Exactly computes 

    function expn(x) {
        return math.exp(x); //1 + x + x*x/2 + x*x*x/6 + x*x*x*x/24 + x*x*x*x*x/120 + x*x*x*x*x*x/720 + x*x*x*x*x*x*x/5040;
    }

    function iterate(t, C1, params) {
        if (Clslider.value == 0) {
            return C1 + params.dt*params.b/params.Vd1;
        }
        else {
            const a1 = -params.Cl/params.Vd1;
            const e = params.b/params.Vd1;
            const dt = params.dt;
            return (expn(a1*dt)*(e + C1*a1) - e)/a1;
        }
    }
    
    let params = { b: [], Cl: [], Vd1: [], tbolus: [], tperiod: [], tfinal: [], dt: []};
    
    params.Vd1 = Vd1slider.value/10000.0;
    params.Cl = Clslider.value/10000.0;
    params.b = bslider.value/tbolusslider.value;
    params.tbolus = tbolusslider.value/10.0;
    params.tperiod = tperiodslider.value/10.0;
    params.tfinal = tfinalslider.value/10.0;
    params.dt = 0.1;

    let t0 = 0;
    let y01 = 0;

    let N = Math.ceil(params.tfinal/params.dt);
    let Nhalf = Math.ceil(params.tbolus/params.dt);
    let Nperiod = Math.ceil(params.tperiod/params.dt);

    let ts = new Array(N + 1);
    let ys1 = new Array(N + 1);

    for (let i = 0; i < N+1; i++) {
       ts[i] = i*params.dt;
    }

    ys1[0] = y01;

    let counter = 0;
    let counterperiod = 0;
    while (counter < N) {
        if (counterperiod < Nhalf) {
            params.b = bslider.value/tbolusslider.value;
        }
        else {
            params.b = 0;
        }
        ys1[counter + 1] = iterate(ts[counter], ys1[counter], params);
        counter = counter + 1;
        counterperiod = counterperiod + 1;
        if (counterperiod > Nperiod) {
            counterperiod = 0;
        }
    }

    let trace_u1 = { x: [], y: [], name: 'Central Compartment'};
    
    for (let i = 0; i < N + 1; i++) {
        trace_u1.x.push(ts[i]);  
        trace_u1.y.push(ys1[i]);
    }

    var layout1 = {
        title: {
            text:'Central Compartment',
            font: {
                family: 'Courier New, monospace',
                size: 24
            },
            //xref: 'paper',
            //x: 0.05,
        },
        xaxis: {
            title: {
                text: 'Time (min)',
                font: {
                    family: 'Courier New, monospace',
                    size: 18,
                    color: '#7f7f7f'
                }
            },
        },
        yaxis: {
            title: {
                text: 'Concentration (mg/L)',
                font: {
                    family: 'Courier New, monospace',
                    size: 18,
                    color: '#7f7f7f'
                }
            }
        }
    };

    Plotly.newPlot('myDiv1', [trace_u1], layout1);
}

dfsolve();

function reset() {
    Vd1slider.value = 700;
    Vd1num.value = 70;
    Clslider.value = 500;
    Clnum.value = 50;
    bslider.value = 2000;
    bnum.value = 200;
    tbolusslider.value = 200;
    tbolusnum.value = 20;
    infusionslider.value = 10*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10;
    tperiodslider.value = 1000;
    tperiodnum.value = 100;
    tfinalslider.value = 10000;
    tfinalnum.value = 1000;
    dfsolve();
}

function continuous() {
    tperiodnum.value = 0.1;
    tperiodslider.value = 1;
    bslider.value = 10*bslider.value/tbolusslider.value;
    bnum.value = bslider.value/10;
    tbolusslider.value = 10;
    tbolusnum.value = 1;
    infusionslider.value = 10*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10;
    dfsolve();
}

function propofol() {
    Vd1slider.value = 700;
    Vd1num.value = 70;
    Clslider.value = 500;
    Clnum.value = 50;
    infusionslider.value = 10*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10;
    dfsolve();
}

//function is called when slider value changes
Vd1slider.addEventListener("change", function() {
    Vd1num.value = Vd1slider.value/10;
    dfsolve();    
});

Vd1num.addEventListener("change", function() {
    Vd1slider.value = Vd1num.value*10;
    dfsolve();    
});

Clslider.addEventListener("change", function() {
    Clnum.value = Clslider.value/10;
    dfsolve();    
});

Clnum.addEventListener("change", function() {
    Clslider.value = Clnum.value*10;
    dfsolve();    
});

bslider.addEventListener("change", function() {
    bnum.value = bslider.value/10;
    infusionslider.value = 10.0*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10.0;
    dfsolve();    
});

bnum.addEventListener("change", function() {
    bslider.value = bnum.value*10;
    infusionslider.value = 10.0*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10.0;
    dfsolve();    
});

tbolusslider.addEventListener("change", function() {
    tbolusnum.value = tbolusslider.value/10;
    infusionslider.value = 10.0*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10.0;
    dfsolve();    
});

tbolusnum.addEventListener("change", function() {
    tbolusslider.value = tbolusnum.value*10;
    infusionslider.value = 10.0*bslider.value/tbolusslider.value;
    infusionnum.value = infusionslider.value/10.0;
    dfsolve();    
});

infusionslider.addEventListener("change", function() {
    infusionnum.value = infusionslider.value/10;
    bslider.value = infusionslider.value*tbolusslider.value/10.0;
    bnum.value = bslider.value/10;
    dfsolve();
});

infusionnum.addEventListener("change", function() {
    infusionslider.value = infusionnum.value*10;
    bslider.value = infusionslider.value*tbolusslider.value/10.0;
    bnum.value = bslider.value/10;
    dfsolve();
});

tperiodslider.addEventListener("change", function() {
    tperiodnum.value = tperiodslider.value/10;
    dfsolve();    
});

tperiodnum.addEventListener("change", function() {
    tperiodslider.value = tperiodnum.value*10;
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