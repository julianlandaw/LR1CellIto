var inafacnum = document.getElementById("inafac");
var itofacnum = document.getElementById("itofac");
var pclnum = document.getElementById("pcl");
var beatsnum = document.getElementById("beats");
var tauXfacnum = document.getElementById("tauXfac"); //5
var icalfacnum = document.getElementById("icalfac"); //1.15
var ikfacnum = document.getElementById("ikfac"); //1.0
var ikifacnum = document.getElementById("ikifac"); //2.2
var yshiftnum = document.getElementById("yshift");

var Cell = new LR1CellIto();

reset();
pacecell();

function pacecell() {
    Cell.inafac = parseFloat(inafacnum.value);
    Cell.itofac = parseFloat(itofacnum.value);
    Cell.icalfac = parseFloat(icalfacnum.value);
    Cell.ikifac = parseFloat(ikifacnum.value);
    Cell.ikfac = parseFloat(ikfacnum.value);
    Cell.tauXfac = parseFloat(tauXfacnum.value);
    Cell.yshift = parseFloat(yshiftnum.value);
    let pcl = parseFloat(pclnum.value);
    let beats = parseFloat(beatsnum.value);
    for (let i = 0; i < Math.ceil(10000/Cell.dt); i++) {
        Cell.stepdt(0);
    }
    let t = 0;
    let nextstim = 100;
    for (t = Cell.dt; t < pcl*10; t = t + Cell.dt) {
        if (t > nextstim - Cell.dt/2 & t < nextstim + Cell.stimduration - Cell.dt/2) {
            Cell.stepdt(Cell.stimulus);
        }
        else {
            Cell.stepdt(0);
        }
        if (t > nextstim + Cell.stimduration + Cell.dt/2) {
            nextstim = nextstim + pcl;
        }
    }
    let ts = [];
    let vs = [];
    let xs = [];
    let apds = [];
    let xsinit = [];
    let apdcounter = [];
    let apdnum = 0;
    let t0 = 0;
    let v0 = Cell.v;
    let x0 = Cell.xr;
    t = t0;
    ts.push(t0);
    vs.push(v0);
    xs.push(x0);
    nextstim = 100;
    let tsave = 5;
    let oldv = v0;
    let newv = v0;
    let startapdtime = 0;
    let minapd = 0;
    let maxapd = 0;
    let thisapd = 0;
    let minX = 1;
    let maxX = 0;
    let thisX = 0;
    for (t = Cell.dt; t < pcl*beats; t = t + Cell.dt) {
        oldv = Cell.v
        if (t > nextstim - Cell.dt/2 & t < nextstim + Cell.stimduration - Cell.dt/2) {
            Cell.stepdt(Cell.stimulus);
        }
        else {
            Cell.stepdt(0);
        }
        newv = Cell.v
        if (newv > Cell.threshold & oldv < Cell.threshold) {
            startapdtime = t;
            thisX = Cell.xr;
            xsinit.push(thisX);
            if (thisX > maxX) {
                maxX = thisX;
            }
            if (thisX < minX) {
                minX = thisX;
            }
        }
        if (newv < Cell.threshold & oldv > Cell.threshold & startapdtime > 1) {
            thisapd = t - startapdtime;
            apds.push(thisapd);
            apdnum = apdnum + 1;
            apdcounter.push(apdnum);
            if (apdnum == 1) {
                minapd = thisapd;
                maxapd = thisapd;
            }
            if (thisapd < minapd) {
                minapd = thisapd;
            }
            if (thisapd > maxapd) {
                maxapd = thisapd;
            }
        }
        if (t > nextstim + Cell.stimduration + Cell.dt/2) {
            nextstim = nextstim + pcl;
        }
        if (t > tsave - Cell.dt/2) {
            ts.push(t/1000);
            vs.push(Cell.v);
            xs.push(Cell.xr);
            tsave = tsave + 5;
        }
    }
    const trace1 = {
        x: ts,
        y: vs,
        type: 'scatter',
        mode: 'lines',
        name: 'Action Potential'
    };

    const trace2 = {
        x: ts,
        y: xs,
        xlabel: 'Time',
        type: 'scatter',
        mode: 'lines',
        xaxis: 'x2',
        yaxis: 'y2',
        name: 'X'
    };
    const trace3 = {
        x: apdcounter,
        y: apds,
        type: 'scatter',
        mode: 'markers',
        xaxis: 'x3',
        yaxis: 'y3',
        name: 'APDs'
    };

    const trace4 = {
        x: xsinit,
        y: apds,
        type: 'scatter',
        mode: 'markers',
        xaxis: 'x4',
        yaxis: 'y4',
        name: 'APD vs X'
    };
                
    var layout1 = {
        title: {
            text:'The LR1 Model',
            font: {
                family: 'Courier New, monospace',
                size: 24
            },
            //xref: 'paper',
            //x: 0.05,
        },
        xaxis: {
            title: 'Time (s)'
        },
        xaxis2: {
            title: 'Time (s)'
        },
        xaxis3: {
            title: 'Beat Num'
        },
        xaxis4: {
            title: 'X',
            range: [minX - 0.01, maxX + 0.01]
        },
        yaxis: {
            title: 'Voltage (mV)'
        },
        yaxis2: {
            title: 'X'
        },
        yaxis3: {
            title: 'APD (ms)',
            range: [minapd - 100, maxapd + 100]
        },
        yaxis4: {
            title: 'APD (ms)',
            range: [minapd - 100, maxapd + 100]
        },
        showlegend: false,
        grid: {
            rows: 4,
            columns: 1,
            pattern: 'independent',
            roworder: 'top to bottom',
            ygap: 0.5
        },
        height: 600
    };
        
    Plotly.newPlot('myDiv1', [trace1, trace2, trace3, trace4], layout1);
}

function reset() {
    inafacnum.value = 1.0;
    itofacnum.value = 0;
    pclnum.value = 500;
    beatsnum.value = 10;
    tauXfacnum.value = 1;
    icalfacnum.value = 1;
    ikifacnum.value = 1;
    ikfacnum.value = 1;
    yshiftnum.value = 0;
}

function chaos() {
    inafacnum.value = 1.0;
    itofacnum.value = 1.05;
    pclnum.value = 378; //507, 880
    beatsnum.value = 50;
    tauXfacnum.value = 5;
    icalfacnum.value = 1.15;
    ikifacnum.value = 2.2;
    ikfacnum.value = 1;
    yshiftnum.value = 8.0;
}

function eads() {
    inafacnum.value = 1.0;
    itofacnum.value = 0.0;
    icalfacnum.value = 1.0;
    ikfacnum.value = 1.0;
    ikifacnum.value = 1.0;
    tauXfacnum.value = 10.0;
    yshiftnum.value = 0.0;
    beatsnum.value = 10;
    pclnum.value = 1575;
}

/*
inafacnum.addEventListener("change", function() {
    pacecell();
});

itofacnum.addEventListener("change", function() {
    pacecell();
});

pclnum.addEventListener("change", function() {
    pacecell();
});

beatsnum.addEventListener("change", function() {
    pacecell();
});

tauXfacnum.addEventListener("change", function() {
    pacecell();
});

icalfacnum.addEventListener("change", function() {
    pacecell();
});

ikifacnum.addEventListener("change", function() {
    pacecell();
});

ikfacnum.addEventListener("change", function() {
    pacecell();
});

yshiftnum.addEventListener("change", function() {
    pacecell();
});
*/
