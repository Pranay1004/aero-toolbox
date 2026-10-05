// Reentry Aerothermodynamics Calculator - Web Frontend (Tkinter-matched)
const API = '/api';

// Matplotlib-style plot layout matching Tkinter "clam" theme
const PL = {
  paper_bgcolor: '#ffffff',
  plot_bgcolor:  '#ffffff',
  font: {family: 'Segoe UI, Arial, sans-serif', color: '#000000', size: 11},
  margin: {l: 58, r: 22, t: 32, b: 48},
  xaxis: {gridcolor: '#e0e0e0', gridwidth: 1, zerolinecolor: '#cccccc'},
  yaxis: {gridcolor: '#e0e0e0', gridwidth: 1, zerolinecolor: '#cccccc'}
};
const PC = {responsive: true, displayModeBar: true,
  modeBarButtonsToRemove: ['lasso2d','select2d','autoScale2d']};

// Tab switching
document.querySelectorAll('.tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.querySelectorAll('.tab-page').forEach(p => p.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById('page-' + btn.dataset.tab).classList.add('active');
  });
});

// TPS source toggle
document.getElementById('tps-source').addEventListener('change', function() {
  document.getElementById('tps-custom-q').disabled = (this.value !== 'custom');
  document.getElementById('tps-custom-t').disabled = (this.value !== 'custom');
});

// Substructure material auto-fill
const MAT_DB = {
  "Al-7xxx (Default)": {rho:2810,cp:960,T:450},
  "Al-2024": {rho:2780,cp:875,T:450},
  "Al-6061": {rho:2700,cp:896,T:450},
  "Inconel 718": {rho:8190,cp:435,T:950},
  "SS-304": {rho:8000,cp:500,T:1100},
  "Titanium (Ti-6Al-4V)": {rho:4430,cp:526,T:650},
};
document.getElementById('tps-sub-mat').addEventListener('change', function() {
  const m = MAT_DB[this.value];
  if (m) { document.getElementById('tps-sub-rho').value=m.rho; document.getElementById('tps-sub-cp').value=m.cp; document.getElementById('tps-target-T').value=m.T; }
});

function showLoad(id) { document.getElementById(id).innerHTML = '<div style="text-align:center;padding:40px;color:#888">Computing...</div>'; }
function showErr(id, m) { document.getElementById(id).textContent = 'ERROR: ' + m; }

// ================================================================
// Tab 1: Stagnation (4x2 grid — matches Tkinter figsize=(7.5,10))
// ================================================================
async function calcStagnation() {
  showLoad('s-plot'); document.getElementById('s-output').value = 'Computing...';
  try {
    const r = await fetch(API+'/stagnation',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({
      planet:gv('s-planet'), h0_km:gf('s-h0'), V0:gf('s-V0'), gamma0:gf('s-gamma'),
      RN:gf('s-RN'), BC:gf('s-BC'), LD:gf('s-LD'), eps:gf('s-eps'), t_max:gf('s-tmax')
    })});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    const rp = d.report;
    document.getElementById('s-output').value =
      ` STAGNATION HEATING — FULL TRAJECTORY REPORT\n` +
      ` -----------------------------------------------------\n` +
      `  Entry Conditions:  h0 = ${gf('s-h0')} km  V0 = ${gf('s-V0')} m/s  gamma0 = ${gf('s-gamma')}°\n` +
      `  Vehicle Specs:     RN = ${gf('s-RN')} m  BC = ${gf('s-BC')} kg/m²  L/D = ${gf('s-LD')}  eps = ${gf('s-eps')}\n` +
      ` -----------------------------------------------------\n` +
      `  Simulation Results:\n` +
      `    Stop Reason: ${rp.stop_reason}\n` +
      `    Flight Time: ${rp.flight_time.toFixed(1)} s\n` +
      `    End State:   h = ${rp.end_h.toFixed(1)} km, V = ${rp.end_V.toFixed(1)} m/s, Mach = ${rp.end_Mach.toFixed(2)}\n` +
      `    Max Mach:    ${rp.max_mach.toFixed(1)}\n` +
      ` -----------------------------------------------------\n` +
      `  FAY-RIDDELL (Eq. 15.38):\n` +
      `    Peak Heat Flux qt: ${rp.fr_peak_q.toFixed(2)} W/cm²\n` +
      `      at h = ${rp.fr_peak_h.toFixed(1)} km, V = ${rp.fr_peak_V.toFixed(0)} m/s, t = ${rp.fr_peak_t.toFixed(1)} s\n` +
      ` -----------------------------------------------------\n` +
      `  ALLEN-EGGERS (Eq. 15.39):\n` +
      `    Peak Heat Flux qt: ${rp.ae_peak_q.toFixed(2)} W/cm²\n` +
      ` -----------------------------------------------------\n` +
      `  FULL FAY-RIDDELL (High-Fidelity Real Gas Iterative):\n` +
      `    Peak Heat Flux qt: ${rp.full_peak_q.toFixed(2)} W/cm²\n` +
      `    Peak Wall Temp Tw: ${rp.full_peak_Tw.toFixed(0)} K\n` +
      ` -----------------------------------------------------\n`;

    const t=d.t, h=d.h_km, V=d.V;
    const el = document.getElementById('s-plot');

    // Tkinter colors: blue solid, red dashed, green solid
    const cBlu='#0072B2', cRed='#D55E00', cGrn='#009E73', cMag='#CC79A7';

    // 4x2 subplots with independent domains
    const traces = [
      // Row 0: Heat flux vs altitude | Heat flux vs time
      {x:d.q_fr,y:h,name:'Fay-Riddell',xaxis:'x',yaxis:'y',line:{color:cBlu,width:2}},
      {x:d.q_ae,y:h,name:'Allen-Eggers',xaxis:'x',yaxis:'y',line:{color:cRed,width:2,dash:'dash'}},
      {x:d.q_full,y:h,name:'Full F-R',xaxis:'x',yaxis:'y',line:{color:cGrn,width:2}},
      {x:t,y:d.q_fr,xaxis:'x2',yaxis:'y2',line:{color:cBlu,width:2}},
      {x:t,y:d.q_ae,xaxis:'x2',yaxis:'y2',line:{color:cRed,width:2,dash:'dash'}},
      {x:t,y:d.q_full,xaxis:'x2',yaxis:'y2',line:{color:cGrn,width:2}},
      // Row 1: Wall temp vs altitude | Wall temp vs time
      {x:d.Tw_fr,y:h,name:'Fay-Riddell',xaxis:'x3',yaxis:'y3',line:{color:cBlu,width:2}},
      {x:d.Tw_ae,y:h,name:'Allen-Eggers',xaxis:'x3',yaxis:'y3',line:{color:cRed,width:2,dash:'dash'}},
      {x:d.Tw_full,y:h,name:'Full F-R',xaxis:'x3',yaxis:'y3',line:{color:cGrn,width:2}},
      {x:t,y:d.Tw_full,xaxis:'x4',yaxis:'y4',line:{color:cGrn,width:2}},
      {x:t,y:d.Tw_ae,xaxis:'x4',yaxis:'y4',line:{color:cRed,width:2,dash:'dash'}},
      // Row 2: Velocity vs altitude | Mach vs altitude
      {x:V,y:h,xaxis:'x5',yaxis:'y5',line:{color:cBlu,width:2}},
      {x:d.Mach,y:h,xaxis:'x6',yaxis:'y6',line:{color:cMag,width:2}},
      // Row 3: Bar chart (peak q) | Load factor
      {x:d.V_bar,y:d.q_fr_bar,name:'F-R',xaxis:'x7',yaxis:'y7',type:'bar',marker:{color:cBlu}},
      {x:d.V_bar,y:d.q_ae_bar,name:'A-E',xaxis:'x7',yaxis:'y7',type:'bar',marker:{color:cRed}},
      {x:t,y:d.n_z,xaxis:'x8',yaxis:'y8',line:{color:cRed,width:2}},
    ];

    const axBase = {...PL.xaxis, showgrid:true, gridcolor:'#e0e0e0', gridwidth:1};
    const ayBase = {...PL.yaxis, showgrid:true, gridcolor:'#e0e0e0', gridwidth:1};

    Plotly.newPlot(el, traces, {
      ...PL, showlegend:false, height:880,
      grid: {rows:4, columns:2, pattern:'independent', xgap:0.06, ygap:0.06},
      xaxis:  {...axBase, title:{text:'q_t (W/cm²)',font:{size:10}}, domain:[0,0.46], anchor:'y'},
      yaxis:  {...ayBase, title:{text:'h (km)',font:{size:10}}, domain:[0.765,1], anchor:'x'},
      xaxis2: {...axBase, title:{text:'t (s)',font:{size:10}}, domain:[0.54,1], anchor:'y2'},
      yaxis2: {...ayBase, title:{text:'q_t (W/cm²)',font:{size:10}}, domain:[0.765,1], anchor:'x2'},
      xaxis3: {...axBase, title:{text:'T_w (K)',font:{size:10}}, domain:[0,0.46], anchor:'y3'},
      yaxis3: {...ayBase, title:{text:'h (km)',font:{size:10}}, domain:[0.515,0.735], anchor:'x3'},
      xaxis4: {...axBase, title:{text:'t (s)',font:{size:10}}, domain:[0.54,1], anchor:'y4'},
      yaxis4: {...ayBase, title:{text:'T_w (K)',font:{size:10}}, domain:[0.515,0.735], anchor:'x4'},
      xaxis5: {...axBase, title:{text:'V (m/s)',font:{size:10}}, domain:[0,0.46], anchor:'y5'},
      yaxis5: {...ayBase, title:{text:'h (km)',font:{size:10}}, domain:[0.265,0.485], anchor:'x5'},
      xaxis6: {...axBase, title:{text:'Mach',font:{size:10}}, domain:[0.54,1], anchor:'y6'},
      yaxis6: {...ayBase, title:{text:'h (km)',font:{size:10}}, domain:[0.265,0.485], anchor:'x6'},
      xaxis7: {...axBase, title:{text:'Entry velocity (m/s)',font:{size:10}}, domain:[0,0.46], anchor:'y7'},
      yaxis7: {...ayBase, title:{text:'q_t (W/cm²)',font:{size:10}}, domain:[0,0.235], anchor:'x7'},
      xaxis8: {...axBase, title:{text:'t (s)',font:{size:10}}, domain:[0.54,1], anchor:'y8'},
      yaxis8: {...ayBase, title:{text:'Load factor (g)',font:{size:10}}, domain:[0,0.235], anchor:'x8'},
    }, PC);
  } catch(e) { showErr('s-output',e.message); }
}

// ================================================================
// Tab 2: BC Calc
// ================================================================
function calcBC() {
  const m=+document.getElementById('bc-m').value, cd=+document.getElementById('bc-cd').value, s=+document.getElementById('bc-s').value;
  if(m>0&&cd>0&&s>0) document.getElementById('bc-result').textContent = (m/(cd*s)).toFixed(4)+' kg/m²';
  else document.getElementById('bc-result').textContent = 'Enter positive values';
}
function calcBCInv() {
  const mode=document.querySelector('input[name="bc-inv"]:checked').value;
  const bc=+document.getElementById('bc-inv-bc').value, p1=+document.getElementById('bc-inv-p1').value, p2=+document.getElementById('bc-inv-p2').value;
  let r='';
  if(mode==='mass') r='m = '+(bc*p1*p2).toFixed(4)+' kg';
  else if(mode==='cd') r='Cd = '+(p1/(bc*p2)).toFixed(4);
  else r='S = '+(p1/(bc*p2)).toFixed(4)+' m²';
  document.getElementById('bc-inv-result').textContent = r;
}

// ================================================================
// Tab 3: Ballistic (single plot with twin y-axis — matches Tkinter figsize=(6.2,4.6))
// ================================================================
async function calcBallistic() {
  showLoad('bl-plot');
  try {
    const r = await fetch(API+'/ballistic',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({
      Ve:gf('bl-Ve'),BC:gf('bl-BC'),gamma:gf('bl-gamma'),RN:gf('bl-RN')
    })});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    document.getElementById('bl-output').value =
      ` BALLISTIC REENTRY ANALYSIS\n -----------------------------------------------------\n`+
      `  Entry Velocity: ${gf('bl-Ve')} m/s\n  Entry Angle: ${gf('bl-gamma')}°\n  Nose Radius: ${gf('bl-RN')} m\n  BC: ${gf('bl-BC')} kg/m²\n`+
      ` -----------------------------------------------------\n  Total Heat Load Q: ${d.Q.toFixed(2)} J/cm²\n  Peak Heat Rate qs: ${d.qmax.toFixed(2)} W/cm²\n -----------------------------------------------------\n`;

    Plotly.newPlot('bl-plot',[
      {x:d.gs,y:d.Qs,name:'Heat load (J/cm²)',line:{color:'#0072B2',width:2}},
      {x:d.gs,y:d.qms,name:'Peak rate (W/cm²)',yaxis:'y2',line:{color:'#D55E00',width:2,dash:'dash'}},
    ],{
      ...PL, height:380,
      title:{text:'Ballistic body vs entry angle',font:{size:12}},
      xaxis:{...PL.xaxis, title:{text:'Flight-path angle γₑ (deg)',font:{size:11}}},
      yaxis:{...PL.yaxis, title:{text:'Heat load (J/cm²)',font:{size:11,color:'#0072B2'}}, titlefont:{color:'#0072B2'}},
      yaxis2:{...PL.yaxis, title:{text:'Peak heat rate (W/cm²)',font:{size:11,color:'#D55E00'}}, overlaying:'y', side:'right', titlefont:{color:'#D55E00'}},
      showlegend:true, legend:{x:0.02,y:0.98,font:{size:10}},
    }, PC);
  } catch(e) { showErr('bl-output',e.message); }
}

// ================================================================
// Tab 4: Lifting (single plot — matches Tkinter figsize=(6.5,4.4))
// ================================================================
async function calcLifting() {
  showLoad('lt-plot');
  try {
    const angles = document.getElementById('lt-angles').value.trim();
    const body = {planet:gv('lt-planet'),Ve:gf('lt-Ve'),BC:gf('lt-BC'),LD:gf('lt-LD'),gamma0:gf('lt-gamma'),RN:gf('lt-RN')};
    if(angles) body.gamma_angles = angles.split(',').map(Number).filter(x=>!isNaN(x));
    const r = await fetch(API+'/lifting',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    let txt = ` LIFTING REENTRY ANALYSIS (${gv('lt-planet').toUpperCase()})\n -----------------------------------------------------\n`+
      `  Base Heat Load Q: ${d.Q.toFixed(2)} J/cm²\n  Base Peak Heat Rate: ${d.qmax.toFixed(2)} W/cm²\n -----------------------------------------------------\n  VARIATION SUMMARY:\n`;
    d.trajectories.forEach(t => { txt += `  γ = ${t.gamma.toFixed(1)}° | Peak q_conv: ${t.q_peak.toFixed(2)} W/cm² | Duration: ${t.duration.toFixed(1)} s\n`; });
    txt += ' -----------------------------------------------------\n';
    document.getElementById('lt-output').value = txt;

    const colors=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd'];
    Plotly.newPlot('lt-plot',
      d.trajectories.map((tr,i)=>({x:tr.t,y:tr.q_dot_conv,name:`γ = ${tr.gamma.toFixed(1)}°`,line:{color:colors[i%5],width:2}})),
      {
        ...PL, height:380,
        title:{text:'Lifting Reentry: Convective Heat Flux vs Entry Angle ('+gv('lt-planet')+')',font:{size:12}},
        xaxis:{...PL.xaxis, title:{text:'Flight Time (s)',font:{size:11}}},
        yaxis:{...PL.yaxis, title:{text:'Convective Heat Flux q_conv (W/cm²)',font:{size:11}}},
        showlegend:true, legend:{x:0.7,y:0.98,font:{size:9}},
      }, PC);
  } catch(e) { showErr('lt-output',e.message); }
}

// ================================================================
// Tab 5: Trajectory (4x2 grid — matches Tkinter figsize=(7.2,9.0))
// ================================================================
let _lastTraj = null;
async function loadPreset(n) {
  const r = await fetch(API+'/presets'); const p = (await r.json())[n]; if(!p)return;
  sv('tr-planet',p.planet); sv('tr-V0',p.V0); sv('tr-gamma',p.gamma0); sv('tr-h0',p.h0_km);
  sv('tr-RN',p.RN); sv('tr-BC',p.BC); sv('tr-LD',p.LD); sv('tr-bank',p.bank);
  calcTrajectory();
}
async function calcTrajectory() {
  showLoad('tr-plot'); document.getElementById('tr-output').value = 'Running trajectory...';
  try {
    const r = await fetch(API+'/trajectory',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({
      planet:gv('tr-planet'),V0:gf('tr-V0'),gamma0:gf('tr-gamma'),h0_km:gf('tr-h0'),
      RN:gf('tr-RN'),BC:gf('tr-BC'),LD:gf('tr-LD'),bank:gf('tr-bank'),dt:gf('tr-dt'),t_max:gf('tr-tmax')
    })});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    _lastTraj = d;
    const rp=d.report;
    document.getElementById('tr-output').value =
      ` REENTRY TRAJECTORY REPORT\n -----------------------------------------------------\n`+
      `  Initial State: V0 = ${gf('tr-V0')} m/s  gamma0 = ${gf('tr-gamma')}°  h0 = ${gf('tr-h0')} km\n`+
      `  Vehicle:       L/D = ${gf('tr-LD')}  BC = ${gf('tr-BC')} kg/m²  RN = ${gf('tr-RN')} m  sigma = ${gf('tr-bank')}°\n`+
      ` -----------------------------------------------------\n  Simulation Results:\n`+
      `    Stop Reason: ${rp.stop_reason}\n    Flight Time: ${rp.flight_time.toFixed(1)} s\n    Range:       ${rp.range_km.toFixed(1)} km\n`+
      `  Peak Conditions:\n    Peak Heat Rate:     ${rp.peak_q.toFixed(2)} W/cm²  (at t=${rp.peak_q_t.toFixed(1)} s, h=${rp.peak_q_h.toFixed(1)} km)\n`+
      `    Peak Surface Temp:  ${rp.T_surf_max_K.toFixed(1)} K (${rp.T_surf_max_C.toFixed(1)} °C)  [ε=0.85]\n`+
      `    Peak Load Factor:   ${rp.peak_nz.toFixed(3)} g\n    Peak Dyn. Pressure: ${rp.peak_qdyn.toFixed(2)} kPa\n`+
      `    Max Mach:           ${rp.max_mach.toFixed(1)}\n  Integrated Heating:\n    Total Heat Load Q:  ${rp.total_Q.toFixed(2)} J/cm²\n -----------------------------------------------------\n`;

    const t=d.t, h=d.h_km, V=d.V;
    const cBlu='#0072B2', cRed='#D55E00', cGrn='#009E73', cMag='#CC79A7';
    const axB={...PL.xaxis,showgrid:true,gridcolor:'#e0e0e0'};
    const ayB={...PL.yaxis,showgrid:true,gridcolor:'#e0e0e0'};

    Plotly.newPlot('tr-plot',[
      // Row 0: V-h | t-h
      {x:V,y:h,xaxis:'x',yaxis:'y',line:{color:cBlu,width:2}},
      {x:t,y:h,xaxis:'x2',yaxis:'y2',line:{color:cBlu,width:2}},
      // Row 1: Mach-h | t-Mach
      {x:d.Mach,y:h,xaxis:'x3',yaxis:'y3',line:{color:cMag,width:2}},
      {x:t,y:d.Mach,xaxis:'x4',yaxis:'y4',line:{color:cMag,width:2}},
      // Row 2: Heat rate | Load factor
      {x:t,y:d.q_dot_conv,name:'Convective',xaxis:'x5',yaxis:'y5',line:{color:cBlu,width:2}},
      {x:t,y:d.q_dot_rad,name:'Radiative',xaxis:'x5',yaxis:'y5',line:{color:cRed,width:2,dash:'dash'}},
      {x:t,y:d.n_z,xaxis:'x6',yaxis:'y6',line:{color:cRed,width:2}},
      // Row 3: Pressure | Heat load & range (twin y)
      {x:t,y:d.q_dyn.map(v=>v/1000),name:'Dynamic',xaxis:'x7',yaxis:'y7',line:{color:cGrn,width:2}},
      {x:t,y:d.p_atm.map(v=>v/1000),name:'Static',xaxis:'x7',yaxis:'y7',line:{color:cBlu,width:2,dash:'dash'}},
      {x:t,y:d.Q_cumul,name:'Heat load',xaxis:'x8',yaxis:'y8',line:{color:cBlu,width:2}},
      {x:t,y:d.range_km,name:'Range',xaxis:'x8',yaxis:'y9',line:{color:cRed,width:2,dash:'dash'}},
    ],{
      ...PL, showlegend:false, height:880,
      grid:{rows:4,columns:2,pattern:'independent',xgap:0.06,ygap:0.06},
      xaxis:  {...axB,title:{text:'V (m/s)',font:{size:10}},domain:[0,0.46],anchor:'y',autorange:'reversed'},
      yaxis:  {...ayB,title:{text:'h (km)',font:{size:10}},domain:[0.765,1],anchor:'x'},
      xaxis2: {...axB,title:{text:'t (s)',font:{size:10}},domain:[0.54,1],anchor:'y2'},
      yaxis2: {...ayB,title:{text:'h (km)',font:{size:10}},domain:[0.765,1],anchor:'x2'},
      xaxis3: {...axB,title:{text:'Mach',font:{size:10}},domain:[0,0.46],anchor:'y3',autorange:'reversed'},
      yaxis3: {...ayB,title:{text:'h (km)',font:{size:10}},domain:[0.515,0.735],anchor:'x3'},
      xaxis4: {...axB,title:{text:'t (s)',font:{size:10}},domain:[0.54,1],anchor:'y4'},
      yaxis4: {...ayB,title:{text:'Mach',font:{size:10}},domain:[0.515,0.735],anchor:'x4'},
      xaxis5: {...axB,title:{text:'t (s)',font:{size:10}},domain:[0,0.46],anchor:'y5'},
      yaxis5: {...ayB,title:{text:'Heat rate (W/cm²)',font:{size:10}},domain:[0.265,0.485],anchor:'x5'},
      xaxis6: {...axB,title:{text:'t (s)',font:{size:10}},domain:[0.54,1],anchor:'y6'},
      yaxis6: {...ayB,title:{text:'Load factor (g)',font:{size:10}},domain:[0.265,0.485],anchor:'x6'},
      xaxis7: {...axB,title:{text:'t (s)',font:{size:10}},domain:[0,0.46],anchor:'y7'},
      yaxis7: {...ayB,title:{text:'Pressure (kPa)',font:{size:10}},domain:[0,0.235],anchor:'x7'},
      xaxis8: {...axB,title:{text:'t (s)',font:{size:10}},domain:[0.54,1],anchor:'y8'},
      yaxis8: {...ayB,title:{text:'Heat load (J/cm²)',font:{size:10}},domain:[0,0.235],anchor:'x8'},
      yaxis9:{...ayB,title:{text:'Range (km)',font:{size:10}},domain:[0,0.235],anchor:'x8',overlaying:'y8',side:'right'},
    }, PC);
  } catch(e) { showErr('tr-output',e.message); }
}

// ================================================================
// Tab 6: TPS (1x3 grid — matches Tkinter figsize=(8.0,4.2))
// ================================================================
async function calcTPS() {
  showLoad('tps-plot'); document.getElementById('tps-output').value='Computing...';
  try {
    const body = {planet:gv('tps-planet'),eps:gf('tps-eps'),d_mm:gf('tps-d'),rho_tile:gf('tps-rho'),k_tile:gf('tps-k'),
      cp_tile:gf('tps-cp'),Cv:gf('tps-Cv'),hv:gf('tps-hv'),glass_d_um:gf('tps-glass-d'),glass_k:gf('tps-glass-k'),
      sub_d_mm:gf('tps-sub-d'),sub_rho:gf('tps-sub-rho'),sub_cp:gf('tps-sub-cp'),
      target_T:gf('tps-target-T'),use_tdep:document.getElementById('tps-tdep').checked,
      heat_source:document.getElementById('tps-source').value==='custom'?'custom':'trajectory',
      custom_q_peak:gf('tps-custom-q'),custom_duration:gf('tps-custom-t')};
    const r = await fetch(API+'/tps',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    document.getElementById('tps-output').value =
      ` 3-LAYER COMPOSITE TPS TRANSIENT THERMAL REPORT (${gv('tps-planet').toUpperCase()})\n`+
      ` -----------------------------------------------------\n  STATIC MASS & THERMAL ESTIMATES:\n`+
      `  - Heat-Sink Mass: ${d.m_sink.toFixed(4)} g/cm²\n  - Radiative Eq. Temp Tw: ${d.Tw_eq.toFixed(1)} K (${(d.Tw_eq-273.15).toFixed(1)} °C)\n`+
      `  - Glass Coating Status: ${d.glass_status} (${d.T_surf_max_C.toFixed(1)} °C)\n`+
      ` -----------------------------------------------------\n  TRANSIENT THERMAL SOAK RESULTS:\n`+
      `   - Peak Glass Surface T_coat: ${d.T_surface[d.T_surface.length-1]?.toFixed(1)} °C\n`+
      `   - Peak Backwall Metal T: ${d.T_backwall_max.toFixed(1)} K (${(d.T_backwall_max-273.15).toFixed(1)} °C)\n`+
      `   - Thermal Soak Time: ${d.t_soak.toFixed(1)} s\n -----------------------------------------------------\n`;

    const cBlu='#0072B2', cRed='#D55E00', cGrn='#009E73';

    Plotly.newPlot('tps-plot',[
      {x:d.t_bal,y:d.q_bal,name:'Ballistic',line:{color:cRed,width:1.5}},
      {x:d.t_glide,y:d.q_glide,name:'Glide',line:{color:cBlu,width:1.5,dash:'dash'}},
      {x:d.tps_t,y:d.T_surface,name:'Glass Surface',xaxis:'x2',yaxis:'y2',line:{color:cRed,width:1.5}},
      {x:d.tps_t,y:d.T_tile_surface,name:'Tile Surface',xaxis:'x2',yaxis:'y2',line:{color:'#FF8C00',width:1.5,dash:'dash'}},
      {x:d.tps_t,y:d.T_center,name:'Tile Center',xaxis:'x2',yaxis:'y2',line:{color:cGrn,width:1.5,dash:'dot'}},
      {x:d.tps_t,y:d.T_backwall,name:'Backwall',xaxis:'x3',yaxis:'y3',line:{color:cBlu,width:2}},
    ],{
      ...PL, showlegend:true, legend:{x:0.5,y:1.02,xanchor:'center',orientation:'h',font:{size:8}}, height:340,
      grid:{rows:1,columns:3,pattern:'independent',xgap:0.06},
      xaxis:  {...PL.xaxis,title:{text:'Time (s)',font:{size:9}},domain:[0,0.30],anchor:'y'},
      yaxis:  {...PL.yaxis,title:{text:'q_total (W/cm²)',font:{size:9}},anchor:'x'},
      xaxis2: {...PL.xaxis,title:{text:'Time (s)',font:{size:9}},domain:[0.35,0.65],anchor:'y2'},
      yaxis2: {...PL.yaxis,title:{text:'Temp (°C)',font:{size:9}},anchor:'x2'},
      xaxis3: {...PL.xaxis,title:{text:'Time (s)',font:{size:9}},domain:[0.70,1],anchor:'y3'},
      yaxis3: {...PL.yaxis,title:{text:'Temp (°C)',font:{size:9}},anchor:'x3'},
    }, PC);
  } catch(e) { showErr('tps-output',e.message); }
}

async function solveTPS() {
  showLoad('tps-plot'); document.getElementById('tps-output').value='Solving tile thickness...';
  try {
    const body = {planet:gv('tps-planet'),eps:gf('tps-eps'),rho_tile:gf('tps-rho'),k_tile:gf('tps-k'),cp_tile:gf('tps-cp'),
      glass_d_um:gf('tps-glass-d'),glass_k:gf('tps-glass-k'),sub_d_mm:gf('tps-sub-d'),sub_rho:gf('tps-sub-rho'),
      sub_cp:gf('tps-sub-cp'),target_T:gf('tps-target-T'),use_tdep:document.getElementById('tps-tdep').checked,
      method:document.getElementById('tps-method').value};
    const r = await fetch(API+'/tps/solve',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    document.getElementById('tps-output').value =
      ` OPTIMAL TILE SIZING RESULTS\n -----------------------------------------------------\n`+
      `  Attainable: ${d.is_attainable?'YES':'NO'}\n  Required Tile Thickness: ${d.d_req_mm.toFixed(1)} mm\n`+
      `  Peak Backwall Temp: ${d.T_backwall_max.toFixed(1)} K (${(d.T_backwall_max-273.15).toFixed(1)} °C)\n`+
      `  Tile Areal Mass: ${d.mass_areal_g_cm2.toFixed(4)} g/cm²\n -----------------------------------------------------\n`;

    Plotly.newPlot('tps-plot',[
      {x:d.tps_t,y:d.T_surface,name:'Glass Surface',line:{color:'#D55E00',width:2}},
      {x:d.tps_t,y:d.T_backwall,name:'Backwall',line:{color:'#0072B2',width:2}},
    ],{
      ...PL, height:340,
      title:{text:'Sized TPS Thermal Response',font:{size:12}},
      xaxis:{...PL.xaxis,title:{text:'Time (s)',font:{size:11}}},
      yaxis:{...PL.yaxis,title:{text:'Temp (°C)',font:{size:11}}},
      showlegend:true, legend:{font:{size:10}},
    }, PC);
  } catch(e) { showErr('tps-output',e.message); }
}
function solveTPSTarget() { openSolver('temp'); }

// ================================================================
// Tab 7: Parametric Sweep (2x2 grid — matches Tkinter figsize=(8.5,7.5))
// ================================================================
async function calcSinglePoint() {
  showLoad('sw-plot'); document.getElementById('sw-output').value='Computing...';
  try {
    const body = getSweepBody();
    const r = await fetch(API+'/parametric-sweep',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    const sp=d.single_point,traj=d.traj;
    document.getElementById('sw-output').value =
      ` AEROTHERMAL COURSE — TRAJECTORY & HEATING REPORT\n -----------------------------------------------------\n`+
      `  Entry: V_atm = ${gv('sw-Vatm')} m/s  gamma_e = ${gv('sw-gamma')}°\n`+
      `  Vehicle: R_n = ${gv('sw-RN')} m  BC = ${gv('sw-BC')} kg/m²  L/D = ${gv('sw-LD')}\n`+
      ` -----------------------------------------------------\n`+
      `  PARAMETER COMPARISON: ANALYTICAL vs NUMERICAL\n`+
      `  Peak Heat Flux [W/cm²]:   Analytical=${sp.q_max_an.toFixed(2)}  Numerical=${sp.q_max_num.toFixed(2)}\n`+
      `  Total Heat Load [J/cm²]:  Analytical=${sp.Q_an.toFixed(2)}  Numerical=${sp.Q_num.toFixed(2)}\n`+
      `  Distributed q at θ=${gv('sw-theta')}°: ${sp.q_local.toFixed(2)} W/cm²\n`+
      `  Hemisphere Total: ${sp.Q_hemi.toFixed(0)} W\n -----------------------------------------------------\n`;

    const t=traj.t, h=traj.h_km, V=traj.V;
    const axB={...PL.xaxis,showgrid:true,gridcolor:'#e0e0e0'};
    const ayB={...PL.yaxis,showgrid:true,gridcolor:'#e0e0e0'};

    Plotly.newPlot('sw-plot',[
      {x:t,y:traj.q_dot_conv,name:'Convective',xaxis:'x',yaxis:'y',line:{color:'#0072B2',width:2}},
      {x:t,y:traj.q_dot_rad,name:'Radiative',xaxis:'x',yaxis:'y',line:{color:'#D55E00',width:2,dash:'dash'}},
      {x:t,y:h,xaxis:'x2',yaxis:'y2',line:{color:'#0072B2',width:2}},
      {x:traj.Q_cumul,y:h,xaxis:'x3',yaxis:'y3',line:{color:'#D55E00',width:2}},
      {x:V,y:h,xaxis:'x4',yaxis:'y4',line:{color:'#009E73',width:2}},
    ],{
      ...PL, showlegend:true, legend:{font:{size:9}}, height:650,
      grid:{rows:2,columns:2,pattern:'independent',xgap:0.06,ygap:0.06},
      xaxis:  {...axB,title:{text:'Time (s)',font:{size:10}},domain:[0,0.46],anchor:'y'},
      yaxis:  {...ayB,title:{text:'Stagnation heat rate (W/cm²)',font:{size:10}},domain:[0.53,1],anchor:'x'},
      xaxis2: {...axB,title:{text:'Time (s)',font:{size:10}},domain:[0.54,1],anchor:'y2'},
      yaxis2: {...ayB,title:{text:'Altitude (km)',font:{size:10}},domain:[0.53,1],anchor:'x2'},
      xaxis3: {...axB,title:{text:'Accumulated Heat Load (J/cm²)',font:{size:10}},domain:[0,0.46],anchor:'y3'},
      yaxis3: {...ayB,title:{text:'Altitude (km)',font:{size:10}},domain:[0,0.47],anchor:'x3'},
      xaxis4: {...axB,title:{text:'Velocity (m/s)',font:{size:10}},domain:[0.54,1],anchor:'y4'},
      yaxis4: {...ayB,title:{text:'Altitude (km)',font:{size:10}},domain:[0,0.47],anchor:'x4',autorange:'reversed'},
    }, PC);
  } catch(e) { showErr('sw-output',e.message); }
}

async function calcSweep() {
  showLoad('sw-plot'); document.getElementById('sw-output').value='Running parametric sweep...';
  try {
    const body = getSweepBody();
    const r = await fetch(API+'/parametric-sweep',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    const sv = document.querySelector('input[name="sw-var"]:checked').value;
    const varLabel = sv==='gamma'?'FPA (deg)':sv==='BC'?'BC (kg/m²)':'V_atm (m/s)';
    let txt = `==========================================================\n PARAMETRIC SWEEP RESULTS (Sweep: ${varLabel})\n==========================================================\n`;
    txt += ` ${varLabel.padStart(12)}  ${'Max q_dot'.padStart(12)}  ${'Total Q'.padStart(12)}  ${'Time'.padStart(10)}\n`;
    txt += ` ${''.padStart(12)}  ${(  'W/cm²').padStart(12)}  ${(  'J/cm²').padStart(12)}  ${'(s)'.padStart(10)}\n`;
    txt += ' '+'-'.repeat(12)+'  '+'-'.repeat(12)+'  '+'-'.repeat(12)+'  '+'-'.repeat(10)+'\n';
    d.sweep_vals.forEach((v,i) => { txt += ` ${v.toFixed(2).padStart(12)}  ${d.q_maxs[i].toFixed(2).padStart(12)}  ${d.Q_totals[i].toFixed(2).padStart(12)}  ${d.flight_times[i].toFixed(1).padStart(10)}\n`; });
    txt += '==========================================================\n';
    document.getElementById('sw-output').value = txt;

    const axB={...PL.xaxis,showgrid:true,gridcolor:'#e0e0e0'};
    const ayB={...PL.yaxis,showgrid:true,gridcolor:'#e0e0e0'};

    Plotly.newPlot('sw-plot',[
      {x:d.sweep_vals,y:d.q_maxs,name:'Peak q',xaxis:'x',yaxis:'y',line:{color:'#0072B2',width:2,shape:'spline'},mode:'lines+markers',marker:{size:5}},
      {x:d.sweep_vals,y:d.h_peaks,name:'Peak alt',xaxis:'x2',yaxis:'y2',line:{color:'#CC79A7',width:2,shape:'spline'},mode:'lines+markers',marker:{size:5}},
      {x:d.sweep_vals,y:d.Q_totals,name:'Total Q',xaxis:'x3',yaxis:'y3',line:{color:'#D55E00',width:2,shape:'spline'},mode:'lines+markers',marker:{size:5}},
      {x:d.sweep_vals,y:d.flight_times,name:'Time',xaxis:'x4',yaxis:'y4',line:{color:'#009E73',width:2,shape:'spline'},mode:'lines+markers',marker:{size:5}},
    ],{
      ...PL, showlegend:false, height:650,
      grid:{rows:2,columns:2,pattern:'independent',xgap:0.06,ygap:0.06},
      xaxis:  {...axB,title:{text:varLabel,font:{size:10}},domain:[0,0.46],anchor:'y'},
      yaxis:  {...ayB,title:{text:'Peak Heat Flux (W/cm²)',font:{size:10}},domain:[0.53,1],anchor:'x'},
      xaxis2: {...axB,title:{text:varLabel,font:{size:10}},domain:[0.54,1],anchor:'y2'},
      yaxis2: {...ayB,title:{text:'Peak Heating Altitude (km)',font:{size:10}},domain:[0.53,1],anchor:'x2'},
      xaxis3: {...axB,title:{text:varLabel,font:{size:10}},domain:[0,0.46],anchor:'y3'},
      yaxis3: {...ayB,title:{text:'Total Heat Load (J/cm²)',font:{size:10}},domain:[0,0.47],anchor:'x3'},
      xaxis4: {...axB,title:{text:varLabel,font:{size:10}},domain:[0.54,1],anchor:'y4'},
      yaxis4: {...ayB,title:{text:'Reentry Time (s)',font:{size:10}},domain:[0,0.47],anchor:'x4'},
    }, PC);
  } catch(e) { showErr('sw-output',e.message); }
}

function getSweepBody() {
  return {Vatm:gf('sw-Vatm'),gamma:gf('sw-gamma'),BC:gf('sw-BC'),RN:gf('sw-RN'),LD:gf('sw-LD'),
    eps:gf('sw-eps'),k_sg:gf('sw-ksg'),theta:gf('sw-theta'),
    sweep_var:document.querySelector('input[name="sw-var"]:checked').value,
    sweep_start:gf('sw-start'),sweep_end:gf('sw-end'),sweep_steps:gf('sw-steps')};
}

// ================================================================
// Tab 8: Corridor (single plot — matches Tkinter figsize=(6.5,5.2))
// ================================================================
async function calcCorridor() {
  showLoad('co-plot');
  try {
    const r = await fetch(API+'/corridor',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({
      RN:gf('co-RN'),BC:gf('co-BC'),LD:gf('co-LD'),qtmax:gf('co-qtmax'),nmax:gf('co-nmax'),qmax:gf('co-qmax')
    })});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    const wmax = Math.max(...d.corridor_width.filter(x=>isFinite(x)));
    document.getElementById('co-output').value =
      `==============================================================\n REENTRY CORRIDOR\n==============================================================\n`+
      ` Rₙ=${gv('co-RN')}m  BC=${gv('co-BC')}kg/m²  L/D=${gv('co-LD')}\n q̇_max=${gv('co-qtmax')} W/cm²  n_max=${gv('co-nmax')}g  q_max=${gv('co-qmax')} Pa\n`+
      `==============================================================\n ${'V'.padStart(8)}  ${'h_upper'.padStart(10)}  ${'h_lower'.padStart(10)}  ${'width'.padStart(10)}  status\n`+
      ` ${'m/s'.padStart(8)}  ${'km'.padStart(10)}  ${'km'.padStart(10)}  ${'km'.padStart(10)}\n`;
    for(let i=0;i<d.V.length;i+=40){
      const hu=d.h_upper[i],hl=d.h_lower[i],w=d.corridor_width[i];
      document.getElementById('co-output').value +=
        ` ${d.V[i].toFixed(0).padStart(8)}  ${(hu?.toFixed(2)||'NaN').padStart(10)}  ${(hl?.toFixed(2)||'NaN').padStart(10)}  ${(w?.toFixed(2)||'NaN').padStart(10)}  ${(isFinite(w)&&w>0?'OK':'NO').padStart(6)}\n`;
    }
    document.getElementById('co-output').value += `==============================================================\n Max corridor width: ${wmax.toFixed(2)} km\n==============================================================\n`;

    Plotly.newPlot('co-plot',[
      {x:d.V,y:d.h_heat,name:`Heat flux (${gv('co-qtmax')} W/cm²)`,line:{color:'#D55E00',width:2,dash:'dash'}},
      {x:d.V,y:d.h_load,name:`Load (${gv('co-nmax')} g)`,line:{color:'#E69F00',width:2,dash:'dot'}},
      {x:d.V,y:d.h_q,name:`Dyn. press. (${(gv('co-qmax')/1000).toFixed(0)} kPa)`,line:{color:'#009E73',width:2,dash:'dot'}},
      {x:d.V,y:d.h_upper,name:'Eq. glide (σ=0°)',line:{color:'#666666',width:1.5}},
    ],{
      ...PL, height:440,
      title:{text:'Entry corridor — altitude vs velocity',font:{size:12}},
      xaxis:{...PL.xaxis,title:{text:'Relative velocity (m/s)',font:{size:11}},autorange:'reversed'},
      yaxis:{...PL.yaxis,title:{text:'Altitude (km)',font:{size:11}}},
      showlegend:true, legend:{font:{size:9}},
    }, PC);
  } catch(e) { showErr('co-output',e.message); }
}

async function calcCorridorDrag() {
  showLoad('co-plot');
  try {
    const r = await fetch(API+'/corridor',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({
      RN:gf('co-RN'),BC:gf('co-BC'),LD:gf('co-LD'),qtmax:gf('co-qtmax'),nmax:gf('co-nmax'),qmax:gf('co-qmax')
    })});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    document.getElementById('co-output').value = 'Drag acceleration plot displayed.\n';

    Plotly.newPlot('co-plot',[
      {x:d.V,y:d.D_heat,name:`Heat flux (${gv('co-qtmax')} W/cm²)`,line:{color:'#D55E00',width:2,dash:'dash'}},
      {x:d.V,y:d.D_load,name:`Load (${gv('co-nmax')} g)`,line:{color:'#E69F00',width:2,dash:'dot'}},
      {x:d.V,y:d.D_q,name:`Dyn. press. (${(gv('co-qmax')/1000).toFixed(0)} kPa)`,line:{color:'#009E73',width:2,dash:'dot'}},
      {x:d.V,y:d.D_u,name:'Eq. glide (σ=0°)',line:{color:'#666666',width:1.5}},
      {x:d.V,y:d.D_des,name:'Desired drag',line:{color:'#0072B2',width:3}},
    ],{
      ...PL, height:440,
      title:{text:'Entry corridor — drag acceleration',font:{size:12}},
      xaxis:{...PL.xaxis,title:{text:'Relative velocity (m/s)',font:{size:11}},autorange:'reversed'},
      yaxis:{...PL.yaxis,title:{text:'Drag acceleration (m/s²)',font:{size:11}},range:[0,28]},
      showlegend:true, legend:{font:{size:9}},
    }, PC);
  } catch(e) { showErr('co-output',e.message); }
}
function calcCorridorAlt() { calcCorridor(); }

// ================================================================
// Tab 9: OpenFOAM / PATO
// ================================================================
let _lastOFText = '';
async function calcOpenFOAM() {
  try {
    const r = await fetch(API+'/openfoam',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({
      V0:gf('tr-V0'),gamma0:gf('tr-gamma'),h0_km:gf('tr-h0'),RN:gf('tr-RN'),BC:gf('tr-BC'),
      LD:gf('tr-LD'),bank:gf('tr-bank'),dt:gf('tr-dt'),t_max:gf('tr-tmax'),
      of_RN:gf('of-RN'),eps:gf('of-eps'),title:document.getElementById('of-title').value,dt_out:gf('of-dt')
    })});
    const d = await r.json(); if(d.detail) throw new Error(d.detail);
    _lastOFText = d.full_text;
    document.getElementById('of-output').value = d.full_text;
  } catch(e) { document.getElementById('of-output').value = 'ERROR: '+e.message; }
}
function downloadOF() {
  if(!_lastOFText){alert('Generate first.');return;}
  const b=new Blob([_lastOFText],{type:'text/plain'}),a=document.createElement('a');
  a.href=URL.createObjectURL(b);a.download='boundaryConditions.txt';a.click();
}

// ================================================================
// Solver Modal
// ================================================================
let _svType='';
function openSolver(type) {
  _svType=type;
  const modal=document.getElementById('solver-modal'),title=document.getElementById('sv-title'),form=document.getElementById('sv-form');
  document.getElementById('sv-output').value=''; document.getElementById('sv-plot').innerHTML='';
  if(type==='duration'){
    title.textContent='Target Flight Duration Solver';
    form.innerHTML=`<span class="form-label">Target Duration (s):</span><input id="sv-target" class="form-entry" type="text" value="1200">`+
      `<span class="form-label">Parameter to Vary:</span><select id="sv-vary" class="form-combo"><option value="gamma">Flight Path Angle (γ₀)</option><option value="BC">Ballistic Coeff (BC)</option><option value="L_D">Lift-to-Drag (L/D)</option><option value="bank">Bank Angle (σ)</option></select>`+
      `<span class="form-label">Optimizer Method:</span><select id="sv-method" class="form-combo"><option>bounded</option><option>golden</option><option>brent</option></select>`;
  } else if(type==='heating'){
    title.textContent='Target Heating Conditions Solver';
    form.innerHTML=`<span class="form-label">Target Type:</span><select id="sv-type" class="form-combo"><option value="q_max">Peak Heat Flux (W/cm²)</option><option value="Q_total">Total Heat Load (J/cm²)</option></select>`+
      `<span class="form-label">Target Value:</span><input id="sv-target" class="form-entry" type="text" value="100.0">`+
      `<span class="form-label">Parameter to Vary:</span><select id="sv-vary" class="form-combo"><option value="gamma">Flight Path Angle (γ₀)</option><option value="BC">Ballistic Coeff (BC)</option><option value="L_D">Lift-to-Drag (L/D)</option><option value="bank">Bank Angle (σ)</option></select>`+
      `<span class="form-label">Optimizer Method:</span><select id="sv-method" class="form-combo"><option>bounded</option><option>golden</option><option>brent</option></select>`;
  } else {
    title.textContent='Target Surface Temperature Solver (T_surf)';
    form.innerHTML=`<span class="form-label">Target T_surf (K):</span><input id="sv-target" class="form-entry" type="text" value="1533.15">`+
      `<span class="form-label">Emissivity ε:</span><input id="sv-eps" class="form-entry" type="text" value="0.85">`+
      `<span class="form-label">Solver Method:</span><select id="sv-method" class="form-combo"><option>bisect</option><option>brentq</option><option>secant</option></select>`;
  }
  modal.style.display='block';
}
function closeSolver() { document.getElementById('solver-modal').style.display='none'; }

async function runSolver() {
  const out=document.getElementById('sv-output'),pd=document.getElementById('sv-plot');
  out.value='Solving...'; pd.innerHTML='';
  try {
    const base={V0:gf('tr-V0'),gamma0:gf('tr-gamma'),h0_km:gf('tr-h0'),RN:gf('tr-RN'),BC:gf('tr-BC'),LD:gf('tr-LD'),bank:gf('tr-bank')};
    let url,body;
    if(_svType==='duration'){
      url=API+'/trajectory/solve-duration';
      body={...base,target_t:gv('sv-target'),vary:document.getElementById('sv-vary').value,method:document.getElementById('sv-method').value};
    } else if(_svType==='heating'){
      url=API+'/trajectory/solve-heating';
      body={...base,target_val:+document.getElementById('sv-target').value,target_type:document.getElementById('sv-type').value,
        vary:document.getElementById('sv-vary').value,method:document.getElementById('sv-method').value};
    } else {
      url=API+'/trajectory/solve-temp';
      body={...base,T_surf_target_K:+document.getElementById('sv-target').value,eps:gf('sv-eps'),
        method:document.getElementById('sv-method').value,planet:gv('tr-planet')};
    }
    const r=await fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d=await r.json(); if(d.detail) throw new Error(d.detail);
    const st=d.is_attainable?'SUCCESSFULLY MATCHED':'UNATTAINABLE TARGET';
    if(_svType==='temp'){
      out.value=` TARGET SURFACE TEMPERATURE SOLVER RESULTS (${gv('tr-planet').toUpperCase()})\n =========================================================\n`+
        `  Status: ${st}\n  Target: ${d.T_surf_target_K?.toFixed(1)} K (${d.T_surf_target_C?.toFixed(1)} °C)\n`+
        `  Parameter: ${d.vary_param?.toUpperCase()} = ${d.matched_value?.toFixed(4)}\n  Achieved: ${d.T_surf_peak_K?.toFixed(1)} K (${d.T_surf_peak_C?.toFixed(1)} °C)\n =========================================================\n`;
    } else {
      out.value=` SOLVER RESULTS\n =========================================================\n  Status: ${st}\n`+
        `  Matched Value: ${d.matched_value?.toFixed(4)}\n  Achieved: ${d.achieved_t_s?.toFixed(1)||d.achieved_val?.toFixed(2)}\n`+
        `  Error: ${d.error_s?.toFixed(1)||d.error_val?.toFixed(2)}\n =========================================================\n`;
    }
    if(d.trajectory){
      const tr=d.trajectory;
      Plotly.newPlot('sv-plot',[{x:tr.t,y:tr.q_dot_conv,name:'Conv',line:{color:'#0072B2',width:2}},
        {x:tr.t,y:tr.q_dot_rad,name:'Rad',line:{color:'#D55E00',width:2,dash:'dash'}}],
        {...PL,height:200,title:{text:'Solved Trajectory',font:{size:10}},
         xaxis:{...PL.xaxis,title:{text:'t (s)',font:{size:9}}},yaxis:{...PL.yaxis,title:{text:'q (W/cm²)',font:{size:9}}},
         showlegend:true,legend:{font:{size:8}}},PC);
    }
  } catch(e) { out.value='ERROR: '+e.message; }
}

// ================================================================
// Monte Carlo
// ================================================================
function runMC() {
  _svType='mc';
  const modal=document.getElementById('solver-modal');
  document.getElementById('sv-title').textContent='Monte Carlo Stochastic Simulation (3σ)';
  document.getElementById('sv-form').innerHTML=
    `<span class="form-label">Number of Runs:</span><input id="sv-runs" class="form-entry" type="text" value="100">`+
    `<span class="form-label">Std. Dev. γ₀ (deg):</span><input id="sv-gstd" class="form-entry" type="text" value="0.3">`+
    `<span class="form-label">Std. Dev. BC (kg/m²):</span><input id="sv-bcstd" class="form-entry" type="text" value="15.0">`+
    `<span class="form-label">Std. Dev. L/D:</span><input id="sv-ldstd" class="form-entry" type="text" value="0.02">`+
    `<span class="form-label">Distribution:</span><select id="sv-dist" class="form-combo"><option value="normal">Normal</option><option value="uniform">Uniform</option></select>`;
  document.getElementById('sv-output').value='';document.getElementById('sv-plot').innerHTML='';
  modal.style.display='block';
}
const _origRun = runSolver;
window.runSolver = async function() {
  if(_svType!=='mc') return _origRun();
  const out=document.getElementById('sv-output'),pd=document.getElementById('sv-plot');
  out.value='Running Monte Carlo...'; pd.innerHTML='';
  try {
    const body={V0:gf('tr-V0'),gamma0:gf('tr-gamma'),h0_km:gf('tr-h0'),RN:gf('tr-RN'),BC:gf('tr-BC'),LD:gf('tr-LD'),
      n_runs:+document.getElementById('sv-runs').value, gamma_std:+document.getElementById('sv-gstd').value,
      BC_std:+document.getElementById('sv-bcstd').value, LD_std:+document.getElementById('sv-ldstd').value,
      distribution:document.getElementById('sv-dist').value};
    const r=await fetch(API+'/trajectory/monte-carlo',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
    const d=await r.json(); if(d.detail) throw new Error(d.detail);
    out.value=
      ` MONTE CARLO STOCHASTIC DISPERSION ANALYSIS\n =========================================================\n`+
      `  Successful Runs: ${d.n_runs}\n =========================================================\n`+
      `  ${'Metric'.padEnd(20)}  ${'Mean'.padStart(10)}  ${'Std'.padStart(10)}  ${'Min'.padStart(10)}  ${'Max'.padStart(10)}\n`+
      `  ${'Flight Time (s)'.padEnd(20)}  ${d.flight_time.mean.toFixed(1).padStart(10)}  ${d.flight_time.std.toFixed(1).padStart(10)}  ${d.flight_time.min.toFixed(1).padStart(10)}  ${d.flight_time.max.toFixed(1).padStart(10)}\n`+
      `  ${'Peak q (W/cm²)'.padEnd(20)}  ${d.peak_q.mean.toFixed(2).padStart(10)}  ${d.peak_q.std.toFixed(2).padStart(10)}  ${d.peak_q.min.toFixed(2).padStart(10)}  ${d.peak_q.max.toFixed(2).padStart(10)}\n`+
      `  ${'Total Q (J/cm²)'.padEnd(20)}  ${d.total_Q.mean.toFixed(1).padStart(10)}  ${d.total_Q.std.toFixed(1).padStart(10)}  ${d.total_Q.min.toFixed(1).padStart(10)}  ${d.total_Q.max.toFixed(1).padStart(10)}\n`+
      `  ${'Peak Load (g)'.padEnd(20)}  ${d.peak_g.mean.toFixed(2).padStart(10)}  ${d.peak_g.std.toFixed(2).padStart(10)}  ${d.peak_g.min.toFixed(2).padStart(10)}  ${d.peak_g.max.toFixed(2).padStart(10)}\n`+
      `  ${'Range (km)'.padEnd(20)}  ${d.range_km.mean.toFixed(1).padStart(10)}  ${d.range_km.std.toFixed(1).padStart(10)}  ${d.range_km.min.toFixed(1).padStart(10)}  ${d.range_km.max.toFixed(1).padStart(10)}\n`+
      ` =========================================================\n`+
      `  3σ CONFIDENCE BOUNDS:\n`+
      `    Flight Time:  ${(d.flight_time.mean-3*d.flight_time.std).toFixed(1)} s  to  ${(d.flight_time.mean+3*d.flight_time.std).toFixed(1)} s\n`+
      `    Peak q:       ${Math.max(d.peak_q.mean-3*d.peak_q.std,0).toFixed(2)}  to  ${(d.peak_q.mean+3*d.peak_q.std).toFixed(2)} W/cm²\n`+
      `    Total Q:      ${Math.max(d.total_Q.mean-3*d.total_Q.std,0).toFixed(1)}  to  ${(d.total_Q.mean+3*d.total_Q.std).toFixed(1)} J/cm²\n`+
      `    Peak g-load:  ${Math.max(d.peak_g.mean-3*d.peak_g.std,0).toFixed(2)}  to  ${(d.peak_g.mean+3*d.peak_g.std).toFixed(2)} g\n`+
      `    Range:        ${Math.max(d.range_km.mean-3*d.range_km.std,0).toFixed(1)}  to  ${(d.range_km.mean+3*d.range_km.std).toFixed(1)} km\n`+
      ` =========================================================\n`;
  } catch(e) { out.value='ERROR: '+e.message; }
};

function pushToTPS() { document.querySelector('[data-tab="tps"]').click(); calcTPS(); }
function pushToOF() { document.querySelector('[data-tab="openfoam"]').click(); calcOpenFOAM(); }

// ================================================================
// Reference Tab
// ================================================================
document.getElementById('ref-content').innerHTML = `<h2 style="color:#000">REENTRY AEROTHERMODYNAMICS — EQUATION REFERENCE</h2>
<p>Source: Re-entry Missions, Chapter 15, Sections 15.4.4 &amp; 15.5</p>
<h3>HEAT RATE &amp; HEAT LOAD DEFINITIONS</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.28: q_dot = dQ/dt — Instantaneous heat-transfer rate</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.29: q_dot = f · ρ₀ · σ · C_D · S · V³ / 2 — Energy-conversion form</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.30: q_dot ~ σ · V³ — Scaling form</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.31: σ = ρ_inf / ρ₀ — Density ratio</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.32: Q = ∫q_dot dt = q_dot_ave · Δt — Total heat load</div>
<h3>SCALING LAWS (Eqs. 15.33 – 15.37)</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.33: Q ~ Vₑ² · √(BC / sin(γₑ)) — Ballistic heat load</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.35: q_dot_max ~ V³ · √(BC · sin(γₑ)) — Ballistic peak heat rate</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.36: Q ~ Vₑ² · √((L/D) · BC / sin(γₑ)) — Lifting heat load</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.37: q_dot_max ~ V³ · √(BC · sin(γₑ) / (L/D)) — Lifting peak heat rate</div>
<h3>STAGNATION-POINT HEATING (Eqs. 15.38 – 15.40)</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.38: q_t = 18300 · √(ρ_inf) / √(Rₙ) · (V_inf/1e4)^3.05 — Fay-Riddell [W/cm²]</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.39: q_t = (11030/√Rₙ) · (ρ/ρ₀)^0.5 · (V/V₀)^3.15 — Allen-Eggers</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.40: q_r = 100 · Rₙ · (V/1e4)^8.5 · (ρ/ρ₀)^1.6 — Martin radiative</div>
<h3>TPS SIZING (Eqs. 15.41 – 15.43)</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.41: m = Q / (C_v · ΔT) — Heat-sink TPS mass</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.42: dm = Q / h_v — Ablative TPS mass loss</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.43: T_w = (q_in/(ε·σ_SB))^0.25 — Radiative equilibrium temp</div>
<h3>REENTRY CORRIDOR (Eqs. 15.47 – 15.50)</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.47: (L/m)cos(σ) ≤ [g − V²/r]cos(γ) — Equilibrium-glide constraint</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.48: q_t = (11030/√Rₙ)(ρ/ρ₀)^0.5(V/V_cir)^3.15 ≤ q_t_max</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.49: n_z = (D/m)√(1+(L/D)²)/g ≤ n_max</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Eq. 15.50: q = ½ρV² ≤ q_max</div>
<h3>PHYSICAL CONSTANTS</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">ρ₀ = 1.225 kg/m³  |  V₀ = 7905 m/s  |  g₀ = 9.80665 m/s²  |  σ_SB = 5.670374e-8 W/m²-K⁴</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">R_earth = 6.378e6 m  |  H_scale = 7200 m  |  R_air = 287.05 J/(kg·K)</div>
<h3>VEHICLE PRESETS</h3>
<table style="width:100%;border-collapse:collapse;font-size:9pt">
<tr style="background:#eee"><th style="border:1px solid #ccc;padding:4px">Class</th><th>Rₙ(m)</th><th>BC(kg/m²)</th><th>L/D</th><th>γₑ(deg)</th><th>Vₑ(m/s)</th></tr>
<tr><td style="border:1px solid #ccc;padding:4px">Ballistic</td><td>1.0</td><td>400</td><td>0.0</td><td>20.0</td><td>11000</td></tr>
<tr><td style="border:1px solid #ccc;padding:4px">Lifting</td><td>1.0</td><td>150</td><td>0.3</td><td>1.5</td><td>7900</td></tr>
<tr><td style="border:1px solid #ccc;padding:4px">Winged</td><td>0.3</td><td>100</td><td>1.2</td><td>1.0</td><td>7800</td></tr>
</table>
<h3>UNITS SUMMARY</h3>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Heat rates: W/cm²  |  Heat loads: J/cm²  |  TPS mass: g/cm²  |  Wall temp: K</div>
<div style="padding:2px 4px;margin:1px 0;background:#f8f8f8">Load factor: g's (n_max=2.5g typical)  |  Dynamic pressure: Pa</div>`;

// ================================================================
// Export CSV
// ================================================================
function exportCSV(tab) {
  const map = {stagnation:'s-output', trajectory:'tr-output', corridor:'co-output'};
  const el = document.getElementById(map[tab]);
  if (!el || !el.value.trim()) { alert('No data to export. Run a calculation first.'); return; }
  const blob = new Blob([el.value], {type:'text/plain'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = tab + '_output.txt';
  a.click();
}

// ================================================================
// Helpers
// ================================================================
function gf(id) { return +document.getElementById(id).value; }
function gv(id) { return document.getElementById(id).value; }
function sv(id,v) { document.getElementById(id).value = v; }
