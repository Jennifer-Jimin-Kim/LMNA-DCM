function initSpatial(){
"use strict";
var SP=window.__SPATIAL__;
if(!SP)return;

/* genes with the clearest zonation (Moran's I in LMNA sections, gained vs control) */
var FEATURED=['NPPB','NPPA','POSTN','MYH7','DCN'];
var FEAT_NOTE={
  NPPB:'Natriuretic peptide B — stress-zone marker, sharply regionalised in LMNA-DCM',
  NPPA:'Natriuretic peptide A — co-zonated with NPPB, near-absent in control',
  POSTN:'Periostin — activated-fibroblast zones tracking the fibrotic front',
  MYH7:'β-myosin heavy chain — strongest regional structure of the panel',
  DCN:'Decorin — interstitial ECM, broader zonation than POSTN'
};

var gene='NPPB', mode='gene', showImg=true, view='single', cur=0;
cur=Math.max(0,SP.samples.map(function(s){return s.name;}).indexOf('LMNA-1'));
var imgs={}, Z={};

function css(v){return getComputedStyle(document.documentElement).getPropertyValue(v).trim();}
function hexToRgb(h){h=h.replace('#','');if(h.length===3)h=h[0]+h[0]+h[1]+h[1]+h[2]+h[2];
  return [parseInt(h.slice(0,2),16),parseInt(h.slice(2,4),16),parseInt(h.slice(4,6),16)];}
var SEQ=['#440154','#472D7B','#3B528B','#2C728E','#21918C','#28AE80','#5EC962','#ADDC30','#FDE725'];
var SEQRGB=SEQ.map(hexToRgb), RAMP=new Array(256);
(function(){for(var i=0;i<256;i++){
  var t=i/255*(SEQRGB.length-1), k=Math.min(SEQRGB.length-2,Math.floor(t)), f=t-k;
  var a=SEQRGB[k],b=SEQRGB[k+1];
  RAMP[i]='rgb('+Math.round(a[0]+(b[0]-a[0])*f)+','+Math.round(a[1]+(b[1]-a[1])*f)+','+Math.round(a[2]+(b[2]-a[2])*f)+')';
}})();
function seq(t){return RAMP[Math.max(0,Math.min(255,Math.round(t*255)))];}
function z(id){ if(!Z[id])Z[id]={k:1,tx:0,ty:0}; return Z[id]; }

/* ---------- panels ---------- */
function panelHTML(i,idPrefix){
  var s=SP.samples[i];
  return '<div class="ptitle"><i class="dot" style="background:var('+
    (s.cohort==='Control'?'--coh1':'--coh2')+')"></i><span>'+s.name+'</span>'+
    '<span class="n">'+s.cohort+' · '+s.x.length.toLocaleString()+' spots</span></div>'+
    '<div class="cwrap"><canvas class="umap" id="'+idPrefix+i+'" aria-label="Visium section '+s.name+'"></canvas>'+
    '<div class="tip" id="tip'+idPrefix+i+'"></div>'+
    '<div class="zoomctl"><button type="button" data-zin="'+idPrefix+i+'" title="Zoom in">+</button>'+
    '<button type="button" data-zout="'+idPrefix+i+'" title="Zoom out">−</button>'+
    '<button type="button" data-zreset="'+idPrefix+i+'" title="Reset view">⟳</button></div></div>';
}

var single=document.getElementById('spSingle'), grid=document.getElementById('spGrid');

function buildSingle(){
  single.innerHTML='<div class="panel wide">'+panelHTML(cur,'spS')+'</div>';
  wire('spS'+cur,cur);
}
function buildGrid(){
  grid.innerHTML='';
  SP.samples.forEach(function(s,i){
    var d=document.createElement('div');
    d.className='panel';
    d.innerHTML=panelHTML(i,'spG');
    grid.appendChild(d);
    wire('spG'+i,i);
  });
}

/* ---------- drawing ---------- */
function fit(im,w,h){
  var iw=(im&&im.naturalWidth)||600, ih=(im&&im.naturalHeight)||600;
  var sc=Math.min(w/iw,h/ih);
  return {sc:sc,ox:(w-iw*sc)/2,oy:(h-ih*sc)/2,iw:iw,ih:ih};
}
var frames={};
function draw(cid,i){
  var s=SP.samples[i], cv=document.getElementById(cid);
  if(!cv)return;
  var w=cv.clientWidth,h=cv.clientHeight;
  if(!w||!h)return;
  var dpr=window.devicePixelRatio||1;
  cv.width=w*dpr;cv.height=h*dpr;
  var ctx=cv.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.fillStyle=css('--chart-surface')||'#fff';ctx.fillRect(0,0,w,h);
  var im=imgs[i], F=fit(im,w,h), t=z(cid);
  var sc=F.sc*t.k, ox=F.ox*t.k+t.tx, oy=F.oy*t.k+t.ty;
  var Fv={sc:sc,ox:ox,oy:oy};
  frames[cid]=Fv;
  ctx.save();ctx.beginPath();ctx.rect(0,0,w,h);ctx.clip();
  if(im&&im.naturalWidth&&showImg){
    ctx.globalAlpha=1;
    ctx.drawImage(im,ox,oy,F.iw*sc,F.ih*sc);
  }
  var g=s.genes[gene];
  var r=Math.max(1.1,s.spot_r*sc*0.95);
  var order=[];
  for(var k=0;k<s.x.length;k++)order.push(k);
  if(mode==='gene'&&g)order.sort(function(a,b){return g.v[a]-g.v[b];});
  for(var q=0;q<order.length;q++){
    var j=order[q], lvl;
    if(mode==='gene'&&g)lvl=g.v[j]/255;
    else lvl=Math.min(1,Math.log1p(s.umi[j])/Math.log1p(20000));
    ctx.fillStyle=seq(Math.pow(lvl,0.65));
    ctx.globalAlpha=showImg?0.85:1;
    ctx.beginPath();
    ctx.arc(ox+s.x[j]*sc,oy+s.y[j]*sc,r,0,6.2832);
    ctx.fill();
  }
  ctx.globalAlpha=1;ctx.restore();
}
function render(){
  if(view==='single')draw('spS'+cur,cur);
  else SP.samples.forEach(function(_,i){draw('spG'+i,i);});
  var sg=document.getElementById('spScaleGene');
  if(sg)sg.textContent=mode==='gene'?gene:'total UMI';
  var nt=document.getElementById('spNote');
  if(nt)nt.textContent=(mode==='gene'&&FEAT_NOTE[gene])?FEAT_NOTE[gene]:'';
  var lbl=document.getElementById('spWhich');
  if(lbl)lbl.textContent=SP.samples[cur].name+' · '+(cur+1)+' / '+SP.samples.length;
  document.querySelectorAll('#spDots button').forEach(function(b,i){
    b.setAttribute('aria-pressed',i===cur);
  });
}

/* ---------- interaction ---------- */
function wire(cid,i){
  var cv=document.getElementById(cid), tip=document.getElementById('tip'+cid), s=SP.samples[i];
  if(!cv)return;
  var drag=null;
  cv.addEventListener('wheel',function(ev){
    ev.preventDefault();
    var t=z(cid), rect=cv.getBoundingClientRect();
    var mx=ev.clientX-rect.left, my=ev.clientY-rect.top;
    var f=Math.exp(-ev.deltaY*0.0015), nk=Math.max(1,Math.min(14,t.k*f));
    f=nk/t.k;
    t.tx=mx-(mx-t.tx)*f; t.ty=my-(my-t.ty)*f; t.k=nk;
    if(nk===1){t.tx=0;t.ty=0;}
    render();
  },{passive:false});
  cv.addEventListener('mousedown',function(ev){
    drag={x:ev.clientX,y:ev.clientY,t:z(cid),tx:z(cid).tx,ty:z(cid).ty};
    cv.style.cursor='grabbing';
  });
  window.addEventListener('mouseup',function(){if(drag){drag=null;cv.style.cursor='';}});
  cv.addEventListener('mousemove',function(ev){
    if(drag){
      drag.t.tx=drag.tx+(ev.clientX-drag.x);
      drag.t.ty=drag.ty+(ev.clientY-drag.y);
      render();return;
    }
    var F=frames[cid]; if(!F){tip.style.opacity=0;return;}
    var rect=cv.getBoundingClientRect();
    var mx=ev.clientX-rect.left, my=ev.clientY-rect.top;
    var best=-1,bd=1e9;
    for(var k=0;k<s.x.length;k++){
      var dx=F.ox+s.x[k]*F.sc-mx, dy=F.oy+s.y[k]*F.sc-my, d=dx*dx+dy*dy;
      if(d<bd){bd=d;best=k;}
    }
    var lim=Math.max(400,Math.pow(s.spot_r*F.sc*2.2,2));
    if(best>=0&&bd<lim){
      var g=s.genes[gene];
      tip.textContent=(mode==='gene'&&g)
        ? gene+' '+(g.v[best]/255*g.max).toFixed(2)
        : s.umi[best].toLocaleString()+' UMI';
      tip.style.left=Math.min(rect.width-tip.offsetWidth-8,mx+12)+'px';
      tip.style.top=(my-26)+'px';
      tip.style.opacity=1;
    }else tip.style.opacity=0;
  });
  cv.addEventListener('mouseleave',function(){tip.style.opacity=0;});
  cv.addEventListener('dblclick',function(){Z[cid]={k:1,tx:0,ty:0};render();});
}
document.addEventListener('click',function(ev){
  var b=ev.target.closest&&ev.target.closest('[data-zin],[data-zout],[data-zreset]');
  if(!b)return;
  var id=b.getAttribute('data-zin')||b.getAttribute('data-zout')||b.getAttribute('data-zreset');
  var cv=document.getElementById(id); if(!cv)return;
  var t=z(id), w=cv.clientWidth/2, h=cv.clientHeight/2;
  if(b.hasAttribute('data-zreset')){Z[id]={k:1,tx:0,ty:0};}
  else{
    var f=b.hasAttribute('data-zin')?1.45:1/1.45;
    var nk=Math.max(1,Math.min(14,t.k*f)); f=nk/t.k;
    t.tx=w-(w-t.tx)*f; t.ty=h-(h-t.ty)*f; t.k=nk;
    if(nk===1){t.tx=0;t.ty=0;}
  }
  render();
});

/* ---------- controls ---------- */
function go(n){
  cur=(n+SP.samples.length)%SP.samples.length;
  if(view==='single'){buildSingle();render();}
}
(function(){
  var chips=document.getElementById('spChips');
  FEATURED.forEach(function(g){
    var b=document.createElement('button');b.type='button';b.textContent=g;
    b.setAttribute('aria-pressed',g===gene);
    b.addEventListener('click',function(){
      gene=g;mode='gene';
      chips.querySelectorAll('button').forEach(function(o){o.setAttribute('aria-pressed',o===b);});
      document.getElementById('spGene').value='';
      syncSeg();render();
    });
    chips.appendChild(b);
  });
  var sel=document.getElementById('spGene');
  var o0=document.createElement('option');o0.value='';o0.textContent='All panel genes…';
  sel.appendChild(o0);
  SP.panel.forEach(function(g){
    var o=document.createElement('option');o.value=g;o.textContent=g;sel.appendChild(o);
  });
  sel.addEventListener('change',function(){
    if(!sel.value)return;
    gene=sel.value;mode='gene';
    chips.querySelectorAll('button').forEach(function(o){o.setAttribute('aria-pressed',o.textContent===gene);});
    syncSeg();render();
  });

  var bg=document.getElementById('spModeGene'), bu=document.getElementById('spModeUmi');
  function syncSeg(){
    bg.setAttribute('aria-pressed',mode==='gene');
    bu.setAttribute('aria-pressed',mode==='umi');
  }
  window.__spSync=syncSeg;
  bg.addEventListener('click',function(){mode='gene';syncSeg();render();});
  bu.addEventListener('click',function(){mode='umi';syncSeg();render();});

  var ib=document.getElementById('spImg');
  ib.addEventListener('click',function(){
    showImg=!showImg;ib.setAttribute('aria-pressed',showImg);
    ib.textContent=showImg?'H&E on':'H&E off';render();
  });

  var vs=document.getElementById('spViewSeg');
  vs.querySelectorAll('button').forEach(function(b){
    b.addEventListener('click',function(){
      view=b.getAttribute('data-spview');
      vs.querySelectorAll('button').forEach(function(o){o.setAttribute('aria-pressed',o===b);});
      single.hidden=view!=='single'; grid.hidden=view!=='grid';
      document.getElementById('spNav').hidden=view!=='single';
      if(view==='single')buildSingle(); else buildGrid();
      render();
    });
  });

  document.getElementById('spPrev').addEventListener('click',function(){go(cur-1);});
  document.getElementById('spNext').addEventListener('click',function(){go(cur+1);});
  var dots=document.getElementById('spDots');
  SP.samples.forEach(function(s,i){
    var b=document.createElement('button');b.type='button';b.textContent=s.name;
    b.setAttribute('aria-pressed',i===cur);
    b.addEventListener('click',function(){go(i);});
    dots.appendChild(b);
  });
  syncSeg();
})();

/* ---------- metrics table ---------- */
(function(){
  var b=document.getElementById('visBody');
  if(!b)return;
  SP.metrics.forEach(function(m){
    var tr=document.createElement('tr');
    tr.innerHTML='<td class="id">'+m.name+'</td>'+
      '<td><span class="cond"><i class="dot" style="background:var('+
        (m.cohort==='Control'?'--coh1':'--coh2')+')"></i>'+m.cohort+'</span></td>'+
      '<td class="m">'+m.slide+'</td>'+
      '<td class="num">'+m.spots.toLocaleString()+'</td>'+
      '<td class="num">'+Math.round(m.median_genes).toLocaleString()+'</td>'+
      '<td class="num">'+Math.round(m.median_umi).toLocaleString()+'</td>'+
      '<td class="num">'+m.genes_detected.toLocaleString()+'</td>';
    b.appendChild(tr);
  });
})();

/* ---------- boot ---------- */
SP.samples.forEach(function(s,i){
  var im=new Image();
  im.onload=function(){render();};
  im.onerror=function(){render();};
  im.src=s.img;
  imgs[i]=im;
});
buildSingle(); grid.hidden=true;
window.__renderSpatial=render;
var t;window.addEventListener('resize',function(){clearTimeout(t);t=setTimeout(render,150);});
render();
}
