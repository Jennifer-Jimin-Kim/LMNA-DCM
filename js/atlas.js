function initAtlas(){

"use strict";
var A=window.__ATLAS__, SP=window.__SPATIAL__;
var S=A.stats;

/* ---- decode blob ---- */
function b64ToBytes(b64){
  var bin=atob(b64),len=bin.length,out=new Uint8Array(len);
  for(var i=0;i<len;i++)out[i]=bin.charCodeAt(i);
  return out;
}
var bytes=b64ToBytes(A.blob), N=A.n, G=A.genes.length, off=0;
var xy=new Uint16Array(bytes.buffer,bytes.byteOffset+off,N*2); off+=N*4;
var ctArr=new Uint8Array(bytes.buffer,bytes.byteOffset+off,N); off+=N;
var lnArr=new Uint8Array(bytes.buffer,bytes.byteOffset+off,N); off+=N;
var dsArr=new Uint8Array(bytes.buffer,bytes.byteOffset+off,N); off+=N;
var expr=new Uint8Array(bytes.buffer,bytes.byteOffset+off,N*G);

/* ---- lineage colour map (fixed order by size) ---- */
var LN_ORDER=S.lineages.map(function(d){return d.name;});
var SLOTS=['--d1','--d2','--d3','--d4','--d5','--d6'];
var lnColor={},lnLabel={};
LN_ORDER.forEach(function(name,i){ lnColor[name]= i<6?SLOTS[i]:'--d0'; });
var LN_CATS=A.ln_cats;
function lineageVar(code){return lnColor[LN_CATS[code]]||'--d0';}

var DS=A.ds_cats;
var COHORT_LABEL={};
DS.forEach(function(d){COHORT_LABEL[d]=d;});

/* ---- helpers ---- */
function css(v){return getComputedStyle(document.documentElement).getPropertyValue(v).trim();}
function hexToRgb(h){h=h.replace('#','');if(h.length===3)h=h[0]+h[0]+h[1]+h[1]+h[2]+h[2];
  return [parseInt(h.slice(0,2),16),parseInt(h.slice(2,4),16),parseInt(h.slice(4,6),16)];}
function mix(a,b,t){var A1=hexToRgb(a),B=hexToRgb(b);
  return 'rgb('+Math.round(A1[0]+(B[0]-A1[0])*t)+','+Math.round(A1[1]+(B[1]-A1[1])*t)+','+Math.round(A1[2]+(B[2]-A1[2])*t)+')';}
function fmt(n){return n.toLocaleString();}

/* ---- hero stats ---- */
(function(){
  var el=document.getElementById('heroStats');
  [[fmt(S.n_cells),'Nuclei'],[String(S.n_lineages),'Lineages'],
   [String(S.n_types),'Cell states'],[String(S.sources.length),'Source cohorts']]
  .forEach(function(p){
    var d=document.createElement('div');d.className='hstat';
    d.innerHTML='<div class="v">'+p[0]+'</div><div class="k">'+p[1]+'</div>';
    el.appendChild(d);
  });
})();

/* ---- explorer state ---- */
var mode='expr', gene='LMNA', gIdx=A.genes.indexOf('LMNA'), highlight=null, frames={}, view='split', ptScale=0.85;

var UMAP_ASPECT=A.aspect||1;
function project(w,h){
  var pw=w-16, ph=h-16, s=Math.min(pw/UMAP_ASPECT,ph);
  var dw=s*UMAP_ASPECT, dh=s;
  return {ox:(w-dw)/2, oy:(h-dh)/2, dw:dw, dh:dh};
}
function drawUmap(canvas,dsCode,radius){
  radius=(radius||1.6)*ptScale;
  if(!canvas)return null;
  var w=canvas.clientWidth,h=canvas.clientHeight;
  if(!w||!h)return null;
  var dpr=window.devicePixelRatio||1;
  canvas.width=w*dpr;canvas.height=h*dpr;
  var ctx=canvas.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.fillStyle=css('--chart-surface')||'#fff';ctx.fillRect(0,0,w,h);
  var P=project(w,h);
  var lo=css('--seq-lo'),hi=css('--seq-hi'),dim=css('--line'),base=css('--rule-strong');
  var colCache={};
  function colOf(v){if(!colCache[v])colCache[v]=css(v);return colCache[v];}
  var pts=[],order=[];
  for(var i=0;i<N;i++){
    if(dsCode>=0&&dsArr[i]!==dsCode)continue;
    var x=P.ox+xy[i*2]/65535*P.dw, y=P.oy+(1-xy[i*2+1]/65535)*P.dh;
    pts.push({x:x/w,y:y/h,i:i});
    order.push(i);
  }
  function px(i){return P.ox+xy[i*2]/65535*P.dw;}
  function py(i){return P.oy+(1-xy[i*2+1]/65535)*P.dh;}
  if(mode==='expr'){
    /* grey base so the embedding stays readable */
    ctx.fillStyle=base;ctx.globalAlpha=0.55;
    for(var a=0;a<order.length;a++){
      ctx.beginPath();ctx.arc(px(order[a]),py(order[a]),radius*0.94,0,6.2832);ctx.fill();
    }
    var sorted=order.slice().sort(function(p,q){
      return (gIdx>=0?expr[gIdx*N+p]:0)-(gIdx>=0?expr[gIdx*N+q]:0);
    });
    for(var b=0;b<sorted.length;b++){
      var j=sorted[b];
      var lvl=gIdx>=0?expr[gIdx*N+j]/255:0;
      if(lvl<0.04)continue;
      var on=(highlight===null)||(LN_CATS[lnArr[j]]===highlight);
      if(!on)continue;
      ctx.fillStyle=mix(lo,hi,Math.pow(lvl,0.7));
      ctx.globalAlpha=0.45+lvl*0.55;
      ctx.beginPath();ctx.arc(px(j),py(j),radius*1.18,0,6.2832);ctx.fill();
    }
  }else{
    for(var c=0;c<order.length;c++){
      var k=order[c];
      var onk=(highlight===null)||(LN_CATS[lnArr[k]]===highlight);
      ctx.fillStyle=onk?colOf(lineageVar(lnArr[k])):dim;
      ctx.globalAlpha=onk?0.78:0.4;
      ctx.beginPath();ctx.arc(px(k),py(k),onk?radius*1.12:radius*0.82,0,6.2832);ctx.fill();
    }
  }
  ctx.globalAlpha=1;
  return {w:w,h:h,pts:pts};
}
function renderMain(){
  var split=document.getElementById('panelsSplit'), solo=document.getElementById('panelSingle');
  if(view==='split'){
    split.hidden=false;solo.hidden=true;
    frames.A=drawUmap(document.getElementById('umapA'),0,1.5);
    frames.B=drawUmap(document.getElementById('umapB'),1,1.5);
  }else{
    split.hidden=true;solo.hidden=false;
    var code=view==='all'?-1:parseInt(view,10);
    frames.S=drawUmap(document.getElementById('umapSolo'),code,2.0);
    var n=code<0?S.n_cells:S.datasets.filter(function(d){return d.name===DS[code];})[0].n;
    document.getElementById('soloTitle').textContent=code<0?'All cells':COHORT_LABEL[DS[code]];
    document.getElementById('soloN').textContent=fmt(n)+' nuclei';
    document.getElementById('soloDot').style.background=code<0?'var(--faint)':(code===0?'var(--coh1)':'var(--coh2)');
    var counts={},tot=0;
    for(var q=0;q<N;q++){
      if(code>=0&&dsArr[q]!==code)continue;
      var nm=LN_CATS[lnArr[q]];counts[nm]=(counts[nm]||0)+1;tot++;
    }
    var leg=document.getElementById('soloLeg');leg.innerHTML='';
    S.lineages.forEach(function(l){
      var c=counts[l.name]||0;
      var d=document.createElement('div');d.className='slrow';
      d.innerHTML='<i class="dot" style="background:var('+(lnColor[l.name]||'--d0')+')"></i>'+
        '<span class="nm">'+l.name+'</span><span class="pc">'+(tot?(c/tot*100).toFixed(1):'0.0')+'%</span>';
      leg.appendChild(d);
    });
  }
  var sg=document.getElementById('scaleGene');
  if(sg)sg.textContent=(mode==='expr'&&gIdx>=0)?gene:'';
  document.getElementById('scale').style.visibility=mode==='expr'?'visible':'hidden';
}

/* ---- overview UMAP with state labels ---- */
var CENTROIDS=null;
function centroids(){
  if(CENTROIDS)return CENTROIDS;
  var acc={};
  for(var i=0;i<N;i++){
    var name=A.ct_cats[ctArr[i]];
    if(!acc[name])acc[name]={x:0,y:0,n:0,ln:LN_CATS[lnArr[i]]};
    acc[name].x+=xy[i*2]/65535; acc[name].y+=1-xy[i*2+1]/65535; acc[name].n++;
  }
  CENTROIDS=Object.keys(acc).map(function(k){
    var a=acc[k];return {name:k,x:a.x/a.n,y:a.y/a.n,n:a.n,ln:a.ln};
  }).sort(function(a,b){return b.n-a.n;});
  return CENTROIDS;
}
function renderOverview(){
  var cv=document.getElementById('overviewUmap');
  if(!cv)return;
  var w=cv.clientWidth,h=cv.clientHeight;
  if(!w||!h)return;
  var dpr=window.devicePixelRatio||1;
  cv.width=w*dpr;cv.height=h*dpr;
  var ctx=cv.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.fillStyle=css('--chart-surface')||'#fff';ctx.fillRect(0,0,w,h);
  var P=project(w,h);
  var cache={};
  function colOf(v){if(!cache[v])cache[v]=css(v);return cache[v];}
  for(var i=0;i<N;i++){
    ctx.fillStyle=colOf(lineageVar(lnArr[i]));
    ctx.globalAlpha=0.68;
    ctx.beginPath();
    ctx.arc(P.ox+xy[i*2]/65535*P.dw,P.oy+(1-xy[i*2+1]/65535)*P.dh,1.9*ptScale,0,6.2832);
    ctx.fill();
  }
  ctx.globalAlpha=1;
  var cs=centroids(), ink=css('--ink'), plate=css('--chart-surface');
  var placed=[];
  ctx.font='600 12.5px '+css('--font');
  ctx.textAlign='center';ctx.textBaseline='middle';
  cs.forEach(function(c){
    if(c.n<25)return;
    var x=P.ox+c.x*P.dw, y=P.oy+c.y*P.dh;
    var wpx=ctx.measureText(c.name).width;
    for(var k=0;k<placed.length;k++){
      var q=placed[k];
      if(Math.abs(q.x-x)<(q.w+wpx)/2+8 && Math.abs(q.y-y)<19){ y=q.y+20; }
    }
    x=Math.max(P.ox+wpx/2+4,Math.min(P.ox+P.dw-wpx/2-4,x));
    y=Math.max(14,Math.min(h-10,y));
    ctx.lineWidth=3.6;ctx.strokeStyle=plate;
    ctx.strokeText(c.name,x,y);
    ctx.fillStyle=ink;
    ctx.fillText(c.name,x,y);
    placed.push({x:x,y:y,w:wpx});
  });
  var ov=document.getElementById('ovN');
  if(ov)ov.textContent=fmt(S.n_cells)+' nuclei · '+S.n_types+' states · '+S.n_lineages+' lineages';
  var ol=document.getElementById('olegend');
  if(ol&&!ol.childElementCount){
    S.lineages.forEach(function(l){
      var d=document.createElement('div');d.className='olrow';
      d.innerHTML='<i class="dot" style="background:var('+(lnColor[l.name]||'--d0')+')"></i>'+
        '<span class="nm">'+l.name+'</span><span class="pc">'+(l.n/S.n_cells*100).toFixed(1)+'%</span>';
      ol.appendChild(d);
    });
  }
}

/* ---- lineage cards ---- */
function renderThumbs(){
  S.lineages.forEach(function(l,idx){
    var cv=document.getElementById('thumb'+idx);
    if(!cv)return;
    var w=cv.clientWidth,h=cv.clientHeight;
    if(!w||!h)return;
    var dpr=window.devicePixelRatio||1;
    cv.width=w*dpr;cv.height=h*dpr;
    var ctx=cv.getContext('2d');
    ctx.setTransform(dpr,0,0,dpr,0,0);
    ctx.fillStyle=css('--chart-surface')||'#fff';ctx.fillRect(0,0,w,h);
    var on=css(lnColor[l.name]||'--d0'), dim=css('--line');
    var P=project(w,h);
    for(var i=0;i<N;i+=2){
      var hit=LN_CATS[lnArr[i]]===l.name;
      ctx.fillStyle=hit?on:dim;ctx.globalAlpha=hit?0.9:0.55;
      ctx.beginPath();
      ctx.arc(P.ox+xy[i*2]/65535*P.dw,P.oy+(1-xy[i*2+1]/65535)*P.dh,hit?1.4:0.9,0,6.2832);
      ctx.fill();
    }
    ctx.globalAlpha=1;
  });
}
(function(){
  var wrap=document.getElementById('compCards');
  S.lineages.forEach(function(l,idx){
    var b=document.createElement('button');
    b.className='ccard';b.type='button';b.setAttribute('aria-pressed','false');
    b.innerHTML='<canvas id="thumb'+idx+'" aria-label="'+l.name+'"></canvas>'+
      '<div class="meta"><div class="nm"><i class="dot" style="background:var('+(lnColor[l.name]||'--d0')+')"></i>'+
      l.name+'</div><div class="ct">'+fmt(l.n)+' cells · '+(l.n/S.n_cells*100).toFixed(1)+'%</div></div>';
    b.addEventListener('click',function(){
      var active=b.getAttribute('aria-pressed')==='true';
      wrap.querySelectorAll('.ccard').forEach(function(o){o.setAttribute('aria-pressed','false');});
      highlight=active?null:l.name;
      b.setAttribute('aria-pressed',active?'false':'true');
      renderMain();
      if(!active)document.getElementById('explore').scrollIntoView({behavior:'smooth'});
    });
    wrap.appendChild(b);
  });
  var lg=document.getElementById('legend');
  S.lineages.forEach(function(l){
    var s=document.createElement('span');
    s.innerHTML='<i class="dot" style="background:var('+(lnColor[l.name]||'--d0')+')"></i>'+l.name;
    lg.appendChild(s);
  });
  document.getElementById('cohA').textContent=COHORT_LABEL[DS[0]];
  document.getElementById('cohB').textContent=COHORT_LABEL[DS[1]];
  document.getElementById('cohAn').textContent=fmt(S.datasets.filter(function(d){return d.name===DS[0];})[0].n)+' nuclei';
  document.getElementById('cohBn').textContent=fmt(S.datasets.filter(function(d){return d.name===DS[1];})[0].n)+' nuclei';
})();

/* ---- gene search ---- */
function syncMode(){
  document.getElementById('modeExpr').setAttribute('aria-pressed',mode==='expr');
  document.getElementById('modeType').setAttribute('aria-pressed',mode==='type');
}
function setGene(g){
  gene=(g||'').trim().toUpperCase();
  gIdx=A.genes.indexOf(gene);
  document.getElementById('geneInput').value=gene;
  document.querySelectorAll('#suggest button').forEach(function(b){
    b.setAttribute('aria-pressed',b.textContent===gene);
  });
  mode='expr';syncMode();renderMain();
}
(function(){
  var sug=document.getElementById('suggest');
  ['LMNA','PARP1','POSTN','NPPA','TTN','PECAM1'].forEach(function(g){
    if(A.genes.indexOf(g)<0)return;
    var b=document.createElement('button');b.type='button';b.textContent=g;
    b.setAttribute('aria-pressed','false');
    b.addEventListener('click',function(){setGene(g);});
    sug.appendChild(b);
  });
  var t,input=document.getElementById('geneInput');
  input.value=gene;
  input.addEventListener('input',function(){
    clearTimeout(t);t=setTimeout(function(){setGene(input.value);},180);
  });
  document.getElementById('modeExpr').addEventListener('click',function(){mode='expr';syncMode();renderMain();});
  document.getElementById('modeType').addEventListener('click',function(){mode='type';syncMode();renderMain();});
})();

/* ---- hover ---- */
function wireHover(cid,tid,key){
  var cv=document.getElementById(cid),tip=document.getElementById(tid);
  if(!cv||!tip)return;
  cv.addEventListener('mousemove',function(ev){
    var f=frames[key];if(!f){tip.style.opacity=0;return;}
    var rect=cv.getBoundingClientRect();
    var mx=(ev.clientX-rect.left)/rect.width,my=(ev.clientY-rect.top)/rect.height;
    var best=null,bd=1e9;
    for(var k=0;k<f.pts.length;k++){
      var p=f.pts[k],dx=p.x-mx,dy=p.y-my,d=dx*dx+dy*dy;
      if(d<bd){bd=d;best=p;}
    }
    if(best&&bd<0.0008){
      var i=best.i, label=A.ct_cats[ctArr[i]];
      if(mode==='expr'&&gIdx>=0){
        label+='  ·  '+gene+' '+(expr[gIdx*N+i]/255*A.gmax[gIdx]).toFixed(2);
      }
      tip.textContent=label;
      tip.style.left=Math.min(rect.width-tip.offsetWidth-8,ev.clientX-rect.left+12)+'px';
      tip.style.top=(ev.clientY-rect.top-26)+'px';
      tip.style.opacity=1;
    }else tip.style.opacity=0;
  });
  cv.addEventListener('mouseleave',function(){tip.style.opacity=0;});
}
wireHover('umapA','tipA','A');
wireHover('umapB','tipB','B');
wireHover('umapSolo','tipSolo','S');

/* ---- point size ---- */
(function(){
  var a=document.getElementById('ptSizeA'), b=document.getElementById('ptSizeB');
  var oa=document.getElementById('ptOutA'), ob=document.getElementById('ptOutB');
  function apply(v,src){
    ptScale=parseFloat(v);
    if(a&&src!==a)a.value=v; if(b&&src!==b)b.value=v;
    var txt=ptScale.toFixed(1);
    if(oa)oa.textContent=txt; if(ob)ob.textContent=txt;
    renderAll();
  }
  [a,b].forEach(function(el){
    if(!el)return;
    el.addEventListener('input',function(){apply(el.value,el);});
  });
})();

/* ---- dataset view selector ---- */
(function(){
  var seg=document.getElementById('viewSeg');
  document.getElementById('vbA').textContent=COHORT_LABEL[DS[0]];
  document.getElementById('vbB').textContent=COHORT_LABEL[DS[1]];
  seg.querySelectorAll('button').forEach(function(b){
    b.addEventListener('click',function(){
      view=b.getAttribute('data-view');
      seg.querySelectorAll('button').forEach(function(o){
        o.setAttribute('aria-pressed',o===b?'true':'false');
      });
      renderMain();
    });
  });
})();

/* ---- composition ---- */
(function(){
  var rows=document.getElementById('compRows'),tb=document.getElementById('compTbody');
  var comp=S.composition_lineage, a=DS[0], b=DS[1], max=0;
  comp.forEach(function(r){max=Math.max(max,r[a],r[b]);});
  document.getElementById('thA').textContent=COHORT_LABEL[a]+' %';
  document.getElementById('thB').textContent=COHORT_LABEL[b]+' %';
  comp.forEach(function(r){
    var d=document.createElement('div');d.className='row';
    d.innerHTML='<div class="lb"><i class="dot" style="background:var('+(lnColor[r.name]||'--d0')+')"></i>'+r.name+'</div>'+
      '<div class="bars">'+
      '<div class="barline"><span class="cap">'+COHORT_LABEL[a]+'</span><span class="bartrack">'+
      '<span class="bar" style="width:'+(r[a]/max*100)+'%;background:var(--coh1)"></span></span>'+
      '<span class="val">'+(r[a]*100).toFixed(1)+'%</span></div>'+
      '<div class="barline"><span class="cap">'+COHORT_LABEL[b]+'</span><span class="bartrack">'+
      '<span class="bar" style="width:'+(r[b]/max*100)+'%;background:var(--coh2)"></span></span>'+
      '<span class="val">'+(r[b]*100).toFixed(1)+'%</span></div></div>';
    rows.appendChild(d);
    var ratio=(r[a]>0&&r[b]>0)?Math.log2(r[b]/r[a]):null;
    var tr=document.createElement('tr');
    tr.innerHTML='<td><span class="cond"><i class="dot" style="background:var('+(lnColor[r.name]||'--d0')+')"></i>'+r.name+'</span></td>'+
      '<td class="num">'+(r[a]*100).toFixed(1)+'</td><td class="num">'+(r[b]*100).toFixed(1)+'</td>'+
      '<td class="num">'+(ratio===null?'—':(ratio>0?'+':'')+ratio.toFixed(2))+'</td>';
    tb.appendChild(tr);
  });
  var vc=document.getElementById('viewChart'),vt=document.getElementById('viewTable');
  var chart=document.getElementById('compChart'),table=document.getElementById('compTable');
  function set(isChart){chart.hidden=!isChart;table.hidden=isChart;
    vc.setAttribute('aria-pressed',isChart);vt.setAttribute('aria-pressed',!isChart);}
  vc.addEventListener('click',function(){set(true);});
  vt.addEventListener('click',function(){set(false);});
})();

/* ---- cell states table ---- */
(function(){
  var body=document.getElementById('stateBody'),empty=document.getElementById('stateEmpty');
  var input=document.getElementById('stateSearch');
  var lineageOf={};
  // derive lineage per state from the subsample
  for(var i=0;i<N;i++){ lineageOf[A.ct_cats[ctArr[i]]]=LN_CATS[lnArr[i]]; }
  var rows=S.cell_types.map(function(t){
    return {name:t.name,n:t.n,lineage:lineageOf[t.name]||'—'};
  });
  var sortKey='n',sortDir=-1,maxN=Math.max.apply(null,rows.map(function(r){return r.n;}));
  function render(){
    var q=(input.value||'').trim().toUpperCase();
    var r=rows.filter(function(x){return !q||x.name.toUpperCase().indexOf(q)>-1||x.lineage.toUpperCase().indexOf(q)>-1;})
      .sort(function(x,y){
        var a=x[sortKey],b=y[sortKey];
        if(typeof a==='number')return (a-b)*sortDir;
        return String(a).localeCompare(String(b))*sortDir;
      });
    body.innerHTML='';
    r.forEach(function(x){
      var tr=document.createElement('tr');
      tr.innerHTML='<td class="id">'+x.name+'</td>'+
        '<td><span class="cond"><i class="dot" style="background:var('+(lnColor[x.lineage]||'--d0')+')"></i>'+x.lineage+'</span></td>'+
        '<td class="num">'+fmt(x.n)+'</td>'+
        '<td><span class="minibar" style="width:'+Math.max(3,x.n/maxN*120)+'px;background:var('+(lnColor[x.lineage]||'--d0')+')"></span></td>';
      body.appendChild(tr);
    });
    empty.hidden=r.length>0;
  }
  input.addEventListener('input',render);
  document.querySelectorAll('#states th.sortable').forEach(function(th){
    th.addEventListener('click',function(){
      var k=th.getAttribute('data-sort');
      if(k===sortKey)sortDir=-sortDir;else{sortKey=k;sortDir=k==='n'?-1:1;}
      document.querySelectorAll('#states th.sortable .ar').forEach(function(a){a.textContent='';});
      th.querySelector('.ar').textContent=sortDir>0?'↑':'↓';
      render();
    });
  });
  render();
})();

/* ---- donors + facets ---- */
(function(){
  var body=document.getElementById('sampleBody'),empty=document.getElementById('sampleEmpty');
  var count=document.getElementById('resultCount');
  var rows=S.donors.slice();
  var sortKey='n',sortDir=-1;
  var groups=document.querySelectorAll('#samples .fgroup');
  groups.forEach(function(g){
    var key=g.getAttribute('data-facet');
    var vals={};
    rows.forEach(function(r){vals[r[key]]=(vals[r[key]]||0)+1;});
    var box=g.querySelector('[data-opts]');
    Object.keys(vals).sort().forEach(function(v){
      var l=document.createElement('label');l.className='fopt';
      l.innerHTML='<input type="checkbox" value="'+v+'"> '+v+
        ' <span class="cnt">'+vals[v]+'</span>';
      box.appendChild(l);
    });
  });
  function render(){
    var f={};
    groups.forEach(function(g){
      var k=g.getAttribute('data-facet'),sel=[];
      g.querySelectorAll('input:checked').forEach(function(i){sel.push(i.value);});
      if(sel.length)f[k]=sel;
    });
    var r=rows.filter(function(d){
      return Object.keys(f).every(function(k){return f[k].indexOf(d[k])>-1;});
    }).sort(function(x,y){
      var a=x[sortKey],b=y[sortKey];
      if(typeof a==='number')return (a-b)*sortDir;
      return String(a).localeCompare(String(b))*sortDir;
    });
    body.innerHTML='';
    r.forEach(function(d){
      var col=d.group===DS[0]?'--coh1':'--coh2';
      var tr=document.createElement('tr');
      tr.innerHTML='<td class="id">'+d.name+'</td>'+
        '<td><span class="cond"><i class="dot" style="background:var('+col+')"></i>'+d.group+'</span></td>'+
        '<td class="m">'+d.species+'</td>'+
        '<td class="m">'+d.dataset+'</td>'+
        '<td class="num">'+fmt(d.n)+'</td>';
      body.appendChild(tr);
    });
    empty.hidden=r.length>0;
    count.textContent=r.length+' of '+rows.length+' samples';
  }
  document.querySelectorAll('#samples .fgroup input').forEach(function(i){i.addEventListener('change',render);});
  document.getElementById('resetFacets').addEventListener('click',function(){
    document.querySelectorAll('#samples .fgroup input').forEach(function(i){i.checked=false;});render();
  });
  document.querySelectorAll('#samples th.sortable').forEach(function(th){
    th.addEventListener('click',function(){
      var k=th.getAttribute('data-sort');
      if(k===sortKey)sortDir=-sortDir;else{sortKey=k;sortDir=k==='n'?-1:1;}
      document.querySelectorAll('#samples th.sortable .ar').forEach(function(a){a.textContent='';});
      th.querySelector('.ar').textContent=sortDir>0?'↑':'↓';
      render();
    });
  });
  render();
})();

/* ---- markers ---- */
(function(){
  var body=document.getElementById('markerBody');
  var mbt=S.mean_by_type, types=Object.keys(mbt);
  var lineageOf={};
  for(var i=0;i<N;i++){ lineageOf[A.ct_cats[ctArr[i]]]=LN_CATS[lnArr[i]]; }
  var rows=S.genes.map(function(g,gi){
    var best=null,bv=-1;
    types.forEach(function(t){ if(mbt[t][gi]>bv){bv=mbt[t][gi];best=t;} });
    return {gene:g,type:best,val:bv};
  }).sort(function(a,b){return b.val-a.val;});
  var maxV=rows[0].val||1;
  rows.forEach(function(r){
    var lin=lineageOf[r.type]||'—';
    var tr=document.createElement('tr');
    tr.style.cursor='pointer';
    tr.innerHTML='<td class="id">'+r.gene+'</td>'+
      '<td><span class="cond"><i class="dot" style="background:var('+(lnColor[lin]||'--d0')+')"></i>'+r.type+'</span></td>'+
      '<td class="num">'+r.val.toFixed(2)+'</td>'+
      '<td><span class="minibar" style="width:'+Math.max(3,r.val/maxV*120)+'px;background:var('+(lnColor[lin]||'--d0')+')"></span></td>';
    tr.addEventListener('click',function(){
      setGene(r.gene);
      document.getElementById('explore').scrollIntoView({behavior:'smooth'});
    });
    body.appendChild(tr);
  });
})();

/* ---- about meta ---- */
(function(){
  var dl=document.getElementById('aboutMeta');
  [['Nuclei',fmt(S.n_cells)],['Genes',fmt(S.n_genes)],['Samples',String(S.n_donors)],
   ['Cell states',String(S.n_types)],['Lineages',String(S.n_lineages)],
   ['Species',S.species.map(function(c){return c.name;}).join(' · ')],
   ['Source cohorts',S.sources.map(function(c){return c.name;}).join(' · ')],
   ['Groups',S.datasets.map(function(d){return d.name+' ('+fmt(d.n)+')';}).join(' · ')],
   ['Spatial','10x Visium · 4 sections'],
   ['Integration','Seurat v5 · Harmony'],['Licence','CC BY 4.0']]
  .forEach(function(p){
    dl.insertAdjacentHTML('beforeend','<dt>'+p[0]+'</dt><dd>'+p[1]+'</dd>');
  });
})();

/* ---- home link ---- */
(function(){
  var a=document.getElementById('homeLink');
  if(!a)return;
  a.addEventListener('click',function(ev){
    ev.preventDefault();
    window.scrollTo({top:0,behavior:'smooth'});
    if(history.replaceState)history.replaceState(null,'',location.pathname+location.search);
  });
})();

/* ---- theme ---- */
(function(){
  var root=document.documentElement, btn=document.getElementById('themeBtn');
  function stored(){try{return localStorage.getItem('lmna-theme');}catch(e){return null;}}
  function save(v){try{localStorage.setItem('lmna-theme',v);}catch(e){}}
  var t=stored();
  if(t==='dark'||t==='light')root.setAttribute('data-theme',t);
  function current(){
    var a=root.getAttribute('data-theme');
    if(a)return a;
    return (window.matchMedia&&window.matchMedia('(prefers-color-scheme: dark)').matches)?'dark':'light';
  }
  btn.addEventListener('click',function(){
    var next=current()==='dark'?'light':'dark';
    root.setAttribute('data-theme',next);save(next);renderAll();
  });
})();

/* ---- lifecycle ---- */
function renderAll(){renderMain();renderThumbs();renderOverview();}
syncMode();
renderAll();
var rt;window.addEventListener('resize',function(){clearTimeout(rt);rt=setTimeout(renderAll,150);});
if(window.matchMedia){var mq=window.matchMedia('(prefers-color-scheme: dark)');
  if(mq.addEventListener)mq.addEventListener('change',renderAll);}
if(document.fonts&&document.fonts.ready)document.fonts.ready.then(renderAll);

window.__renderAll=renderAll;
}
