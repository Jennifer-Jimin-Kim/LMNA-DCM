function initSpatial(){
"use strict";
var SP=window.__SPATIAL__;
if(!SP)return;
var gene='POSTN', mode='gene', showImg=true, alpha=0.85;
var imgs={}, loaded=0;

function css(v){return getComputedStyle(document.documentElement).getPropertyValue(v).trim();}
function hexToRgb(h){h=h.replace('#','');if(h.length===3)h=h[0]+h[0]+h[1]+h[1]+h[2]+h[2];
  return [parseInt(h.slice(0,2),16),parseInt(h.slice(2,4),16),parseInt(h.slice(4,6),16)];}
function mix(a,b,t){var A=hexToRgb(a),B=hexToRgb(b);
  return 'rgb('+Math.round(A[0]+(B[0]-A[0])*t)+','+Math.round(A[1]+(B[1]-A[1])*t)+','+Math.round(A[2]+(B[2]-A[2])*t)+')';}

/* build panels */
var wrap=document.getElementById('spatialPanels');
SP.samples.forEach(function(s,i){
  var d=document.createElement('div');
  d.className='panel';
  d.innerHTML='<div class="ptitle"><i class="dot" style="background:var('+
    (s.cohort==='Control'?'--coh1':'--coh2')+')"></i><span>'+s.name+'</span>'+
    '<span class="n">'+s.cohort+' · '+s.x.length.toLocaleString()+' spots</span></div>'+
    '<canvas class="umap" id="sp'+i+'" aria-label="Visium section '+s.name+'"></canvas>'+
    '<div class="tip" id="sptip'+i+'"></div>';
  wrap.appendChild(d);
  var im=new Image();
  im.onload=function(){loaded++;render();};
  im.onerror=function(){loaded++;render();};
  im.src=s.img;
  imgs[i]=im;
});

function fit(s,im,w,h){
  var iw=(im&&im.naturalWidth)||600, ih=(im&&im.naturalHeight)||600;
  var sc=Math.min(w/iw,h/ih);
  return {sc:sc,ox:(w-iw*sc)/2,oy:(h-ih*sc)/2,iw:iw,ih:ih};
}
var frames={};
function drawOne(i){
  var s=SP.samples[i], cv=document.getElementById('sp'+i);
  if(!cv)return;
  var w=cv.clientWidth,h=cv.clientHeight;
  if(!w||!h)return;
  var dpr=window.devicePixelRatio||1;
  cv.width=w*dpr;cv.height=h*dpr;
  var ctx=cv.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.fillStyle=css('--chart-surface')||'#fff';ctx.fillRect(0,0,w,h);
  var im=imgs[i], F=fit(s,im,w,h);
  frames[i]=F;
  if(im&&im.naturalWidth&&showImg){
    ctx.globalAlpha=1;
    ctx.drawImage(im,F.ox,F.oy,F.iw*F.sc,F.ih*F.sc);
  }
  var lo=css('--seq-lo'),hi=css('--seq-hi');
  var g=s.genes[gene];
  var r=Math.max(1.1,s.spot_r*F.sc*0.95);
  var order=[];
  for(var k=0;k<s.x.length;k++)order.push(k);
  if(mode==='gene'&&g)order.sort(function(a,b){return g.v[a]-g.v[b];});
  for(var q=0;q<order.length;q++){
    var j=order[q], lvl;
    if(mode==='gene'&&g)lvl=g.v[j]/255;
    else lvl=Math.min(1,Math.log1p(s.umi[j])/Math.log1p(20000));
    if(mode==='gene'&&lvl<0.03)continue;
    ctx.fillStyle=mix(lo,hi,Math.pow(lvl,0.7));
    ctx.globalAlpha=showImg?(0.25+lvl*0.75)*alpha:0.9;
    ctx.beginPath();
    ctx.arc(F.ox+s.x[j]*F.sc,F.oy+s.y[j]*F.sc,r,0,6.2832);
    ctx.fill();
  }
  ctx.globalAlpha=1;
}
function render(){
  SP.samples.forEach(function(_,i){drawOne(i);});
  var sg=document.getElementById('spScaleGene');
  if(sg)sg.textContent=mode==='gene'?gene:'total UMI';
}

/* hover */
SP.samples.forEach(function(s,i){
  var cv=document.getElementById('sp'+i), tip=document.getElementById('sptip'+i);
  if(!cv)return;
  cv.addEventListener('mousemove',function(ev){
    var F=frames[i]; if(!F){tip.style.opacity=0;return;}
    var rect=cv.getBoundingClientRect();
    var mx=ev.clientX-rect.left, my=ev.clientY-rect.top;
    var best=-1,bd=1e9;
    for(var k=0;k<s.x.length;k++){
      var dx=F.ox+s.x[k]*F.sc-mx, dy=F.oy+s.y[k]*F.sc-my, d=dx*dx+dy*dy;
      if(d<bd){bd=d;best=k;}
    }
    if(best>=0&&bd<400){
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
});

/* controls */
(function(){
  var sel=document.getElementById('spGene');
  SP.panel.forEach(function(g){
    var o=document.createElement('option');o.value=g;o.textContent=g;
    if(g===gene)o.selected=true;
    sel.appendChild(o);
  });
  sel.addEventListener('change',function(){gene=sel.value;mode='gene';syncSeg();render();});
  var bg=document.getElementById('spModeGene'), bu=document.getElementById('spModeUmi');
  function syncSeg(){
    bg.setAttribute('aria-pressed',mode==='gene');
    bu.setAttribute('aria-pressed',mode==='umi');
  }
  bg.addEventListener('click',function(){mode='gene';syncSeg();render();});
  bu.addEventListener('click',function(){mode='umi';syncSeg();render();});
  var ib=document.getElementById('spImg');
  ib.addEventListener('click',function(){
    showImg=!showImg;ib.setAttribute('aria-pressed',showImg);
    ib.textContent=showImg?'H&E on':'H&E off';render();
  });
  syncSeg();
})();

/* metrics table */
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

window.__renderSpatial=render;
var t;window.addEventListener('resize',function(){clearTimeout(t);t=setTimeout(render,150);});
render();
}
