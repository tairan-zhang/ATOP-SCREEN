export const views=['density','removal','scatter','components','signal'];
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const number=v=>v===null?'—':Number(v.toFixed(3)).toString();
export function renderNcsSvg(data,view,labels,options={}){
 if(!views.includes(view))throw Error('Unknown NCS view');
 const W=900,H=590,L=104,R=35,T=82,B=83,PW=W-L-R,PH=H-T-B;
 const pale=color=>'#'+color.slice(1).match(/../g).map(c=>Math.round(parseInt(c,16)*.25+255*.75).toString(16).padStart(2,'0')).join('');
 const blue=options.blue||'#00609D',red=options.red||'#C02635',font=Number(options.font)||14;
 const title=options.title||labels['ncs_'+view],text=(x,y,s,size=font,anchor='start',color='#1d1d1f')=>`<text x="${x}" y="${y}" font-size="${size}" text-anchor="${anchor}" fill="${color}">${esc(s)}</text>`;
 const line=(x1,y1,x2,y2,color='#666',dash='')=>`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="1.5" ${dash?'stroke-dasharray="'+dash+'"':''}/>`;
 let content=text(L,32,title,font+5)+text(L,57,`${labels.ncs_threshold}: ${number(data.params.threshold)} · ${labels.ncs_removed}: ${data.removed} / ${data.total} (${data.removedPercent.toFixed(1)}%)`,font-1);
 const X=(x,max=1)=>L+x/max*PW,Y=(y,max=1)=>T+PH-y/max*PH;
 function axes(xmax,ymax,xlabel,ylabel){content+=line(L,T,L,T+PH)+line(L,T+PH,W-R,T+PH);for(let i=0;i<=5;i++){const x=i*xmax/5,y=i*ymax/5;content+=line(X(x,xmax),T+PH,X(x,xmax),T+PH+5)+text(X(x,xmax),T+PH+24,number(x),font,'middle');content+=line(L-5,Y(y,ymax),L,Y(y,ymax))+text(L-12,Y(y,ymax)+5,number(y),font,'end');}content+=text(L+PW/2,H-22,options.xLabel||xlabel,font,'middle')+`<text transform="translate(16 ${T+PH/2}) rotate(-90)" font-size="${font}" text-anchor="middle">${esc(options.yLabel||ylabel)}</text>`;}
 function path(points,xmax,ymax,color=blue){let d='',valid=false;for(const[x,y]of points){if(y===null||!Number.isFinite(y)){valid=false;continue;}d+=`${valid?'L':'M'}${X(x,xmax).toFixed(2)},${Y(y,ymax).toFixed(2)} `;valid=true;}return `<path d="${d}" fill="none" stroke="${color}" stroke-width="2.4"/>`;}
 const clip=`<defs><clipPath id="ncs-plot-clip"><rect x="${L}" y="${T}" width="${PW}" height="${PH}"/></clipPath></defs>`;
 if(view==='density'){
  const ymax=Math.max(...data.density.map(r=>r[1]),1)*1.1;axes(1,ymax,'NCS',labels.ncs_density_axis);
  const d=`M${L},${T+PH} `+data.density.map(([x,y])=>`L${X(x)},${Y(y,ymax)}`).join(' ')+` L${X(1)},${T+PH} Z`;
  content+=`<path d="${d}" fill="${pale(blue)}"/>`+`<rect x="${L}" y="${T}" width="${PW*data.params.threshold}" height="${PH}" fill="${red}" opacity=".10"/>`+path(data.density,1,ymax)+line(X(data.params.threshold),T,X(data.params.threshold),T+PH,red);
  content+=text(W-R-12,T+25,`${labels.ncs_median}: ${number(data.median)} · P10: ${number(data.q10)}`,font,'end');
 } else if(view==='removal'){
  axes(1,100,labels.ncs_threshold,labels.ncs_removed_percent);
  content+=path(data.curve.map(r=>[r.threshold,r.removed/data.total*100]),1,100)+line(X(data.params.threshold),T,X(data.params.threshold),T+PH,red)+`<circle cx="${X(data.params.threshold)}" cy="${Y(data.removedPercent,100)}" r="6" fill="white" stroke="${red}" stroke-width="2"/>`;
 } else if(view==='scatter'){
  const max=Math.max(...data.scatter.map(r=>r[0]),.1)*1.03;axes(max,1,'|mean LFC|','NCS');
  content+=clip+`<g clip-path="url(#ncs-plot-clip)"><rect x="${L}" y="${Y(data.params.threshold)}" width="${PW}" height="${PH*data.params.threshold}" fill="${pale(red)}" opacity=".55"/>`;
  for(const[i,[x,y,guide,gene]]of data.scatter.entries())if(y!==null)content+=`<circle data-ncs-point="${i}" tabindex="0" role="button" aria-label="${esc(gene)} · ${esc(guide)} · NCS ${number(y)}" cx="${X(x,max).toFixed(2)}" cy="${Y(y).toFixed(2)}" r="1.8" fill="#96989D" opacity=".45"/>`;
  content+='</g>'+line(L,Y(data.params.threshold),W-R,Y(data.params.threshold),red);
  const trend=Array.from({length:20},(_,i)=>{const p=data.scatter.filter(r=>r[0]>=i*max/20&&r[0]<(i+1)*max/20&&r[1]!==null);return [(i+.5)*max/20,p.length>=5?p.reduce((s,r)=>s+r[1],0)/p.length:null];});content+=path(trend,max,1);
  content+=text(W-R-12,T+22,`${labels.ncs_displayed}: ${data.scatter.length} · ${labels.ncs_trend}`,font-2,'end');
 } else if(view==='components'){
  const edges=[0,.05,.1,.2,.4,1];
  const bins=['<0.05','0.05–0.10','0.10–0.20','0.20–0.40','≥0.40'];
  for(const[c,component]of ['dir','mag','stab'].entries()){
   const left=L+c*PW/3,width=PW/3-24;
   content+=text(left+width/2,T-12,labels['ncs_'+component],font,'middle')+line(left,T,left,T+PH)+line(left,T+PH,left+width,T+PH);
   for(let k=0;k<=4;k++)content+=text(left-6,Y(k/4)+4,number(k/4),font-3,'end');
   for(const box of data.boxes.filter(r=>r.component===component)){const x=left+(box.bin+.5)*width/5,bw=width/5*.55;content+=`<text transform="translate(${x},${T+PH+15}) rotate(-35)" font-size="${font-3}" text-anchor="end">${esc(bins[box.bin])}</text>`;if(!box.stats)continue;const[min,q1,med,q3,max]=box.stats,color=edges[box.bin+1]<=data.params.threshold?red:blue;content+=line(x,Y(min),x,Y(max),color)+`<rect x="${x-bw/2}" y="${Y(q3)}" width="${bw}" height="${Math.max(1,(q3-q1)*PH)}" fill="${pale(color)}" stroke="${color}"/>`+line(x-bw/2,Y(med),x+bw/2,Y(med),color);}
   const cutoff=data.params.threshold,bin=Math.min(4,edges.findIndex((v,i)=>i<5&&cutoff>=v&&cutoff<edges[i+1]));
   const fraction=cutoff>=1?5:Math.max(0,bin)+(cutoff-edges[Math.max(0,bin)])/(edges[Math.max(0,bin)+1]-edges[Math.max(0,bin)]);
   content+=line(left+fraction*width/5,T,left+fraction*width/5,T+PH,red,'4 5');
  }content+=text(L+PW/2,H-9,options.xLabel||labels.ncs_bins,font,'middle')+`<text transform="translate(16 ${T+PH/2}) rotate(-90)" font-size="${font}" text-anchor="middle">${esc(options.yLabel||labels.ncs_component_score)}</text>`;
 } else {
  const ymax=Math.max(120,...data.curve.flatMap(r=>[r.positive??0,r.negative??0]))*1.06;axes(1,ymax,labels.ncs_threshold,labels.ncs_signal_axis);
  content+=line(L,Y(100,ymax),W-R,Y(100,ymax),'#999','3 5')+path(data.curve.map(r=>[r.threshold,r.positive]),1,ymax,red)+path(data.curve.map(r=>[r.threshold,r.negative]),1,ymax,blue)+line(X(data.params.threshold),T,X(data.params.threshold),T+PH,'#999','4 4');
  content+=text(W-R-12,T+20,`${labels.ncs_positive} · Top ${data.params.topN}`,font-1,'end',red)+text(W-R-12,T+41,`${labels.ncs_negative} · Top ${data.params.topN}`,font-1,'end',blue);
 }
 return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(title)}"><rect x="-900" y="-590" width="2700" height="1770" fill="white"/><g font-family="Arial, sans-serif" fill="#1d1d1f">${content}</g></svg>`;
}
