export const quantile=(sorted,p)=>{if(!sorted.length)return null;const i=(sorted.length-1)*p,a=Math.floor(i);return sorted[a]+(sorted[Math.ceil(i)]-sorted[a])*(i-a);};
const median=v=>quantile([...v].sort((a,b)=>a-b),.5);
export function scoreGenes(rows,threshold,minGuides){const groups=new Map();for(const r of rows){if(r.ncs!==null&&r.ncs<threshold)continue;if(!groups.has(r.gene))groups.set(r.gene,[]);groups.get(r.gene).push(r.mean);}const scores=[];for(const[gene,v]of groups){if(v.length<minGuides)continue;v.sort((a,b)=>a-b);const k=Math.ceil(2*v.length/3);scores.push({gene,n:v.length,positive:v.slice(-k).reduce((a,b)=>a+b,0)/k,negative:v.slice(0,k).reduce((a,b)=>a+b,0)/k});}return scores;}
const signal=(scores,n)=>({positive:median(scores.map(r=>r.positive).sort((a,b)=>b-a).slice(0,n)),negative:median(scores.map(r=>r.negative).sort((a,b)=>a-b).slice(0,n)),genes:scores.length});
const lowerBound=(a,v)=>{let l=0,r=a.length;while(l<r){const m=(l+r)>>1;if(a[m]<v)l=m+1;else r=m;}return l;};
export function prepareModel(rows, minGuides, params) {
 const values=rows.map(r=>r.ncs).filter(x=>x!==null).sort((a,b)=>a-b),hist=Array(256).fill(0);for(const x of values)hist[Math.min(255,Math.floor(x*256))]++;
 // Gaussian KDE on a fixed grid; the count curve always uses exact observations.
 const bw=.02,density=hist.map((_,i)=>{const x=i/255;let sum=0;for(let j=0;j<256;j++){const z=(x-(j+.5)/256)/bw;if(Math.abs(z)<5)sum+=hist[j]*Math.exp(-.5*z*z);}return [x,sum/(Math.max(1,values.length)*bw*Math.sqrt(2*Math.PI))];});
 const bins=[0,.05,.1,.2,.4,1.000001],boxes=[];
 for(let b=0;b<5;b++){const bin=rows.filter(r=>r.ncs!==null&&r.ncs>=bins[b]&&r.ncs<bins[b+1]);for(const component of ['dir','mag','stab']){const v=bin.map(r=>r[component]).sort((a,b)=>a-b);boxes.push({bin:b,component,n:v.length,stats:v.length?[v[0],quantile(v,.25),quantile(v,.5),quantile(v,.75),v.at(-1)]:null});}}
 const base=signal(scoreGenes(rows,0,minGuides),params.topN);const curve=Array.from({length:41},(_,i)=>{const threshold=i/40,s=signal(scoreGenes(rows,threshold,minGuides),params.topN);return {threshold,removed:lowerBound(values,threshold),genes:s.genes,positive:base.positive===null||Math.abs(base.positive)<1e-12||s.positive===null?null:Math.abs(s.positive)/Math.abs(base.positive)*100,negative:base.negative===null||Math.abs(base.negative)<1e-12||s.negative===null?null:Math.abs(s.negative)/Math.abs(base.negative)*100};});
 const step=Math.max(1,Math.ceil(rows.length/4000)),scatter=rows.filter((_,i)=>i%step===0).map(r=>[Math.abs(r.mean),r.ncs,r.guide,r.gene,r.mean,r.dir,r.mag,r.stab]);
 return {rows,values,density,boxes,curve,scatter,total:rows.length,unknown:rows.length-values.length,median:quantile(values,.5),q10:quantile(values,.1),baseGenes:base.genes,base,minGuides:minGuides};
}
export function snapshot(model,params){const threshold=params.threshold,removed=lowerBound(model.values,threshold),genes=scoreGenes(model.rows,threshold,model.minGuides).length;const {rows,values,...data}=model;return {...data,params,removed,retained:model.total-removed,genes,removedPercent:100*removed/model.total};}
