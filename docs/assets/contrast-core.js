/* Synthetic teaching data only. Does not run limma or DESeq2. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.ContrastCore=api;})(typeof globalThis!=='undefined'?globalThis:this,function(){
'use strict';
function bh(values){if(values.some(p=>!Number.isFinite(p)||p<0||p>1))throw Error('P-values must be finite and between zero and one.');const order=values.map((_,i)=>i).sort((a,b)=>values[a]-values[b]),q=Array(values.length);let next=1;for(let k=order.length-1;k>=0;k--){const i=order[k];next=Math.min(next,values[i]*values.length/(k+1));q[i]=next;}return q;}
function createResults(seed=408){let n=seed>>>0;const u=()=>{n=(Math.imul(n,1664525)+1013904223)>>>0;return(n+1)/4294967297;},normal=()=>Math.sqrt(-2*Math.log(u()))*Math.cos(2*Math.PI*u());
 const rows=Array.from({length:320},(_,i)=>{const fc=Math.max(-5.2,Math.min(5.2,normal()*.77+(i%5===0?2.5:i%7===0?-2.4:0))),mean=10**(1.15+u()*3.55),p=10**(-Math.max(.003,Math.abs(fc)*1.25+u()*1.5-.8));
 const ratio=2**fc,targets=[2*mean*ratio/(ratio+1),2*mean/(ratio+1)];
 const profile=targets.flatMap(target=>{const values=Array.from({length:6},()=>.65+u()*.7),sum=values.reduce((a,b)=>a+b,0);return values.map(v=>v*target*6/sum);});
 return{id:`SIM-${String(i+1).padStart(4,'0')}`,fc,mean,p,q:0,annotated:i%11!==0,profile};});const qs=bh(rows.map(row=>row.p));rows.forEach((r,i)=>r.q=qs[i]);return rows;}
function classify(row,alpha=.05,lfc=1,reverse=false){const fc=reverse?-row.fc:row.fc;if(row.q>=alpha||Math.abs(fc)<lfc)return'neutral';return fc>=0?'up':'down';}
function summarize(rows,alpha,lfc,reverse=false){const counts={up:0,down:0,neutral:0};for(const row of rows)counts[classify(row,alpha,lfc,reverse)]++;return counts;}
function project(point,camera,w,h){const[x,y,z]=point;const yaw=camera.mode==='3d'?camera.yaw:0,pitch=camera.mode==='3d'?camera.pitch:0,rx=x*Math.cos(yaw)+z*Math.sin(yaw),rz=-x*Math.sin(yaw)+z*Math.cos(yaw),ry=y*Math.cos(pitch)-rz*Math.sin(pitch),depth=y*Math.sin(pitch)+rz*Math.cos(pitch),scale=camera.mode==='3d'?7/(7-depth):1,unit=Math.min(w,h)*.205*camera.zoom;return{x:w*.5+rx*unit*scale,y:h*.51-ry*unit*scale,depth,scale};}
function pick(points,x,y,r=16){let found=null,best=r;for(const p of points){const d=Math.hypot(p.x-x,p.y-y);if(d<best){best=d;found=p;}}return found;}
function route({assay='rnaseq',platforms='single',organisms='one',counts=true,confounded=false}){
 if(organisms==='many')return{state:'stop',title:'Analyze one organism at a time.',engine:'Separate analyses',detail:'The app selects an organism before modeling. It does not combine different species in one expression model.'};
 if(assay==='mixed')return{state:'stop',title:'Separate the assay types.',engine:'Joint analysis blocked',detail:'You can view and group these samples together, but mixed microarray and RNA-seq data are not analyzed jointly.'};
 if(assay==='microarray'&&platforms==='multiple')return{state:'stop',title:'Choose a single microarray platform.',engine:'limma · one GPL',detail:'Multiple GPLs can be viewed and exported together. Their probe matrices are not directly merged for limma analysis.'};
 if(assay==='microarray')return{state:'ready',title:'Follow the microarray workflow.',engine:'limma',detail:'GEO Series Matrix expression values, optional log2/normalization/vooma settings, and moderated t- or F-tests.'};
 if(!counts)return{state:'stop',title:'Provide validated raw counts.',engine:'DESeq2 requires counts',detail:'Use non-negative integer counts with sample columns matched to GSM accessions. TPM, FPKM, CPM, and log-normalized values are not valid raw-count inputs.'};
 if(platforms==='multiple'&&confounded)return{state:'stop',title:'The platform and group effects cannot be separated.',engine:'Confounded design',detail:'A platform-adjusted comparison is blocked when condition and platform are confounded. Review your sample-group design.'};
 if(platforms==='multiple')return{state:'review',title:'Check the shared count matrix and design.',engine:'DESeq2 · ~ Platform + condition',detail:'All selected GPLs must be represented in one validated count matrix. The model must have full rank. Exact official-GEO2R parity is not claimed for this extended model.'};
 return{state:'ready',title:'Follow the RNA-seq workflow.',engine:'DESeq2',detail:'Two groups use a Wald test with poscounts size factors. Three or more use an overall LRT. At least two matched samples per analyzed group are required.'};
}
return{bh,createResults,classify,summarize,project,pick,route};
});
