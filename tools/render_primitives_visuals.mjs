// README schematics of explicitly specified examples, not simulation measurements.
// Run: node tools/render_primitives_visuals.mjs
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import assert from 'node:assert/strict';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const C = {bg:'#0b1221', panel:'#121e32', edge:'#425775', text:'#f1f5fc', muted:'#b3c2d7', cyan:'#55e0d2', violet:'#b69aff', amber:'#ffc46b'};
const esc = s => String(s).replaceAll('&','&amp;').replaceAll('<','&lt;');
const txt = (x,y,s,size=22,color=C.text,anchor='start') => `<text x="${x}" y="${y}" font-size="${size}" fill="${color}" text-anchor="${anchor}">${esc(s)}</text>`;
const box = (x,y,w,h,fill=C.panel,stroke=C.edge) => `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="18" fill="${fill}" stroke="${stroke}"/>`;
const line = (x,y,u,v,color=C.edge,width=3,more='') => `<line x1="${x}" y1="${y}" x2="${u}" y2="${v}" stroke="${color}" stroke-width="${width}" ${more}/>`;
const arrow = (d,color='cyan') => `<path d="${d}" fill="none" stroke="${C[color]}" stroke-width="3" marker-end="url(#${color})"/>`;
function node(x,y,label,color=C.text,r=23) {
  return `<circle cx="${x}" cy="${y}" r="${r}" fill="${C.bg}" stroke="${color}" stroke-width="3"/>`+txt(x,y+7,label,21,color,'middle');
}
function graph(edges,positions,edgeColor=C.cyan,fresh=[]) {
  let b='';
  for(const [a,c] of edges)b+=line(...positions[a],...positions[c],edgeColor,4);
  for(const [id,p] of Object.entries(positions))b+=node(...p,id,fresh.includes(Number(id))?C.cyan:C.text);
  return b;
}
function write(name,height,title,desc,b) {
  const markers=['cyan','violet','amber'].map(k=>`<marker id="${k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto"><path d="M0 0L10 5L0 10Z" fill="${C[k]}"/></marker>`).join('');
  fs.writeFileSync(path.join(root,'docs/images',name),`<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="${height}" viewBox="0 0 1200 ${height}" role="img" aria-labelledby="title desc"><title id="title">${esc(title)}</title><desc id="desc">${esc(desc)}</desc><defs>${markers}</defs><rect width="1200" height="${height}" rx="24" fill="${C.bg}"/><g font-family="Helvetica,Arial,sans-serif">${b}</g></svg>\n`);
}

// This is exactly the labeled rule used by the README's export command.
const lhs=[[0,1],[0,2]], rhs=[[0,2],[0,3],[1,3],[2,3]];
let b=txt(36,42,'THE NETWORK ITSELF CHANGES',16,C.cyan);
b+=txt(36,87,'A rewrite replaces relationships—not positions.',34);
b+=box(26,120,508,350)+box(666,120,508,350);
b+=txt(54,160,'BEFORE · matched pattern',22,C.amber);
b+=txt(694,160,'AFTER · replacement',22,C.cyan);
b+=graph(lhs,{0:[360,280],1:[145,220],2:[145,367]},C.amber);
b+=graph(rhs,{0:[1020,285],1:[735,205],2:[735,380],3:[850,285]},C.cyan,[3]);
b+=txt(280,440,'2 relationships consumed',21,C.muted,'middle');
b+=txt(920,440,'4 produced · node 3 is new',21,C.muted,'middle');
b+=arrow('M552 285H644');
b+=txt(600,332,'apply',19,C.muted,'middle');
b+=txt(600,358,'the rule',19,C.muted,'middle');
b+=txt(36,515,'Before: (0,1), (0,2)',23,C.amber);
b+=txt(666,515,'After: (0,2), (0,3), (1,3), (2,3)',23,C.cyan);
b+=txt(36,562,'The numbers are names. The lines are the data. The drawing has no physical scale.',23);
b+=txt(36,599,'Even (0,2) is consumed and produced again: same endpoints, a new edge occurrence.',20,C.muted);
write('readme-rewrite-primitives.svg',628,'A local rewrite changes the network','Specified rule: consume edges (0,1) and (0,2); produce (0,2), (0,3), (1,3), (2,3). Node 3 is fresh. The two panels are graph layouts, not physical coordinates. This is the rule in the README export command.',b);

// A, B subdivide distinct edge occurrences. C joins one output of each.
const initial=[[0,1],[1,2]];
const A={take:[[0,1]],put:[[0,3],[3,1]]};
const B={take:[[1,2]],put:[[1,4],[4,2]]};
const E={take:[[3,1],[1,4]],put:[[3,4]]};
const key=e=>JSON.stringify(e);
function apply(state,event) {
  const out=state.map(e=>[...e]);
  for(const edge of event.take){const i=out.findIndex(e=>key(e)===key(edge));assert(i>=0);out.splice(i,1);}
  return out.concat(event.put);
}
const afterA=apply(initial,A),afterB=apply(initial,B),joined=apply(afterA,B);
assert.deepEqual(joined.map(key).sort(),apply(afterB,A).map(key).sort());
assert(E.take.every(e=>[...A.put,...B.put].some(p=>key(p)===key(e))));
const afterC=apply(joined,E);
function chain(cx,cy,labels,caption,color=C.text) {
  let q=box(cx-143,cy-64,286,132);
  q+=txt(cx,cy-31,caption,22,color,'middle');
  const gap=236/(labels.length-1),pts=labels.map((v,i)=>[cx-118+i*gap,cy+16]);
  for(let i=1;i<pts.length;i++)q+=line(...pts[i-1],...pts[i],C.edge,3);
  pts.forEach((p,i)=>{q+=node(...p,labels[i],labels[i]===3?C.cyan:labels[i]===4?C.violet:C.text,17);});
  return q;
}
b=txt(36,42,'TWO DIFFERENT GRAPHS',16,C.cyan);
b+=txt(36,87,'Possible histories are not causal dependencies.',34);
b+=txt(36,134,'MULTIWAY · each box is a whole labeled network state',21,C.muted);
b+=chain(180,333,[0,1,2],'start');
b+=chain(600,227,[0,3,1,2],'after A',C.cyan);
b+=chain(600,439,[0,1,4,2],'after B',C.violet);
b+=chain(1020,333,[0,3,1,4,2],'same final state');
b+=arrow('M330 310L448 252','cyan')+arrow('M751 249L868 309','violet');
b+=arrow('M330 357L448 416','violet')+arrow('M751 416L868 357','cyan');
b+=txt(369,265,'A',25,C.cyan)+txt(810,265,'B',25,C.violet);
b+=txt(369,425,'B',25,C.violet)+txt(810,425,'A',25,C.cyan);
b+=txt(600,558,'A splits edge (0,1). B splits edge (1,2). Either can happen first.',22,C.text,'middle');
b+=line(36,592,1164,592);
b+=txt(36,636,'CAUSAL · each circle is an event, not a spatial node',21,C.muted);
b+=node(175,718,'A',C.cyan,28)+node(175,842,'B',C.violet,28)+node(535,780,'C',C.amber,30);
b+=arrow('M210 722L492 770','cyan')+arrow('M210 836L492 791','violet');
b+=txt(327,719,'produces (3,1)',20,C.cyan,'middle');
b+=txt(327,855,'produces (1,4)',20,C.violet,'middle');
b+=txt(581,751,'C needs both new edges.',23,C.amber);
b+=txt(581,786,'It joins (3,1) and (1,4)',21);
b+=txt(581,820,'into (3,4).',21);
b+=txt(36,922,'No A → B dependency: A and B do not consume one another’s results.',23);
b+=txt(36,960,'This small example uses subdivision and joining rules. It is not a claim about quantum interference.',19,C.muted);
write('readme-history-vs-causality.svg',987,'Branching histories versus causal dependencies','Two independent subdivisions A and B can happen in either order and reach the same graph. The multiway diagram branches and merges. A later joining event C consumes an edge produced by A and an edge produced by B, so its causal graph has A to C and B to C, but no A to B edge. Specified illustrative example, not a sampled physical process.',b);

const pathEdges=[[0,1],[1,2],[2,3]],cycleEdges=[...pathEdges,[3,0]];
// Both candidates have the same named outer endpoints and connectivity.
function connected(edges,start,end){const seen=new Set([start]);let more=true;while(more){more=false;for(const [a,c] of edges){if(seen.has(a)&&!seen.has(c)){seen.add(c);more=true;}if(seen.has(c)&&!seen.has(a)){seen.add(a);more=true;}}}return seen.has(end);}
assert(connected(pathEdges,0,3)&&connected(cycleEdges,0,3));
b=txt(36,42,'THE MISSING FIBER CONSTRUCTION',16,C.cyan);
b+=txt(36,87,'Which details disappear in a coarser description?',33);
b+=box(26,128,485,456,C.panel,C.violet);
b+=txt(53,171,'TWO DIFFERENT DETAILED STATES',20,C.violet);
b+=graph(pathEdges,{0:[95,255],1:[212,222],2:[328,282],3:[440,255]},C.violet);
b+=txt(268,334,'an open chain',21,C.muted,'middle');
b+=graph(cycleEdges,{0:[95,440],1:[212,390],2:[328,390],3:[440,440]},C.violet);
b+=txt(268,552,'a closed loop',21,C.muted,'middle');
b+=arrow('M525 278L748 338')+arrow('M525 468L748 408');
b+=txt(637,211,'CHOSEN OBSERVER',18,C.cyan,'middle');
b+=txt(637,244,'keeps endpoints',20,C.text,'middle');
b+=txt(637,513,'and whether they',20,C.text,'middle');
b+=txt(637,541,'are connected',20,C.text,'middle');
b+=box(770,240,404,268);
b+=txt(972,288,'ONE COARSE DESCRIPTION x',19,C.cyan,'middle');
b+=line(852,362,1087,362,C.cyan,4)+node(852,362,0)+node(1087,362,3);
b+=txt(972,443,'“0 and 3 are connected”',23,C.text,'middle');
b+=txt(36,633,'Both detailed states belong to the fiber over x under this projection.',24);
b+=txt(36,675,'But they need not behave alike. Grouping them does not make them gauge-equivalent.',22,C.amber);
b+=txt(36,712,'The research task: derive a useful projection and transport from the rewriting process itself.',21,C.muted);
write('readme-fiber-projection.svg',742,'A fiber as alternatives over one coarse description','A chosen coarse observer keeps endpoint identities 0 and 3 and whether they are connected. It maps both a four-node path and a four-node cycle to the same description. Both lie in its set-theoretic fiber. This does not imply dynamical or gauge equivalence. The projection is an explicit illustrative assumption, not a derived physical observer.',b);
console.log('Rendered three specified primitive and projection examples.');
