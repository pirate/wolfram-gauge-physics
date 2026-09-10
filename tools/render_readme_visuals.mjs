// README-only illustrations. No simulation state is changed.
// Run from the repository root: node tools/render_readme_visuals.mjs
import fs from 'node:fs';
import assert from 'node:assert/strict';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const out = path.join(root, 'docs/images');
const C = {bg:'#0b1221', panel:'#121e32', edge:'#354764', text:'#f1f5fc', muted:'#a8b8d0', cyan:'#55e0d2', violet:'#b69aff', amber:'#ffc46b'};
const esc = s => String(s).replaceAll('&','&amp;').replaceAll('<','&lt;');
const text = (x,y,s,size=18,fill=C.text,anchor='start',extra='') => `<text x="${x}" y="${y}" fill="${fill}" font-size="${size}" text-anchor="${anchor}" ${extra}>${esc(s)}</text>`;
const line = (x1,y1,x2,y2,color=C.edge,width=2,extra='') => `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="${width}" ${extra}/>`;
const rect = (x,y,w,h,fill=C.panel,stroke=C.edge) => `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="18" fill="${fill}" stroke="${stroke}"/>`;
function svg(name,height,title,description,body) {
  const head = `<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="${height}" viewBox="0 0 1200 ${height}" role="img" aria-labelledby="title desc"><title id="title">${esc(title)}</title><desc id="desc">${esc(description)}</desc><defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="${C.cyan}"/></marker></defs><g font-family="Helvetica,Arial,sans-serif"><rect width="1200" height="${height}" rx="24" fill="${C.bg}"/>`;
  fs.writeFileSync(path.join(out,name), head+body+'</g></svg>');
}

// Exploded schematic: projection onto the base is not an extra spatial axis.
let fiber = text(36,43,'WHAT IS A FIBER?',15,C.cyan,'start','letter-spacing="2"');
fiber += text(36,83,'One location. A separate space of internal possibilities.',29);
fiber += text(36,116,'The base graph connects locations. A fiber graph describes the internal structure at one location.',17,C.muted);
fiber += rect(24,144,1152,250);
fiber += text(48,179,'FIBERS',16,C.cyan,'start','letter-spacing="1.5"');
const fiberXs=[250,600,950];
for(let i=0;i<3;i++) {
  const x=fiberXs[i],pts=[[x,220],[x+53,308],[x-53,308]];
  for(let k=0;k<3;k++) fiber+=line(...pts[k],...pts[(k+1)%3],C.violet,3);
  pts.forEach(([px,py],k)=>{fiber+=`<circle cx="${px}" cy="${py}" r="17" fill="${C.violet}"/>`;fiber+=text(px,py+6,k,17,C.bg,'middle');});
  fiber+=text(x,359,`fiber at ${'ABC'[i]}`,18,C.text,'middle');
  fiber+=line(x,376,x,473,C.muted,2,'stroke-dasharray="5 6"');
  fiber+=text(x+12,428,'belongs to',14,C.muted);
}
fiber+=line(275,504,575,504,C.cyan,4);fiber+=line(625,504,925,504,C.cyan,4);
fiber+=text(48,510,'BASE GRAPH',15,C.cyan,'start','letter-spacing="1"');
fiberXs.forEach((x,i)=>{fiber+=`<circle cx="${x}" cy="504" r="26" fill="${C.bg}" stroke="${C.cyan}" stroke-width="3"/>`;fiber+=text(x,512,'ABC'[i],23,C.text,'middle');});
fiber+=text(425,536,'spatial relationship',15,C.muted,'middle');fiber+=text(775,536,'spatial relationship',15,C.muted,'middle');
fiber+=rect(24,569,1152,94,'#152a36','#28736d');
fiber+=text(46,603,'A fiber node is not another location in ordinary space.',23);
fiber+=text(46,634,'Vertical separation is for the drawing only. Each entire triangle belongs to the single base node below it.',17,C.muted);
svg('readme-what-is-a-fiber.svg',687,'A fiber is an internal graph attached to one base node','Three connected base nodes A, B, and C each have a three-vertex internal triangle graph. Dashed projection lines associate every triangle with its base node. The vertical separation is diagrammatic, not an extra spatial dimension.',fiber);

// A specified, legal C3 connection: two identity maps and one +1 rotation.
// The last fiber repeats the first base location; it is not a fourth vertex.
let a = text(36,43,'TRANSPORT AROUND A LOOP',15,C.cyan,'start','letter-spacing="2"');
a += text(36,82,'Return to the same place. Keep a different internal alignment.',28);
a += text(36,113,'A chosen three-link connection with a triangle fiber; not a recorded trajectory.',17,C.muted);
const centers = [146,449,752,1055], names=['A · start','B','C','A · return'];
const colors=[C.cyan,C.violet,C.amber];
for(let j=0;j<4;j++) {
  const x=centers[j];
  a+=rect(x-107,147,214,256);
  a+=text(x,182,names[j],21,C.text,'middle');
  const pts=[[x,229],[x+65,337],[x-65,337]];
  for(let k=0;k<3;k++) a+=line(...pts[k],...pts[(k+1)%3],C.edge,4);
  for(let k=0;k<3;k++) {
    const token=j===3?(k+2)%3:k, [px,py]=pts[k];
    a+=`<circle cx="${px}" cy="${py}" r="21" fill="${colors[token]}"/>`;
    a+=text(px,py+7,token,21,C.bg,'middle','font-weight="700"');
    a+=text(px,py+(k===0?-31:39),`slot ${k}`,13,C.muted,'middle');
  }
  if(j<3) {
    a+=line(x+111,277,centers[j+1]-113,277,C.cyan,3,'marker-end="url(#arrow)"');
    a+=text((x+centers[j+1])/2,254,j===2?'+1':'identity',15,C.cyan,'middle');
  }
}
a+=rect(36,428,1128,83,'#152a36','#28736d');
a+=text(58,461,'Follow marker 0:',20,C.cyan);
a+=text(248,461,'it starts in slot 0 and returns in slot 1.',20);
a+=text(58,488,'The loop’s net map is a rotation. This net map is its holonomy.',18,C.muted);
svg('readme-fiber-holonomy.svg',535,'Holonomy in a triangle fiber','Four panels show A, B, C, then A again. Two identity transports followed by a rotation move marker 0 from slot 0 to slot 1 at the same base vertex.',a);

// Read the saved raw-link witness and reconstruct its documented mesh ordering.
const visualData=JSON.parse(fs.readFileSync(path.join(out,'readme-visual-data.json'),'utf8'));
const witness=visualData.reaction;
const side=4, vertex=(x,y)=>((y+side)%side)*side+(x+side)%side;
const faces=[], edgeSet=new Map();
for(let y=0;y<side;y++) for(let x=0;x<side;x++) {
  const [a,b,c,d]=[vertex(x,y),vertex(x+1,y),vertex(x+1,y+1),vertex(x,y+1)];
  for(let f of [[a,b,c],[a,c,d]]) {
    const k=f.indexOf(Math.min(...f)); f=[...f.slice(k),...f.slice(0,k)]; faces.push(f);
    for(let k=0;k<3;k++) {const e=[f[k],f[(k+1)%3]].sort((a,b)=>a-b);edgeSet.set(e.join(','),e);}
  }
}
const edges=[...edgeSet.values()].sort((a,b)=>a[0]-b[0]||a[1]-b[1]);
const ids=new Map(edges.map((e,i)=>[e.join(','),i]));
const n=witness.n, mod=x=>(x%n+n)%n;
const mul=(a,b)=>mod(a%n+(a<n?1:-1)*(b%n))+( (a<n)===(b<n)?0:n );
const inv=a=>a<n?mod(-a):a;
const kind=a=>a===0?0:a<n?2:1;
const raws=[witness.before_raw_links,witness.between_raw_links,witness.after_raw_links];
const qs=[witness.before_charge,witness.between_charge,witness.after_charge];
assert.equal(raws.length,3);
for(let t=0;t<3;t++) {
  const q=faces.map(f=> {let h=0; for(let k=0;k<3;k++){const u=f[k],v=f[(k+1)%3];let g=raws[t][ids.get([u,v].sort((a,b)=>a-b).join(','))];if(u>v)g=inv(g);h=mul(g,h);} return kind(h);});
  assert.deepEqual(q,qs[t],'face colors must be reconstructed from raw links');
  assert.equal(q.reduce((a,b)=>a+b),witness.conserved_charge);
}
const changed=[[],edges.filter((e,i)=>raws[0][i]!==raws[1][i]),edges.filter((e,i)=>raws[1][i]!==raws[2][i])];
const selected=[6,1,0,3,2], stateColors=['#203641',C.violet,C.amber], type=['E','R','Z'];
let b=text(36,43,'TWO REACTIONS · ONE NET TRANSFER',15,C.cyan,'start','letter-spacing="2"');
b+=text(36,83,'A shared face connects two local reactions.',30);
b+=text(36,115,'Saved C₅ shared-link states • five-face patch unwrapped from a periodic mesh',17,C.muted);
for(let t=0;t<3;t++) {
  const px=24+t*392; b+=rect(px,143,368,332);
  b+=text(px+20,181,['01  Before','02  After first reaction','03  After second reaction'][t],22);
  b+=text(px+20,210,['Neighboring reaction is inactive','Shared face becomes flat','Shared face returns to R'][t],16,C.muted);
  const point=v=>{let x=v%side,y=Math.floor(v/side); if(x===3)x=-1;return [px+130+(x-.5*y)*82,275+y*Math.sqrt(3)/2*82];};
  for(const f of selected) {
    const pts=faces[f].map(point), q=qs[t][f];
    b+=`<polygon points="${pts.map(p=>p.join(',')).join(' ')}" fill="${stateColors[q]}" fill-opacity="${q===0?1:.84}" stroke="${C.bg}" stroke-width="3"/>`;
    const cx=pts.reduce((s,p)=>s+p[0],0)/3,cy=pts.reduce((s,p)=>s+p[1],0)/3;
    b+=text(cx,cy+3,`${type[q]} · ${q}`,17,q===0?C.cyan:C.bg,'middle','font-weight="700"');
    b+=text(cx,cy+19,`f${f}`,12,q===0?C.muted:'#26334d','middle');
  }
  for(const e of changed[t]) b+=line(...point(e[0]),...point(e[1]),C.cyan,4);
  const ledger=[[3,'source'],[0,'shared'],[6,'destination']];
  ledger.forEach(([f,label],i)=>{const x=px+65+i*117;b+=text(x,405,label,14,C.muted,'middle');b+=text(x,440,qs[t][f],30,stateColors[qs[t][f]]==='#203641'?C.cyan:stateColors[qs[t][f]],'middle','font-weight="700"');});
}
b+=rect(24,497,1152,90,'#152a36','#28736d');
b+=text(45,534,'NET CHANGE',14,C.cyan,'start','letter-spacing="1.5"');
b+=text(224,534,'source −1',22,C.amber); b+=line(353,527,505,527,C.cyan,3,'marker-end="url(#arrow)"');
b+=text(526,534,'destination +1',22,C.amber);b+=text(1150,534,'Total Q = 52 throughout',19,C.text,'end');
b+=text(45,565,'Cyan lines: changed raw links. E / R / Z: identity / reflection / rotation. The arrow shows net balance, not a path.',16,C.muted);
svg('readme-reaction-transfer.svg',612,'Charge transfer through two neighboring reactions','The saved raw-link witness keeps total charge 52. Face 3 changes 2 to 2 to 1; shared face 0 changes 1 to 0 to 1; face 6 changes 1 to 2 to 2. Cyan edges are the actual changed links.',b);
// Two orders of the recorded circuit actions, not invented spatial paths.
const [A,B]=visualData.transport.generator_permutations;
assert.equal(visualData.transport.same_face_classes,true);
assert.equal(visualData.transport.orders_gauge_equivalent,false);
const sequences=[[0,A[0],B[A[0]]],[0,B[0],A[B[0]]]];
assert.deepEqual(sequences,[[0,0,1],[0,1,2]]);
let c=text(36,43,'ORDER OF TRANSPORT',15,C.cyan,'start','letter-spacing="2"');
c+=text(36,83,'Same circuits. Different order. Different connection.',30);
c+=text(36,115,'Recorded action on complete gauge states, starting from state 0 in both cases.',17,C.muted);
for(let row=0;row<2;row++) {
  const y=224+166*row; c+=rect(24,y-72,1152,139);
  c+=text(48,y-32,row===0?'A THEN B':'B THEN A',16,C.muted,'start','letter-spacing="1.5"');
  const xs=[238,614,990];
  sequences[row].forEach((state,k)=>{
    const color=k===2?(row===0?C.violet:C.amber):C.text;
    c+=`<circle cx="${xs[k]}" cy="${y}" r="36" fill="${C.bg}" stroke="${color}" stroke-width="2"/>`;
    c+=text(xs[k],y+11,state,31,color,'middle');
    c+=text(xs[k],y+58,k===2?'final connection':'gauge state',14,C.muted,'middle');
    if(k<2){const isA=(row===0)===(k===0),color=isA?C.cyan:C.amber;
      c+=line(xs[k]+50,y,xs[k+1]-54,y,color,3);
      c+=`<path d="M ${xs[k+1]-54} ${y-7} L ${xs[k+1]-40} ${y} L ${xs[k+1]-54} ${y+7}" fill="${color}"/>`;
      c+=text((xs[k]+xs[k+1])/2,y-17,isA?'circuit A':'circuit B',19,color,'middle');
    }
  });
}
c+=rect(24,485,1152,78,'#152a36','#28736d');
c+=text(45,518,'Every defect returns to its face. Every face keeps the same class.',21);
c+=text(45,547,'The difference is in the relationships between loops. A face-color map does not show it.',18,C.muted);
svg('readme-transport-memory.svg',588,'Transport order changes the complete connection','Circuit A then B takes gauge state 0 through 0 to 1. Circuit B then A takes 0 through 1 to 2. The two final connections have the same face classes but are not gauge equivalent.',c);
console.log('Rendered four README figures; all 96 face charges agree with saved raw links and both circuit orders match the saved action.');
