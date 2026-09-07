import {BUDGET,indexEvolution,graphOf,neighborhood,lineage,causalView,branchialView,metricProjection} from './model.mjs';
const $=id=>document.getElementById(id),canvas=$('canvas'),ctx=canvas.getContext('2d');
let index,state,graph,vertex,view='hypergraph',nodes=[],links=[],hits=[],timer,angle=.4,tilt=.35,zoom=1,drag;
const descriptions={
  hypergraph:'Selected state as an incidence diagram. Click a vertex to inspect its neighborhood. Each diamond is one ordered hyperedge; selection lists its exact slots.',
  causal:'Causal ancestors of the event producing this state. Arrows point from prerequisite to dependent event. Click an event to inspect its output state; Play lineage animates a rewrite history.',
  branchial:'States at the selected rewrite depth. Links are exported branchial event relations mapped to their output states. Click a state to synchronize all views.',
  projection:'Best-effort Euclidean coordinates from shortest-path distances using classical metric MDS. Drag to rotate, scroll to zoom. No springs or molecular shape template.'
};
const notices={
  hypergraph:'Local radius uses the hypergraph 2-section. Hyperedge order is preserved in the inspector. Diagram positions carry no distance meaning.',
  causal:'Causal dependency does not establish spatial locality. This export contains no amplitudes or entanglement measure. Rewrite depth is a scheduler coordinate, not recovered physical time.',
  branchial:'Only relations whose two output states share this depth are displayed. Branchial adjacency alone does not quantify entanglement.',
  projection:'This is a connected local patch in graph-distance units. Cropping can lengthen shortest paths. Stress measures 3-D distance distortion; it does not establish emergent 3-D space.'
};
function option(value,label){const o=document.createElement('option');o.value=value;o.textContent=label;return o;}
function selectState(id){
  state=index.states.get(Number(id));if(!state)return;graph=graphOf(state);
  if(!graph.has(vertex))vertex=graph.keys().next().value;
  $('step').value=state.step;$('step-label').textContent=`${state.step} / ${$('step').max}`;
  const slice=index.slices.get(state.step)||[],shown=slice.slice(0,500);if(!shown.some(s=>s.id===state.id))shown.push(state);
  $('state').replaceChildren(...shown.map(s=>option(s.id,`${s.id} · class ${s.canonical_id}`)));$('state').value=state.id;
  const vs=[...graph.keys()].slice(0,500);if(!vs.includes(vertex))vs.push(vertex);
  $('vertex').replaceChildren(...vs.map(v=>option(v,v)));$('vertex').value=vertex;
  const p=state.intrinsic_geometry||{},dimension=k=>p[k+'_plateau']?`${p[k+'_plateau'].mean.toFixed(2)} (${p[k+'_plateau'].first_scale}–${p[k+'_plateau'].last_scale})`:'no plateau';
  const metrics={State:state.id,'Canonical class':state.canonical_id,Vertices:graph.size,Hyperedges:state.edges.length,'Volume dimension':dimension('volume_dimension'),'Spectral dimension':dimension('spectral_dimension')};
  $('metrics').replaceChildren(...Object.entries(metrics).flatMap(([k,v])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=k;dd.textContent=v;return[dt,dd];}));
  $('selection').textContent=`State ${state.id}, vertex ${vertex}\nNeighbors: ${[...(graph.get(vertex)||[])].slice(0,30).join(', ')}`;
  build();
}
function radial(ids,distance){
  const rings=new Map();for(const id of ids){const r=distance?.get(id)??1;if(!rings.has(r))rings.set(r,[]);rings.get(r).push(id);}
  return ids.map(id=>{const r=distance?.get(id)??1,list=rings.get(r),a=2*Math.PI*list.indexOf(id)/list.length;return{id,x:r*Math.cos(a),y:r*Math.sin(a),z:0,color:id===vertex?'#ffd16d':'#76c7fb',kind:'vertex',label:String(id)};});
}
function edgeDelta(current,parent){
  const remaining=new Map();for(const edge of parent?.edges||[]){const key=JSON.stringify(edge);remaining.set(key,(remaining.get(key)||0)+1);}
  const edges=current.edges.map((edge,i)=>{const key=JSON.stringify(edge),old=remaining.get(key)||0;if(old)remaining.set(key,old-1);return{edge,key:`h${i}`,added:!!parent&&!old};});
  for(const[key,count]of remaining)for(let i=0;i<count;i++)edges.push({edge:JSON.parse(key),key:`r${edges.length}`,removed:true});return edges;
}
function build(){
  nodes=[];links=[];zoom=1;
  $('view-title').textContent=document.querySelector(`[data-view=${view}]`).textContent;$('description').textContent=descriptions[view];$('notice').textContent=notices[view];
  let detail='',truncated=false;
  if(view==='hypergraph'||view==='projection'){
    const local=neighborhood(graph,vertex,Number($('radius').value),view==='projection'?BUDGET.projection:Math.floor(BUDGET.nodes/2));
    truncated=local.truncated;
    if(view==='projection'){
      const projected=metricProjection(graph,local.ids);
      nodes=local.ids.map((id,i)=>({id,x:projected.points[i][0],y:projected.points[i][1],z:projected.points[i][2],kind:'vertex',label:String(id),color:id===vertex?'#ffd16d':'#76c7fb'}));
      detail=`MDS stress ${(projected.stress*100).toFixed(1)}% · ${projected.eigenvalues.length} positive display modes · ${nodes.length}/${graph.size} vertices`;
      // Coordinates carry the geometry; edges are deliberately omitted in this view.
    }else{
      nodes=radial(local.ids,local.distance);const byId=new Map(nodes.map(n=>[n.id,n]));
      const event=index.incoming.get(state.id),parent=event?index.states.get(event.input):undefined;
      for(const e of edgeDelta(state,parent)){
        const included=e.edge.filter(v=>byId.has(v));if(included.length!==e.edge.length)continue;
        if(nodes.length>=BUDGET.nodes||links.length+included.length>BUDGET.edges){truncated=true;break;}
        const center=included.map(v=>byId.get(v)),j=nodes.length;
        const h={id:e.key,x:center.reduce((s,p)=>s+p.x,0)/center.length+.07*Math.sin(j),y:center.reduce((s,p)=>s+p.y,0)/center.length+.07*Math.cos(j),z:0,kind:'hyperedge',label:e.key,color:e.removed?'#ed8b8b':e.added?'#72e6ad':'#bc9bff',edge:e.edge};nodes.push(h);
        included.forEach((v,slot)=>links.push({a:v,b:e.key,color:h.color,dashed:e.removed,slot}));
      }
      detail=`${local.ids.length}/${graph.size} vertices in selected patch · exact ordered hyperedge incidence`;
    }
  }else if(view==='branchial'){
    const b=branchialView(index,state);truncated=b.truncated;nodes=radial(b.ids).map(n=>({...n,kind:'state',color:n.id===state.id?'#ffd16d':'#76c7fb'}));links=b.edges.map(([a,b])=>({a,b}));detail=`${nodes.length}/${b.total} states at depth ${state.step}`;
  }else{
    const event=index.incoming.get(state.id),c=causalView(index,event?.id);truncated=c.truncated;
    const layers=new Map();for(const id of c.ids){const e=index.events.get(id),depth=index.states.get(e.output).step;if(!layers.has(depth))layers.set(depth,[]);layers.get(depth).push(id);}
    for(const[depth,ids]of layers)ids.forEach((id,i)=>nodes.push({id,x:(i-(ids.length-1)/2)*1.3,y:depth*1.5,z:0,kind:'event',label:`e${id}`,color:id===event.id?'#ffd16d':'#b7a06d'}));
    links=c.edges.map(([a,b])=>({a,b,arrow:true,color:'#b7a06d'}));detail=event?`Causal past of event ${event.id} · prerequisite → dependent`:'Initial state: no producing event';
  }
  $('budget').textContent=`${nodes.length}/${BUDGET.nodes} displayed nodes · ${links.length}/${BUDGET.edges} links. Projection cap: ${BUDGET.projection}. ${truncated?'Budget reached: view is truncated.':'Within display budget.'} File limit: 32 MiB; whole export loaded in memory.`;
  $('status').textContent=detail+' · click to inspect'+(view==='projection'?' · drag to rotate':'');draw();
}
function draw(){
  const rect=canvas.getBoundingClientRect(),w=rect.width,h=rect.height,dpr=Math.min(devicePixelRatio,2);canvas.width=w*dpr;canvas.height=h*dpr;ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);hits=[];
  if(!nodes.length)return;
  const points=nodes.map(n=>{if(view!=='projection')return{...n,px:n.x,py:n.y};const x=n.x*Math.cos(angle)+n.z*Math.sin(angle),z=-n.x*Math.sin(angle)+n.z*Math.cos(angle);return{...n,px:x,py:n.y*Math.cos(tilt)-z*Math.sin(tilt)};});
  const xs=points.map(p=>p.px),ys=points.map(p=>p.py),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys),scale=Math.min((w-130)/Math.max(xmax-xmin,2),(h-120)/Math.max(ymax-ymin,2))*zoom;
  for(const p of points){p.sx=w/2+(p.px-(xmin+xmax)/2)*scale;p.sy=h/2-(p.py-(ymin+ymax)/2)*scale;}
  const byId=new Map(points.map(p=>[p.id,p]));
  for(const l of links){const a=byId.get(l.a),b=byId.get(l.b);if(!a||!b)continue;ctx.strokeStyle=l.color||'#45617e';ctx.lineWidth=1;ctx.setLineDash(l.dashed?[4,4]:[]);ctx.beginPath();ctx.moveTo(a.sx,a.sy);ctx.lineTo(b.sx,b.sy);ctx.stroke();ctx.setLineDash([]);
    if(l.arrow){const angle=Math.atan2(b.sy-a.sy,b.sx-a.sx),x=b.sx-10*Math.cos(angle),y=b.sy-10*Math.sin(angle);ctx.beginPath();ctx.moveTo(x-8*Math.cos(angle-.4),y-8*Math.sin(angle-.4));ctx.lineTo(x,y);ctx.lineTo(x-8*Math.cos(angle+.4),y-8*Math.sin(angle+.4));ctx.stroke();}}
  for(const p of points){ctx.fillStyle=p.color;ctx.beginPath();if(p.kind==='hyperedge'){ctx.moveTo(p.sx,p.sy-6);ctx.lineTo(p.sx+6,p.sy);ctx.lineTo(p.sx,p.sy+6);ctx.lineTo(p.sx-6,p.sy);ctx.closePath();}else ctx.arc(p.sx,p.sy,p.color==='#ffd16d'?8:5,0,Math.PI*2);ctx.fill();if(points.length<65||p.color==='#ffd16d'){ctx.fillStyle='#d7e6f9';ctx.font='11px system-ui';ctx.fillText(p.label,p.sx+9,p.sy-8);}hits.push(p);}
}
function stop(){clearInterval(timer);timer=undefined;$('play').textContent='Play lineage';}
document.querySelectorAll('[data-view]').forEach(button=>button.onclick=()=>{view=button.dataset.view;document.querySelectorAll('[data-view]').forEach(b=>b.setAttribute('aria-selected',b===button));if(state)build();});
$('step').oninput=()=>{stop();selectState(index.slices.get(Number($('step').value))?.[0]?.id);};
$('state').onchange=()=>{stop();selectState($('state').value);};
$('vertex').onchange=()=>{vertex=Number($('vertex').value);selectState(state.id);};$('radius').onchange=()=>build();
$('play').onclick=()=>{
  if(timer){stop();return;}let target=state.id;while(index.children.get(target)?.length)target=index.children.get(target)[0];const path=lineage(index,target);let at=0;
  $('play').textContent='Pause';selectState(path[at++]);timer=setInterval(()=>{if(at>=path.length){stop();return;}selectState(path[at++]);},950);
};
canvas.onpointerdown=e=>{drag={x:e.clientX,y:e.clientY,moved:false};canvas.setPointerCapture(e.pointerId);};
canvas.onpointermove=e=>{if(!drag||view!=='projection')return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;if(Math.abs(dx)+Math.abs(dy)>2)drag.moved=true;angle+=dx*.01;tilt+=dy*.01;drag.x=e.clientX;drag.y=e.clientY;draw();};
canvas.onpointerup=e=>{
  if(!drag?.moved){const r=canvas.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top,p=hits.findLast(p=>Math.hypot(p.sx-x,p.sy-y)<13);if(p){stop();if(p.kind==='vertex'){vertex=p.id;selectState(state.id);}else if(p.kind==='state')selectState(p.id);else if(p.kind==='event'){const ev=index.events.get(p.id);selectState(ev.output);$('selection').textContent=`Event ${ev.id}, rule ${ev.rule}\nState ${ev.input} → ${ev.output}`;}else $('selection').textContent=`Hyperedge ${p.id}\nOrdered slots: [${p.edge.join(', ')}]`;}}
  drag=undefined;
};
canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(.3,Math.min(4,zoom*Math.exp(-e.deltaY*.001)));draw();};
new ResizeObserver(draw).observe($('viewport'));
function load(data){stop();index=indexEvolution(data);$('step').max=Math.max(...index.slices.keys());selectState(data.states[0].id);}
$('file').onchange=async()=>{try{const file=$('file').files[0];if(!file)return;if(file.size>BUDGET.fileBytes)throw Error('Export exceeds 32 MiB. Export a bounded window first.');load(JSON.parse(await file.text()));}catch(e){$('status').textContent=e.message;}};
try{const response=await fetch(new URLSearchParams(location.search).get('data')||'../data/example-evolution.json');if(!response.ok)throw Error(`Export HTTP ${response.status}`);const text=await response.text();if(new Blob([text]).size>BUDGET.fileBytes)throw Error('Export exceeds 32 MiB. Export a bounded window first.');load(JSON.parse(text));}catch(e){$('status').textContent=e.message;}
