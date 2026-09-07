// Display adapters only: these never change evolution or infer physical coordinates.
export const BUDGET={nodes:240,edges:1200,projection:96,fileBytes:32*1024*1024};
export function indexEvolution(data){
  if(!Array.isArray(data.states)||!data.states.length||!Array.isArray(data.events))throw Error('Expected an evolution export with states and events.');
  const states=new Map(data.states.map(s=>[s.id,s])),events=new Map(data.events.map(e=>[e.id,e]));
  const incoming=new Map(),children=new Map(),slices=new Map(),causalParents=new Map();
  for(const s of states.values()){if(!slices.has(s.step))slices.set(s.step,[]);slices.get(s.step).push(s);}
  for(const e of events.values()){incoming.set(e.output,e);if(!children.has(e.input))children.set(e.input,[]);children.get(e.input).push(e.output);}
  for(const [a,b]of data.causal_edges||[]){if(!events.has(a)||!events.has(b))continue;if(!causalParents.has(b))causalParents.set(b,[]);causalParents.get(b).push(a);}
  return{data,states,events,incoming,children,slices,causalParents};
}
export function graphOf(state){
  const graph=new Map();for(const edge of state.edges){for(const v of edge)if(!graph.has(v))graph.set(v,new Set());
    for(let i=0;i<edge.length;i++)for(let j=i+1;j<edge.length;j++)if(edge[i]!==edge[j]){graph.get(edge[i]).add(edge[j]);graph.get(edge[j]).add(edge[i]);}}
  return graph;
}
export function neighborhood(graph,root,radius=99,budget=BUDGET.nodes){
  if(!graph.has(root))return{ids:[],distance:new Map(),truncated:false};
  const ids=[root],distance=new Map([[root,0]]);let truncated=false;
  for(let i=0;i<ids.length;i++){const v=ids[i];if(distance.get(v)>=radius)continue;
    for(const w of graph.get(v))if(!distance.has(w)){if(ids.length>=budget){truncated=true;continue;}distance.set(w,distance.get(v)+1);ids.push(w);}}
  return{ids,distance,truncated};
}
export function lineage(index,stateId){
  const result=[],seen=new Set();while(index.states.has(stateId)&&!seen.has(stateId)){seen.add(stateId);result.push(stateId);const event=index.incoming.get(stateId);if(!event)break;stateId=event.input;}return result.reverse();
}
export function causalView(index,eventId){
  const ids=[],seen=new Set(),queue=[eventId];let truncated=false;
  for(let i=0;i<queue.length;i++){const id=queue[i];if(seen.has(id)||!index.events.has(id))continue;if(ids.length>=BUDGET.nodes){truncated=true;break;}seen.add(id);ids.push(id);queue.push(...(index.causalParents.get(id)||[]));}
  const edges=(index.data.causal_edges||[]).filter(([a,b])=>seen.has(a)&&seen.has(b));return{ids,edges:edges.slice(0,BUDGET.edges),truncated:truncated||edges.length>BUDGET.edges};
}
export function branchialView(index,state){
  const slice=index.slices.get(state.step)||[],pairs=[],graph=new Map(slice.map(s=>[s.id,new Set()]));
  // Export endpoints are EVENT IDs. Map to outputs, then restrict both to this slice.
  for(const[a,b]of index.data.branchial_edges||[]){const left=index.events.get(a)?.output,right=index.events.get(b)?.output;if(graph.has(left)&&graph.has(right)&&left!==right){pairs.push([left,right]);graph.get(left).add(right);graph.get(right).add(left);}}
  const local=neighborhood(graph,state.id,99,BUDGET.nodes),ids=local.ids,seen=new Set(ids);
  for(const s of slice)if(ids.length<BUDGET.nodes&&!seen.has(s.id)){ids.push(s.id);seen.add(s.id);}
  const edges=pairs.filter(([a,b])=>seen.has(a)&&seen.has(b));return{ids,edges:edges.slice(0,BUDGET.edges),total:slice.length,truncated:slice.length>ids.length||edges.length>BUDGET.edges};
}
// Classical MDS of shortest-path distances on a bounded connected induced patch.
// Three positive modes supply display coordinates; stress reports distance distortion.
export function metricProjection(graph,ids){
  const n=ids.length;if(!n)return{points:[],stress:0,eigenvalues:[]};const allowed=new Set(ids),dist=[];
  for(const root of ids){const d=new Map([[root,0]]),q=[root];for(let i=0;i<q.length;i++)for(const w of graph.get(q[i])||[])if(allowed.has(w)&&!d.has(w)){d.set(w,d.get(q[i])+1);q.push(w);}
    if(d.size!==n)throw Error('Projection requires one connected patch; choose a vertex in a component.');dist.push(ids.map(v=>d.get(v)));}
  const means=dist.map(row=>row.reduce((s,d)=>s+d*d,0)/n),mean=means.reduce((a,b)=>a+b,0)/n;
  const matrix=dist.map((row,i)=>row.map((d,j)=>-.5*(d*d-means[i]-means[j]+mean)));
  const shift=Math.max(...matrix.map(row=>row.reduce((a,b)=>a+Math.abs(b),0)),1),modes=[],eigenvalues=[];
  for(let k=0;k<Math.min(3,n-1);k++){
    let v=ids.map((_,i)=>Math.sin((i+1)*(k+1)*1.731)+Math.cos((i+1)*.619));
    for(let iteration=0;iteration<250;iteration++){
      const next=matrix.map((row,i)=>row.reduce((s,x,j)=>s+x*v[j],shift*v[i]));
      const center=next.reduce((a,b)=>a+b,0)/n;for(let i=0;i<n;i++)next[i]-=center;
      for(const u of modes){const dot=next.reduce((s,x,i)=>s+x*u[i],0);for(let i=0;i<n;i++)next[i]-=dot*u[i];}
      const norm=Math.hypot(...next);if(norm<1e-12)break;v=next.map(x=>x/norm);
    }
    const norm=Math.hypot(...v);v=v.map(x=>x/norm);
    const eigenvalue=v.reduce((sum,x,i)=>sum+x*matrix[i].reduce((s,b,j)=>s+b*v[j],0),0);
    if(eigenvalue<1e-8)break;modes.push(v);eigenvalues.push(eigenvalue);
  }
  const points=ids.map((_,i)=>[0,1,2].map(k=>modes[k]?modes[k][i]*Math.sqrt(eigenvalues[k]):0));let residual=0,total=0;
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){const embedded=Math.hypot(...points[i].map((x,k)=>x-points[j][k]));residual+=(embedded-dist[i][j])**2;total+=dist[i][j]**2;}
  return{points,stress:Math.sqrt(residual/(total||1)),eigenvalues};
}
