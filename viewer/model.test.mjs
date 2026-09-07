import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {BUDGET,indexEvolution,graphOf,neighborhood,lineage,causalView,branchialView,metricProjection} from './model.mjs';
test('event and state identifiers are distinct in causal and branchial exports',()=>{
  const idx=indexEvolution({states:[{id:10,step:0},{id:20,step:1},{id:30,step:1},{id:40,step:2}],events:[{id:1,input:10,output:20},{id:2,input:10,output:30},{id:3,input:20,output:40}],causal_edges:[[1,3]],branchial_edges:[[1,2],[2,3]]});
  assert.deepEqual(branchialView(idx,idx.states.get(20)).edges,[[20,30]]);
  assert.deepEqual(causalView(idx,3).edges,[[1,3]]);assert.deepEqual(lineage(idx,40),[10,20,40]);
});
test('bounded neighborhoods stay connected and preserve a selected root',()=>{
  const g=graphOf({edges:Array.from({length:1000},(_,i)=>[i,i+1])}),n=neighborhood(g,500,99,20);
  assert.equal(n.ids.length,20);assert.equal(n.ids[0],500);assert.equal(n.truncated,true);assert.ok(metricProjection(g,n.ids).stress<1e-6);
});
test('MDS recovers a line and refuses disconnected distance pairs',()=>{
  const g=graphOf({edges:[[0,1],[1,2],[2,3],[3,4],[8,9]]}),p=metricProjection(g,[0,1,2,3,4]);
  assert.ok(p.stress<1e-6);assert.equal(p.eigenvalues.length,1);assert.throws(()=>metricProjection(g,[0,8]),/connected/);
});
test('real export adapters enforce slice membership and display budgets',()=>{
  const idx=indexEvolution(JSON.parse(readFileSync(new URL('../data/example-evolution.json',import.meta.url)))),s=[...idx.states.values()].at(-1),b=branchialView(idx,s);
  assert.ok(b.ids.length<=BUDGET.nodes);assert.ok(b.edges.length<=BUDGET.edges);
  for(const[a,c]of b.edges){assert.equal(idx.states.get(a).step,s.step);assert.equal(idx.states.get(c).step,s.step);}
  const g=graphOf(s),patch=neighborhood(g,g.keys().next().value,99,BUDGET.projection);assert.ok(Number.isFinite(metricProjection(g,patch.ids).stress));
});
