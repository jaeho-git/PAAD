// Build an editable Nature-inspired scientific deck from content.json and the
// R-generated evidence figures. Tables and prose remain native slide objects.
// Required environment: PAAD_RUNTIME_NODE_MODULES, PAAD_PRESENTATIONS_SKILL,
// PAAD_RUNTIME_PYTHON. Use the workspace dependency runtime, not global installs.
import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
const root = path.resolve(process.argv[2] ?? '');
const modules = process.env.PAAD_RUNTIME_NODE_MODULES;
const skill = process.env.PAAD_PRESENTATIONS_SKILL;
const python = process.env.PAAD_RUNTIME_PYTHON;
if (!modules || !skill || !python) throw new Error('Set PAAD_RUNTIME_NODE_MODULES, PAAD_PRESENTATIONS_SKILL and PAAD_RUNTIME_PYTHON.');
process.env.RUNTIME_NODE_MODULES=modules;
process.env.RUNTIME_PYTHON=python;
const { Presentation, PresentationFile } = await import(pathToFileURL(path.join(modules,'@oai/artifact-tool/dist/artifact_tool.mjs')).href);
const { finalizePresentation } = await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const build=path.join(root,'private',`slide_build_${Date.now()}`);await fs.mkdir(build,{recursive:true});
const content=JSON.parse(await fs.readFile(path.join(root,'private','document_build','content.json'),'utf8'));
const presentation=Presentation.create({slideSize:{width:1280,height:720}});
const family='Pretendard';
const nativeTables=[];const slideInfo=[];
function text(slide,value,x,y,w,h,size=24,bold=false,color='#111111') {
  const s=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  s.text=value;s.text.style={typeface:family,fontSize:size,bold,color,autoFit:'none'};
  return s;
}
function slide(title,notes='') {
  const s=presentation.slides.add();s.background.fill='#FFFFFF';
  text(s,title,52,28,1176,65,37,true);
  text(s,String(presentation.slides.items.length),1193,683,40,22,14,false,'#666666');
  if(content.review_note)text(s,'검토본 · KRAS 출처 확정 전',52,688,750,20,12,false,'#777777');
  s.speakerNotes.textFrame.setText(notes);
  slideInfo.push({number:presentation.slides.items.length,title});
  return s;
}
function paragraphs(s,blocks,top=116) {
  let y=top;
  for(const [label,body] of blocks){
    text(s,label,56,y,1168,30,22,true);y+=34;
    // Full-width Korean explanatory prose, separated from image slides.
    const estimatedLines=Math.max(2,Math.ceil(body.length/73));
    const h=estimatedLines*30+8;
    text(s,body,56,y,1168,h,23);y+=h+16;
  }
}
async function evidence(s,filename,x=48,y=100,w=1184,h=540) {
  s.images.add({blob:new Uint8Array(await fs.readFile(path.join(root,'figures',filename))),contentType:'image/png',alt:filename.replaceAll('_',' '),fit:'contain',position:{left:x,top:y,width:w,height:h}});
}
function cohortDiagram(s,spec) {
  // Native editable shapes reuse the same aggregate geometry as the R figure.
  const dx=140,dy=116;
  for(const e of spec.edges){
    const anchor=(x,y)=>s.shapes.add({geometry:'rect',position:{left:dx+x-.1,top:dy+y-.1,width:.2,height:.2},fill:'none',line:{fill:'none',width:0}});
    const vertical=e.x===e.xend;
    s.shapes.connect(anchor(e.x,e.y),anchor(e.xend,e.yend),{
      kind:'straight',fromSide:vertical?'bottom':'right',toSide:vertical?'top':'left',
      // Artifact Tool uses the OOXML tail end for the destination arrowhead.
      line:{fill:'#555555',width:1.5},tail:{type:e.arrow?'triangle':'none',width:'sm',length:'sm'}});
  }
  for(const n of spec.nodes){
    const box=s.shapes.add({geometry:'rect',name:`Cohort ${n.id}`,position:{left:dx+n.x,top:dy+n.y,width:n.w,height:n.h},fill:n.fill,line:{fill:'#555555',width:1.5}});
    box.text=n.label;
    box.text.style={typeface:family,fontSize:n.font_size,bold:n.id==='analysis',color:'#111111',alignment:'center',verticalAlignment:'middle',autoFit:'none'};
  }
}
function inputs(s,names){
  const strings=names.map((x,i)=>`Input ${i+1}. ${x}`);
  text(s,strings.join('\n'),52,642,1125,44,12,false,'#555555');
}
function nativeTable(s,headers,rows,{top=120,height=470,font=20}={}) {
  const values=[headers,...rows];const n=headers.length;
  const widths=n===6?[335,...Array(5).fill(165)]:n===5?(headers[0]==='Comparison group'?[245,245,285,120,265]:headers[0]==='Model'?[400,70,230,230,230]:headers[0]==='Exposure / contrast'?[440,110,270,170,170]:[175,370,155,230,230]):n===4?(headers[0]==='Gene and outcome'?[300,350,120,390]:headers[0].startsWith('Endpoint')?[460,100,400,200]:[420,340,200,200]):headers[0]==='Study'?[270,330,560]:[575,350,235];
  const headerHeight=n>=5||headers.some(h=>h.length>22)?78:55;
  const t=s.tables.add({rows:values.length,columns:n,left:60,top,width:1160,height,columnWidths:widths,values});
  t.styleOptions={headerRow:true,bandedRows:false};
  const noLine={width:0.2,color:'#FFFFFF'};
  const borderConfig={top:noLine,bottom:noLine,left:noLine,right:noLine};
  t.cells.block({row:0,column:0,rowCount:values.length,columnCount:n}).assign({fill:'#FFFFFF',textStyle:{typeface:family,fontSize:font,color:'#111111'},margins:{left:8,right:8,top:2,bottom:2}});
  t.cells.block({row:0,column:0,rowCount:1,columnCount:n}).textStyle.bold=true;
  // Explicit row heights prevent the default cell height from expanding a
  // 14-row table into the footnotes. Native cells stay editable in PowerPoint.
  t.rows[0].height=headerHeight;
  for(let i=1;i<values.length;i++)t.rows[i].height=(height-headerHeight)/(values.length-1);
  for(let i=0;i<values.length;i++)for(let j=0;j<n;j++){
    t.cells.block({row:i,column:j,rowCount:1,columnCount:1}).assign({borders:borderConfig});
  }
  for(const y of [top,top+headerHeight,top+height])s.shapes.add({geometry:'rect',position:{left:60,top:y,width:1160,height:1},fill:'#111111',line:{fill:'none',width:0}});
  nativeTables.push(presentation.slides.items.length);
}
function tableSlides(t,limit=11) {
  const longRows=t.headers[0].startsWith('Endpoint')||t.headers[0]==='Gene and outcome';
  const rowHeight=longRows?64:t.headers.length>=5?48:44;
  limit=Math.min(limit,Math.floor((450-78)/rowHeight));
  // Balance continuation pages so the final page is never a lone row.
  const chunkSize=Math.ceil(t.rows.length/Math.ceil(t.rows.length/limit));
  const chunks=[];for(let i=0;i<t.rows.length;i+=chunkSize)chunks.push(t.rows.slice(i,i+chunkSize));
  for(let i=0;i<chunks.length;i++){
    const s=slide(`${t.id}. ${t.title}${chunks.length>1?` (${i+1}/${chunks.length})`:''}`,t.note+'\n'+t.inputs.map((x,j)=>`Input ${j+1}. ${x}`).join('\n'));
    nativeTable(s,t.headers,chunks[i],{height:78+chunks[i].length*rowHeight,font:t.headers.length>=5?17:18});
    text(s,t.note,60,588,1160,52,15,false,'#444444');inputs(s,t.inputs.slice(0,3));
  }
}
let s=slide(content.title);
text(s,content.subtitle,56,182,1150,90,43,true);
text(s,'주요 유전자와 병리 소견의 연관성\nKRAS 아형의 예후와 유전체 정보의 추가 가치',56,336,1120,120,31);
if(content.review_note)text(s,content.review_note,56,475,1120,55,18,false,'#555555');
text(s,`Main Figures 1–${content.main_figure_count??4}   Main Tables 1–2\nSupplementary Figures 1–${content.figures.filter(f=>f.id.startsWith('Supplementary')).length}   Supplementary Tables 1–${content.tables.filter(t=>t.id.startsWith('Supplementary')).length}`,56,551,1100,65,23,false,'#555555');
s=slide('연구 질문과 분석 대상',content.methods.join('\n\n'));
paragraphs(s,[['연구 질문','주요 유전자와 KRAS 아형은 병리 소견 및 생존과 어떻게 연관되며, 이미 알려진 임상병리정보에 어떤 예후 정보를 추가하는가?'],['분석 대상',content.methods[0]],['핵심 해석','예후 연관성과 추가 예측 가치를 구분한다. 통계적 유의성이 임상적 유용성 또는 치료 효과를 뜻하지는 않는다.']]);
s=slide('비교군과 기준군에 따른 HR 해석');
paragraphs(s,[['분자와 분모','HR = 비교군의 사건 순간위험률 / 기준군의 사건 순간위험률. 모든 표와 그림에서 Reference를 기준군으로 표시한다. OS의 사건은 전체 사망이다.'],['HR 1.43의 의미','A군을 비교군, B군을 기준군으로 정했다면 A군의 사망 순간위험률이 B군보다 43% 높게 추정된다. HR이 0.70이면 A군에서 30% 낮게 추정된다.'],['해석의 범위','누적 사망확률이나 암 발병확률이 그 비율만큼 변한다는 뜻은 아니다. 방향은 점추정이며 95% CI와 다중비교 보정 p를 함께 확인한다.']]);
s=slide('분석 시점과 통계 모형');
paragraphs(s,[['측정 및 생존 기준',content.methods[1]],['병리 보정','연령 스플라인, 성별, T/N/M category, 선행치료, LVI/PNI 및 절제연을 보정한다. 분화도의 비례위험 위반을 고려해 분화도별 기저위험을 허용한다.'],['검증','같은 환자에서 임상병리 단독과 유전체 추가 모형을 비교하고 환자 bootstrap 500회로 낙관성을 보정한다. OOB 보정도와 성능의 변동성을 함께 제시한다.']]);
for(const f of content.figures.filter(f=>f.id.startsWith('Figure '))) {
  if(f.id==='Figure 1'){
    for(const panel of f.panels){
      s=slide(`${panel.id}. ${panel.title}`,panel.caption+'\n'+f.method+'\n'+f.limitation);
      if(panel.diagram)cohortDiagram(s,panel.diagram);
      else await evidence(s,panel.image,40,100,1200,535);
      inputs(s,panel.inputs);
    }
  }else{
    for(let i=0;i<f.assets.length;i++){
      s=slide(`${f.id}${String.fromCharCode(97+i)}. ${f.panels?.[i]?.title??f.title}`,f.method+'\n'+f.finding+'\n'+f.limitation);
      await evidence(s,f.assets[i]);inputs(s,f.panels?.[i]?.inputs??f.inputs.slice(0,3));
    }
  }
  s=slide(`${f.id} 해석`,[f.method,f.finding,f.interpretation,f.limitation].join('\n\n'));
  paragraphs(s,[['관찰된 결과',f.finding],['임상적 의미',f.interpretation]]);
  s=slide(`${f.id} 분석 방법과 해석 범위`,[f.method,f.limitation].join('\n\n'));
  paragraphs(s,[['분석 방법',f.method],['해석의 한계',f.limitation]]);
  if(f.id==='Figure 1')tableSlides(content.tables[0],14);
  if(f.id==='Figure 4')tableSlides(content.tables[1],13);
}
for(const f of content.figures.filter(f=>f.id.startsWith('Supplementary'))){
  s=slide(`${f.id}. ${f.title}`,[f.method,f.finding,f.interpretation,f.limitation].join('\n\n'));await evidence(s,f.image);inputs(s,f.inputs.slice(0,3));
  s=slide(`${f.id} 설명과 해석`,[f.method,f.finding,f.interpretation,f.limitation].join('\n\n'));
  paragraphs(s,[['관찰된 결과',f.finding],['해석',f.interpretation],['한계',f.limitation]]);
}
// Supplementary tables remain true formatted result tables, not filenames.
for(const t of content.tables.slice(2))tableSlides(t,14);
s=slide('기존 PDAC 연구와의 비교',content.references.slice(0,4).map(r=>r.label+' '+r.url).join('\n'));
nativeTable(s,['Study','Cohort','Interpretation'],[
  ['McIntyre 2020','283명 수술 코호트','주요 유전자와 예후'],
  ['Campbell 2025','508명 수술 코호트','병리·driver gene 수와 OS'],
  ['Varghese 2025','2,336명 전체 코호트','더 큰 기존 PDAC 연구'],
  ['본 분석','1,011명\n사망 631건','상세 병리와 아형별 예후 연결'],
],{height:265,font:22});
text(s,'연구별 병기, 검사 범위, 생존 시작점이 다르므로 절대 생존기간을 단순 비교하지 않는다. 최대 규모 또는 KRAS 아형 예후 차이의 최초 발견으로 주장하지 않는다.',60,443,1160,115,25);
for(const part of content.discussion.slice(0,-1)){
  s=slide(part.title,part.text+'\n'+content.references.map(r=>r.label+' '+r.url).join('\n'));
  text(s,part.text,60,150,1160,400,29);
}
for(let start=0;start<content.references.length;start+=6){
  s=slide('참고문헌',content.references.map(r=>r.url).join('\n'));
  content.references.slice(start,start+6).forEach((r,i)=>text(s,`[${start+i+1}] ${r.label}\n${r.url}`,60,115+i*82,1160,74,20));
}
s=slide('Summary. 임상적 의미와 연구의 기여',content.discussion.at(-1).text);
paragraphs(s,[['최종 결과와 임상적 의미',content.discussion.at(-1).text],['해석 원칙','비교군과 기준군, 통계적 불확실성을 명확히 표시한다. 유전체의 예후 연관성과 실제 추가 예측 가치를 구분하고 개선이 작거나 비유의한 결과도 보고한다.']]);
await fs.writeFile(path.join(build,'slide_manifest.json'),JSON.stringify(slideInfo,null,2));
const candidate=path.join(build,'candidate.pptx');await(await PresentationFile.exportPptx(presentation)).save(candidate);
const final=process.argv[3]?path.resolve(process.argv[3]):path.join(root,'deliverables',`PAAD_Manuscript_Presentation_${new Date().toISOString().slice(0,10).replaceAll('-','')}.pptx`);
await fs.mkdir(path.dirname(final),{recursive:true});
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:final,pythonExecutable:python,
  integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',...nativeTables.flatMap(n=>['--require-native-table-slide',String(n)])],
  requiredNativeTableOwnerSlides:nativeTables,fontPolicy:{basis:'user_request',families:['Pretendard']},verifyArtifactToolImport:true,
  receiptPath:path.join(build,'validation.json')});
// Render the authored slide content; final PPTX receives separate import/render QA.
for(let i=0;i<presentation.slides.items.length;i++){
  const blob=await presentation.export({slide:presentation.slides.items[i],format:'png',scale:1.5});
  await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
console.log(final);
console.log(`QA renders: ${build}`);
