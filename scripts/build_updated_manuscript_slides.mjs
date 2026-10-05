// Nature-inspired PAAD manuscript deck using the finalized R figures.
import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const root=path.resolve(process.argv[2]??'');
const final=path.resolve(process.argv[3]??path.join(root,'deliverables','PAAD_Updated_Manuscript_Presentation_20261005_KR_Explained_v3.pptx'));
const modules=process.env.PAAD_RUNTIME_NODE_MODULES;
const skill=process.env.PAAD_PRESENTATIONS_SKILL;
const python=process.env.PAAD_RUNTIME_PYTHON;
if(!modules||!skill||!python)throw new Error('Required PAAD runtime environment variables are missing.');
process.env.RUNTIME_NODE_MODULES=modules;process.env.RUNTIME_PYTHON=python;
const {Presentation,PresentationFile}=await import(pathToFileURL(path.join(modules,'@oai/artifact-tool/dist/artifact_tool.mjs')).href);
const {finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);

const readTsv=async name=>{
  const text=await fs.readFile(path.join(root,'tables',name),'utf8');
  const rows=text.trim().split(/\r?\n/).map(x=>x.split('\t'));const head=rows.shift();
  return rows.map(row=>Object.fromEntries(head.map((h,i)=>[h,row[i]??''])));
};
const p=x=>{const n=Number(x);return !Number.isFinite(n)?'NE':n<.001?'<0.001':n.toFixed(3)};
const hr=r=>`${Number(r.HR).toFixed(2)} [${Number(r.lower95).toFixed(2)}, ${Number(r.upper95).toFixed(2)}]`;
const lookup=(rows,query)=>rows.find(r=>Object.entries(query).every(([k,v])=>r[k]===v));

const f2=await readTsv('Input_Figure2_global_tests.tsv');
const pg=await readTsv('Input_Preop_platinum_global_comparisons.tsv');
const pc=await readTsv('Input_Preop_platinum_Cox_models.tsv');
const rc=await readTsv('Input_Recurrence_pattern_Cox_models.tsv');
const rlr=(await readTsv('Input_Recurrence_pattern_post_OS_logrank.tsv'))[0];
const surv=await readTsv('Table2_Main_survival_models.tsv');
const klr=lookup(await readTsv('Input_4_Logrank.tsv'),{variable:'KRAS_group',endpoint:'OS'});
const tp=lookup(f2,{gene:'TP53',characteristic:'Differentiation'});
const cd3=lookup(f2,{gene:'CDKN2A',characteristic:'N_stage'});
const cd2=lookup(f2,{gene:'CDKN2A',characteristic:'N_category'});
const kdV=lookup(surv,{variable:'KRAS subtype',comparison_group:'G12D',reference_group:'G12V'});
const kdR=lookup(surv,{variable:'KRAS subtype',comparison_group:'G12D',reference_group:'G12R'});
const lLocal=lookup(rc,{group1:'Local only',group2:'Local + distant'});
const dLocal=lookup(rc,{group1:'Distant only',group2:'Local + distant'});
const repairSets=await readTsv('Input_Preop_platinum_HRD_MMR_set_tests.tsv');
const hrdSet=lookup(repairSets,{gene_set:'HRD',outcome:'Any retained variant'});
const mmrSet=lookup(repairSets,{gene_set:'MMR',outcome:'Any retained variant'});
const distantLR=(await readTsv('Input_Distant_only_pattern_post_OS_logrank.tsv'))[0];
const distantCox=await readTsv('Input_Distant_only_pattern_post_OS_Cox_models.tsv');
const expandedLR=(await readTsv('Input_Expanded_recurrence_post_OS_logrank.tsv'))[0];
const expandedCox=await readTsv('Input_Expanded_recurrence_post_OS_Cox_models.tsv');
const liverLung=lookup(distantCox,{group1:'Liver only',group2:'Lung only'});
const lungOther=lookup(distantCox,{group1:'Lung only',group2:'Other distant'});
const liverOther=lookup(distantCox,{group1:'Liver only',group2:'Other distant'});
const expandedLDLung=lookup(expandedCox,{group1:'Local + distant',group2:'Lung only'});
const signatureAudit=await readTsv('Input_MAF_signature_audit.tsv');
const signatureAuditMap=Object.fromEntries(signatureAudit.map(r=>[r.item,r.value]));
const pooledSignature=await readTsv('Input_Preop_platinum_pooled_signature_sensitivity.tsv');
const individualSignature=await readTsv('Input_Preop_platinum_pooled_signature_by_SBS.tsv');
const requestedSignatureStatus=await readTsv('Input_Signature_requested_metrics_status.tsv');
const assignmentAudit=await readTsv('Input_Signature_target_assignment_audit.tsv');
const platinumPooled=lookup(pooledSignature,{process:'Platinum-associated SBS',metric:'Pooled signature proportion'});
const hrdPooled=lookup(pooledSignature,{process:'HRD-associated SBS',metric:'Pooled signature proportion'});
const significantIndividual=individualSignature.filter(r=>Number(r.p_holm_20_tests)<0.05).length;
const individualProportions=individualSignature.filter(r=>r.metric==='Pooled signature proportion');
const largestIndividual=individualProportions.reduce((a,b)=>Math.abs(Number(a.difference_yes_minus_no))>=Math.abs(Number(b.difference_yes_minus_no))?a:b);
const unavailableMetrics=requestedSignatureStatus.filter(r=>r.status==='Not estimable').length;
const standardTargetTotal=assignmentAudit.reduce((a,r)=>a+Number(r.no_standard_assignment_count)+Number(r.yes_standard_assignment_count),0);
const constrainedTargetTotal=assignmentAudit.reduce((a,r)=>a+Number(r.no_constrained_refit_count)+Number(r.yes_constrained_refit_count),0);

const deck=Presentation.create({slideSize:{width:1280,height:720}});
const family='Pretendard';const ink='#111111', blue='#174A66', pale='#EAF1F4', grey='#666666', orange='#C77732';
const nativeTableSlides=[];const build=path.join(root,'private',`updated_slide_build_${Date.now()}`);await fs.mkdir(build,{recursive:true});
function textbox(slide,value,x,y,w,h,size=24,bold=false,color=ink,align='left'){
  const s=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  s.text=value;s.text.style={typeface:family,fontSize:size,bold,color,alignment:align,verticalAlignment:'top',autoFit:'shrinkText'};return s;
}
function baseSlide(title,section='',notes=''){
  const s=deck.slides.add();s.background.fill='#FFFFFF';
  textbox(s,title,52,25,1170,54,34,true,ink);
  s.shapes.add({geometry:'rect',position:{left:52,top:84,width:1176,height:3},fill:blue,line:{fill:'none',width:0}});
  if(section)textbox(s,section.toUpperCase(),52,94,650,22,13,true,blue);
  textbox(s,String(deck.slides.items.length),1192,684,35,18,13,false,grey,'right');
  s.speakerNotes.textFrame.setText(notes);return s;
}
async function addImage(s,name,x=55,y=115,w=1170,h=500){
  const data=new Uint8Array(await fs.readFile(path.join(root,'figures',name)));
  s.images.add({blob:data,contentType:'image/png',alt:name.replaceAll('_',' '),fit:'contain',position:{left:x,top:y,width:w,height:h}});
}
function inputLine(s,names){textbox(s,names.map((n,i)=>`Input ${i+1}. ${n}`).join('   '),55,654,1110,30,12,false,grey);}
function explanationSlide(title,methods,results,meaning,limits,section){
  const s=baseSlide(title,section,[methods,results,meaning,limits].join('\n\n'));
  const cards=[["분석 방법",methods,55,118,240],["결과",results,650,118,240],["임상적 의미",meaning,55,380,265],["해석 시 주의점",limits,650,380,265]];
  for(const [label,body,x,y,h] of cards){s.shapes.add({geometry:'rect',position:{left:x,top:y,width:575,height:h},fill:'#F7F7F5',line:{fill:'#D8D8D8',width:1}});textbox(s,label,x+22,y+18,530,30,20,true,blue);textbox(s,body,x+22,y+56,530,h-76,20,false,ink);}
  return s;
}
function nativeTable(s,headers,rows,top=150){
  const values=[headers,...rows];const widths=headers.length===5?[320,190,240,180,190]:Array(headers.length).fill(1120/headers.length);
  const table=s.tables.add({rows:values.length,columns:headers.length,left:80,top,width:1120,height:70+rows.length*58,columnWidths:widths,values});
  table.styleOptions={headerRow:true,bandedRows:false};
  table.cells.block({row:0,column:0,rowCount:values.length,columnCount:headers.length}).assign({fill:'#FFFFFF',textStyle:{typeface:family,fontSize:18,color:ink},margins:{left:8,right:8,top:5,bottom:5}});
  table.cells.block({row:0,column:0,rowCount:1,columnCount:headers.length}).assign({fill:'#E6E6E6'});table.cells.block({row:0,column:0,rowCount:1,columnCount:headers.length}).textStyle.bold=true;
  for(let r=1;r<values.length;r++){table.rows[r].height=58;if(r%2===0)table.cells.block({row:r,column:0,rowCount:1,columnCount:headers.length}).assign({fill:'#F5F5F5'});}
  table.rows[0].height=70;nativeTableSlides.push(deck.slides.items.length);return table;
}

let s=deck.slides.add();s.background.fill='#FFFFFF';
textbox(s,'업데이트된 PDAC\n임상-유전체 분석',70,90,780,160,54,true,ink);
s.shapes.add({geometry:'rect',position:{left:70,top:278,width:480,height:5},fill:blue,line:{fill:'none',width:0}});
textbox(s,'분자 아형, 병리, 선행 platinum, 재발양상, 생존',72,310,1080,85,28,false,blue);
textbox(s,'260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx\nKRAS source: KRAS_subtype_MAF · N=1,011 · 전체 병기 포함 · 중복 patient_id 행 제외',72,505,1000,80,22,false,grey);
textbox(s,'2026.10.05',1070,624,145,30,18,true,blue,'right');
s.speakerNotes.textFrame.setText('업데이트 임상자료를 기준으로 전체 논문 흐름을 재구성한 발표자료. 내부 임상/MAF 일치도 비교는 포함하지 않는다.');

s=baseSlide('논문의 전체 흐름과 중심 주장','논문 구성');
const flow=[['1','코호트와 분자적 지형'],['2','Driver와 병리 표현형'],['3','선행 platinum의 임상적 맥락'],['4','재발양상과 재발 후 경과'],['5','유전체 예후 이질성']];
flow.forEach(([n,t],i)=>{const y=135+i*94;s.shapes.add({geometry:'ellipse',position:{left:85,top:y,width:58,height:58},fill:i===3?orange:blue,line:{fill:'none',width:0}});textbox(s,n,85,y+9,58,36,24,true,'#FFFFFF','center');textbox(s,t,170,y+7,720,44,28,i===3,ink);if(i<4)s.shapes.add({geometry:'rect',position:{left:113,top:y+58,width:2,height:36},fill:'#B0B0B0',line:{fill:'none',width:0}});});
textbox(s,'중심 임상적 기여',895,155,310,35,20,true,blue);textbox(s,'Local + distant 재발은 재발 후 경과가 가장 불량한 임상표현형이었다. KRAS는 단순 검출 여부보다 subtype에서 예후 이질성이 더 뚜렷했다.',895,200,300,220,23,false,ink);

s=baseSlide('통계값을 해석하는 기준','통계 해석');
nativeTable(s,['통계값','무엇을 나타내는가','일반적인 해석 기준','읽을 때 주의할 점'],[
  ['p-value','차이 또는 연관성이 없다는 가정과 자료가 얼마나 부합하는지 평가','작을수록 귀무가설과 덜 부합하며 통상 0.05를 기준으로 사용','효과의 크기나 임상적 중요도를 나타내지 않음'],
  ['Adjusted p','여러 검정을 동시에 수행한 영향을 보정한 p-value','Holm은 family-wise error, BH는 false discovery rate를 통제','다중비교가 있으면 raw p보다 adjusted p를 우선 확인'],
  ['95% CI','자료와 모형이 허용하는 효과크기의 불확실성 범위','HR의 CI가 1을 포함하면 차이를 확정하기 어려움','구간의 폭은 추정의 정밀도를 보여줌'],
  ["Cramer's V",'Cramér는 Harald Cramér의 이름, V는 범주형 변수 사이 연관성의 효과크기','0은 연관 없음, 1은 완전한 연관; 약 0.10 small, 0.30 moderate, 0.50 large','표 크기와 임상 맥락에 따라 달라지는 참고 기준'],
  ['Hazard ratio','비교군의 순간 사건 hazard를 기준군과 비교한 비율','1은 차이 없음. 1보다 크면 비교군 hazard 증가, 작으면 감소','누적 사망확률이나 평균 생존기간의 비율이 아님'],
  ['Log-rank p','보정하지 않은 생존곡선 전체가 같은지 비교','낮으면 하나 이상의 곡선이 다르다는 근거','어느 쌍이 다른지는 쌍별 검정이 필요'],
  ['PH test p','Cox HR이 시간에 따라 일정하다는 가정을 점검','0.05 미만이면 proportional hazards 위반 가능성','높은 p가 가정을 입증하는 것은 아님'],
]);

s=baseSlide('분석 코호트 구성','Figure 1a');await addImage(s,'Figure1A_Cohort.png',120,130,1040,470);inputLine(s,['Input_Cohort_flow.tsv']);
s=baseSlide('임상 annotation을 포함한 분자적 지형','Figure 1b');await addImage(s,'Figure1B_Oncoplot.png',25,105,1230,535);inputLine(s,['Input_Figure1B_Oncoplot_gene_counts.tsv','Input_Figure1B_Oncoplot_clinical_tracks.tsv']);
s=baseSlide('KRAS subtype 분포','Figure 1c');await addImage(s,'Figure1C_KRAS_distribution.png',190,115,900,500);inputLine(s,['Input_KRAS_distribution.tsv']);
explanationSlide('Figure 1 결과 해석',
  '중복 ID를 제외한 환자 1명당 1개 열로 구성했다. KRAS를 포함한 유전자 행은 MAF variant class를, 상단 annotation과 KRAS 관련 분석은 엑셀 Main 시트의 KRAS_subtype_MAF를 사용한다.',
  'M0와 M1을 포함한 1,011명이 분석에 포함되었다. Variant count는 기술적 mutation burden으로 제시했으며, 임상정보의 reported TMB와는 별도 측정치로 유지했다.',
  '대규모 코호트에서 분자 아형, 임상병리 특성 및 재발 정보를 같은 환자 단위로 연결해 후속 연관성 분석의 기반을 제시한다.',
  'CNA/LOH, purity, cellularity 및 germline은 추론하지 않았다. “Not detected”는 검사법과 무관하게 확정된 biological wild type을 의미하지 않는다.', 'Figure 1');

s=baseSlide('Driver 변이와 병리 소견의 연관성','Figure 2a');await addImage(s,'Figure2A_Pathology_matrix.png',155,120,970,480);inputLine(s,['Input_2_Pathology_CMH.tsv']);
s=baseSlide('분화도에 따른 TP53 변이 검출률','Figure 2b');await addImage(s,'Figure2_TP53_differentiation.png',180,110,920,500);inputLine(s,['Input_Figure2_global_tests.tsv','Input_Figure2_pairwise_tests.tsv']);
s=baseSlide('림프절 병기에 따른 CDKN2A 변이 검출률','Figure 2c–d');
await addImage(s,'Figure2_CDKN2A_N_stage.png',35,120,585,455);await addImage(s,'Figure2_CDKN2A_N_category.png',660,120,585,455);inputLine(s,['Input_Figure2_counts.tsv','Input_Figure2_pairwise_tests.tsv']);
explanationSlide('Figure 2 결과 해석',
  '변이 검출 여부를 분화도와 N category별로 비교했다. 전체 비교에는 Pearson χ²를 사용했고, 모든 범주 쌍은 Fisher exact test 후 Holm 보정했다.',
  `Cramér는 통계학자 Harald Cramér의 이름이고, V는 범주형 변수의 연관성 크기를 나타내는 효과크기이다. 0은 연관 없음, 1은 완전한 연관을 뜻한다. TP53와 분화도는 p ${p(tp.p)}, V=${Number(tp.cramers_v).toFixed(2)}였다. CDKN2A와 N stage는 p ${p(cd3.p)}, V=${Number(cd3.cramers_v).toFixed(2)}, N category는 p ${p(cd2.p)}, V=${Number(cd2.cramers_v).toFixed(2)}였다.`,
  '관례적 기준에서 V=0.14와 0.11은 모두 small association에 해당한다. 즉 통계적으로 차이는 분명하지만 연관성의 크기는 약하다. TP53 검출률은 WD보다 MD와 PD에서, CDKN2A 검출률은 N0보다 N1에서 높았다.',
  '큰 표본에서는 작은 효과도 낮은 p-value를 보일 수 있다. 따라서 유의한 p-value만으로 강한 생물학적 또는 임상적 연관성을 주장하지 않는다. Cramer\'s V 기준도 표 구조와 맥락에 따라 달라지는 참고값이다.', 'Figure 2');

s=baseSlide('선행 platinum의 임상·유전체 맥락','Figure 3a');await addImage(s,'Figure3_Platinum_context.png',95,112,1090,505);inputLine(s,['Input_Preop_platinum_descriptive.tsv','Input_Preop_platinum_global_comparisons.tsv']);
s=baseSlide('선행 platinum과 OS 및 기록된 재발','Figure 3b–c');await addImage(s,'Figure3_Platinum_OS.png',45,120,650,470);await addImage(s,'Figure3_Platinum_forest.png',710,145,535,400);inputLine(s,['Input_Preop_platinum_OS_logrank.tsv','Input_Preop_platinum_Cox_models.tsv']);
explanationSlide('Figure 3 결과 해석',
  'Yes 301명과 No 706명을 비교하고 Unknown 4명은 제외했다. 전체 KM은 기술적 비교이다. Platinum 노출이 선행치료 여부에 종속되어 있으므로 보정 Cox는 선행치료 환자로 제한했다.',
  `전체 OS log-rank p=0.061은 보정하지 않은 생존곡선 차이가 통상 기준 0.05에 미치지 못했음을 뜻한다. HR은 Yes/No이다. 전체 OS HR ${hr(pc[0])}은 Yes군의 사망 hazard가 17% 높다는 추정이지만 CI가 1을 포함한다. 제한 보정 OS HR ${hr(pc[1])}은 23% 낮다는 추정이지만 역시 CI가 1을 포함한다.`,
  `기록된 재발 HR은 전체 ${hr(pc[2])}, 제한 보정 ${hr(pc[3])}였다. 각각 16% 증가와 33% 감소 방향이지만 두 CI 모두 1을 포함한다. 4개 모형의 Holm 보정 p=${p(pc[0].p_holm_all_models)}이므로 생존 또는 재발 차이를 확정할 통계적 근거가 없다.`,
  '치료효과 추정으로 해석할 수 없다. 전체 재발 모형은 proportional hazards 가정을 위반해 global p=0.003이었고, 선행치료 제한 보정 모형은 p=0.973이었다. Regimen, indication, response 및 dose intensity 정보가 없다.', 'Figure 3');

s=baseSlide('선행 platinum과 HRD/MMR gene-list 변이','Figure 3d');
await addImage(s,'Figure3_Platinum_HRD_MMR.png',70,112,1140,505);
inputLine(s,['Input_Preop_platinum_HRD_MMR_set_tests.tsv','Input_Preop_platinum_HRD_MMR_gene_tests.tsv']);
explanationSlide('Figure 3d 결과 해석',
  'Preop_platinum_exposure Yes/No에서 HRD 18개 및 MMR 4개 gene list의 retained-variant 검출률을 Fisher exact로, 검출 유전자 수를 Wilcoxon으로 비교했다. 4개 set-level 검정은 Holm, 개별 유전자는 list별 BH 보정했다.',
  `HRD-list 검출률은 No ${Number(hrdSet.no_value).toFixed(1)}%, Yes ${Number(hrdSet.yes_value).toFixed(1)}%, OR Yes/No ${Number(hrdSet.effect).toFixed(2)}, Holm p ${p(hrdSet.p_holm)}였다. MMR-list는 No ${Number(mmrSet.no_value).toFixed(1)}%, Yes ${Number(mmrSet.yes_value).toFixed(1)}%, OR ${Number(mmrSet.effect).toFixed(2)}, Holm p ${p(mmrSet.p_holm)}였다. 개별 유전자도 BH 보정 후 유의하지 않았다.`,
  '두 노출군의 baseline HRD/MMR-list retained-variant 분포가 뚜렷하게 다르다는 근거는 없다. 이는 치료군 선택 시 적어도 이 단순 gene-list 지표가 크게 불균형하지 않았음을 보여주는 기술적 결과이다.',
  'Retained variant는 병원성, biallelic loss, LOH, germline, MSI 또는 functional HRD/dMMR 진단이 아니다. 비무작위 치료노출 비교이므로 platinum 감수성이나 치료효과를 뜻하지 않는다.', 'Figure 3');

s=baseSlide('기록된 재발양상의 분포와 특성','Figure 4a–b');await addImage(s,'Figure4_Recurrence_distribution.png',35,120,540,450);await addImage(s,'Figure4_Recurrence_context.png',610,115,630,470);inputLine(s,['Input_Recurrence_pattern_distribution.tsv','Input_Recurrence_pattern_global_comparisons.tsv']);
s=baseSlide('Local + distant 재발에서 가장 불량한 재발 후 생존','Figure 4c–d');await addImage(s,'Figure4_Post_recurrence_OS.png',35,115,650,470);await addImage(s,'Figure4_Post_recurrence_forest.png',700,145,545,405);inputLine(s,['Input_Recurrence_pattern_post_OS_logrank.tsv','Input_Recurrence_pattern_Cox_models.tsv']);
explanationSlide('Figure 4 결과 해석',
  '재발이 기록되고 Recurrence_pattern이 존재하는 659명만 비교했다. 비재발 환자의 공란은 제외했다. 재발 후 시간은 OS months에서 DFS months를 뺀 값으로 정의했다.',
  `전체 log-rank p=${p(rlr.p)}는 세 생존곡선 중 하나 이상이 다르다는 뜻이다. Local only/Local + distant 보정 HR ${hr(lLocal)}은 Local only의 사망 hazard가 42% 낮다는 뜻이며, 반대로 Local + distant는 Local only의 약 ${(1/Number(lLocal.HR)).toFixed(2)}배이다. Holm p ${p(lLocal.p_holm)}였다.`,
  `Distant only/Local + distant HR ${hr(dLocal)}은 Distant only의 사망 hazard가 33% 낮다는 뜻이다. 반대로 Local + distant는 Distant only의 약 ${(1/Number(dLocal.HR)).toFixed(2)}배이며 Holm p ${p(dLocal.p_holm)}였다. 두 CI가 1을 포함하지 않아 Local + distant의 불량한 재발 후 경과를 지지한다.`,
  'Recurrence pattern은 재발 후에 정의되는 post-baseline 변수이다. 진단 시점의 인과적 노출로 사용할 수 없으며, salvage therapy와 추적검사 강도가 재발 후 생존에 영향을 줄 수 있다.', 'Figure 4');

s=baseSlide('Distant only 내부: 원격재발 부위별 생존','Figure 4e–f');
await addImage(s,'Figure4_Distant_only_pattern_KM.png',30,112,645,480);
await addImage(s,'Figure4_Distant_only_pattern_forest.png',690,135,555,430);
inputLine(s,['Input_Distant_only_pattern_post_OS_logrank.tsv','Input_Distant_only_pattern_post_OS_pairwise_logrank.tsv','Input_Distant_only_pattern_post_OS_Cox_models.tsv']);
explanationSlide('Figure 4e–f 결과 해석',
  'Recurrence_pattern=Distant only인 306명을 Distant_pattern으로 Liver only 133명, Lung only 60명, Other distant 113명으로 구분했다. Peritoneum only, Multiple distant 및 그 밖의 원격재발은 Other distant이다. KM/log-rank와 동일 공변량의 Cox를 사용했다.',
  `전체 log-rank p ${p(distantLR.p)}였다. 보정 Cox에서 Liver only/Lung only HR ${hr(liverLung)}, Holm p ${p(liverLung.p_holm)}로 Liver only의 재발 후 사망 hazard가 약 ${Number(liverLung.HR).toFixed(2)}배 높았다. Lung only/Other distant HR ${hr(lungOther)}, Holm p ${p(lungOther.p_holm)}로 Lung only의 hazard가 ${(100*(1-Number(lungOther.HR))).toFixed(1)}% 낮았다.`,
  `Liver only/Other distant HR ${hr(liverOther)}, Holm p ${p(liverOther.p_holm)}로 차이를 지지하지 않았다. 즉 원격재발을 하나로 묶었을 때 가려지는 Lung-only의 상대적으로 양호한 임상경과가 확인되었다.`,
  '재발부위는 post-baseline 정보이다. 전이부담, 구제치료, 영상검사 강도와 같은 미측정 교란이 남아 진단 시점의 예측 또는 인과효과로 해석할 수 없다.', 'Figure 4');

s=baseSlide('재발양상의 5군 확장 분석','Figure 4g–h');
await addImage(s,'Figure4_Expanded_recurrence_KM.png',25,112,650,480);
await addImage(s,'Figure4_Expanded_recurrence_forest.png',690,122,565,465);
inputLine(s,['Input_Expanded_recurrence_post_OS_logrank.tsv','Input_Expanded_recurrence_post_OS_pairwise_logrank.tsv','Input_Expanded_recurrence_post_OS_Cox_models.tsv']);
explanationSlide('Figure 4g–h 결과 해석',
  '기존 Local only와 Local + distant를 유지하고 Distant only를 Liver only, Lung only, Other distant로 분해한 5군을 비교했다. 전체 log-rank 후 10개 쌍별 log-rank 및 Cox 대비를 Holm 보정했다.',
  `전체 log-rank p ${p(expandedLR.p)}였다. Local + distant/Lung only 보정 HR ${hr(expandedLDLung)}, Holm p ${p(expandedLDLung.p_holm)}로 Local + distant의 재발 후 사망 hazard가 Lung only보다 약 ${Number(expandedLDLung.HR).toFixed(2)}배 높았다. 보정 Cox는 공변량 완전사례 652명을 사용했다.`,
  '가장 일관된 분리는 Lung only의 양호한 경과였다. Local + distant는 Lung only보다 불량했지만 Liver only 또는 Other distant와의 직접 Cox 대비는 Holm 보정 후 유의하지 않았다. 따라서 Local + distant가 모든 원격재발보다 항상 불량하다고 주장하지 않는다.',
  '5군 10쌍 비교는 다중검정 부담이 크며 일부 군은 작다. 외부검증 전에는 원격재발 부위를 예후 층화 후보로 제시하는 수준이 적절하다.', 'Figure 4');

s=baseSlide('주요 driver와 보정 OS의 연관성','Figure 5a');await addImage(s,'Figure3A_Gene_OS.png',110,135,1060,420);inputLine(s,['Input_3_Cox_coefficients.tsv','Input_3_PH_tests.tsv']);
s=baseSlide('KRAS subtype별 관찰 및 보정 생존 이질성','Figure 5b–c');await addImage(s,'Figure4A_KRAS_OS.png',25,115,650,480);await addImage(s,'Figure4B_KRAS_pairwise.png',685,120,570,470);inputLine(s,['Input_4_Logrank.tsv','Input_4_Logrank_pairs.tsv','Input_4_Cox_pairs.tsv']);
explanationSlide('Figure 5 결과 해석',
  'KM 및 log-rank는 보정하지 않은 생존곡선을 비교한다. Cox contrast는 연령, 성별, T/N/M, 분화도, 선행치료, LVI/PNI와 절제연을 보정했다. KRAS 10개 쌍은 Holm 보정했다.',
  `KRAS 전체 log-rank p=${p(klr.p)}는 하나 이상의 subtype 생존곡선이 다르다는 뜻이다. G12D/G12V 보정 HR ${hr(kdV)}은 G12D의 순간 사망 hazard가 G12V보다 39% 높다는 뜻이며 Holm p ${p(kdV.adjusted_p)}였다. CI가 1을 포함하지 않는다.`,
  `G12D/G12R HR ${hr(kdR)}은 G12D의 순간 사망 hazard가 G12R보다 67% 높다는 뜻이며 Holm p ${p(kdR.adjusted_p)}였다. KRAS를 단순 검출/미검출로 구분하면 희석되는 예후 이질성이 subtype 수준에서 확인되었다.`,
  'HR은 누적 사망확률이나 평균 생존기간의 비율이 아니다. KM과 forest plot의 p-value는 서로 다른 가설과 보정 방법을 사용한다. 전체 log-rank가 유의하더라도 모든 subtype 쌍이 서로 다르다는 뜻은 아니다.', 'Figure 5');

s=baseSlide('KM과 Cox의 p-value가 다른 이유','통계 해석');
nativeTable(s,['분석','검정 질문','공변량 보정','다중검정 보정','해석'],[
  ['전체 log-rank','하나 이상의 곡선이 다른가?','없음','Omnibus test 1회','차이가 나는 특정 쌍을 알려주지 않음'],
  ['쌍별 log-rank','어느 두 곡선이 다른가?','없음','전체 쌍 내 Holm','보정하지 않은 쌍별 근거'],
  ['Cox forest','비교군/기준군의 보정 HR은?','임상 공변량','Family별 BH 또는 Holm','보정 후 방향과 효과크기'],
  ['PH test','시간에 따라 일정한 HR이 타당한가?','모형별','진단 검정','높은 p가 PH를 입증하지는 않음'],
]);
textbox(s,'KM 곡선이 시각적으로 분리되어도 공변량 또는 다중검정 보정 후 유의하지 않을 수 있다. 반대로 교란을 보정한 Cox contrast에서 유의성이 나타날 수도 있다.',80,545,1120,80,22,false,blue);

s=baseSlide('Table 1. 코호트 특성','주요 표');await addImage(s,'Table1_Cohort_characteristics.png',35,105,1210,535);inputLine(s,['Table1_Cohort_characteristics.tsv']);
s=baseSlide('Table 2. 주요 생존모형과 쌍별 비교','주요 표');await addImage(s,'Table2_Main_survival_models.png',20,105,1240,535);inputLine(s,['Table2_Main_survival_models.tsv']);

s=baseSlide('HRD 및 MMR gene-list 변이 양상','보충 그림');await addImage(s,'Supplementary_Figure1_HRD.png',25,115,610,470);await addImage(s,'Supplementary_Figure2_MMR.png',645,115,610,470);textbox(s,'Variant 분포만을 표시한다. Functional HRD/dMMR, LOH 또는 germline 진단을 의미하지 않는다.',60,610,1120,36,20,true,orange);inputLine(s,['Input_Supplementary_Figure1_HRD_gene_counts.tsv','Input_Supplementary_Figure2_MMR_gene_counts.tsv']);
s=baseSlide('Reported TMB의 탐색적 분석','보충 그림');await addImage(s,'Supplementary_Figure4_Reported_TMB.png',95,120,1090,455);textbox(s,'관측 836명, 결측 175명, 단위 및 검사법 미확인, 외부 TMB-high cutoff 미적용',95,595,1090,42,21,true,orange);inputLine(s,['Input_7_Reported_TMB.tsv','Input_7_TMB_associations.tsv']);

s=baseSlide('Mutational signature 분석 가능성 점검','보충 그림');
await addImage(s,'Supplementary_Figure_Signature_QC.png',55,120,1170,455);
textbox(s,`환자별 분석 결과: ${signatureAuditMap['Cohort samples with >= 10 context-counted SNVs']}/1,011명만 SNV≥10, 최고 cosine 0.899, cosine≥0.90인 환자 0명`,70,590,1140,42,20,true,orange);
inputLine(s,['Input_MAF_signature_audit.tsv','Input_Signature_sample_quality_summary.tsv','Input_Preop_platinum_signature_tests.tsv']);
explanationSlide('Signature 품질 결과 해석',
  'GRCh37 SNV로 SBS96을 만들고 COSMIC v3.6 exome reference에 적합했다. 환자별 기준은 SNV≥10, cosine≥0.90, process count≥5, proportion≥0.20, bootstrap stability≥0.80이었다.',
  `Context-counted SNV≥10은 ${signatureAuditMap['Cohort samples with >= 10 context-counted SNVs']}명이었으나 reconstruction cosine 최고값이 0.899여서 품질기준 통과 환자는 0명이었다. 따라서 SBS31/35, SBS3, MMR SBS의 환자별 positive rate·burden·proportion은 모두 NE이다.`,
  '이것은 signature가 없다는 결론이 아니라, 현재 targeted-panel retained-variant MAF로 환자별 signature를 신뢰성 있게 판정할 수 없다는 분석 가능성 결과이다. 기준을 낮춰 양성률을 만드는 것은 피했다.',
  'Synonymous 변이가 0개이고 환자당 SNV 중앙값은 7개이며 panel BED/callable territory가 없다. 확정적 signature 분석에는 원시 variant call과 panel BED 또는 WES/WGS가 필요하다.', '보충 그림');

s=baseSlide('요청된 환자별 signature 지표의 산출 가능성','보충 그림');
await addImage(s,'Supplementary_Figure_Signature_Requested_Metrics.png',95,125,1090,450);
textbox(s,`Positive rate, patient-level burden, proportion: ${unavailableMetrics}/9개 지표 NE`,95,595,1090,40,21,true,orange);
inputLine(s,['Input_Signature_requested_metrics_status.tsv','Input_Preop_platinum_signature_tests.tsv']);
explanationSlide('환자별 signature 지표 결과 해석',
  'Platinum, HRD, MMR process 각각에 대해 positive rate, attributed retained-SNV burden, proportion을 명시적으로 평가했다. 두 노출군 모두 품질평가 가능 환자가 있어야 비교를 시행하도록 사전 정의했다.',
  `총 9개 요청 지표 중 ${unavailableMetrics}개가 NE였다. No와 Yes 양쪽에서 SNV 수와 cosine 기준을 동시에 충족한 환자가 0명이므로 p-value와 효과크기를 계산하지 않았다.`,
  'NE는 signature가 없다는 뜻이 아니라 현재 targeted-panel retained-variant MAF로 환자 수준의 존재 여부와 크기를 신뢰성 있게 추정할 수 없다는 뜻이다.',
  '결과를 본 뒤 cosine 또는 변이수 기준을 낮추면 post-hoc thresholding이 된다. 따라서 기준을 낮춰 양성률을 만드는 분석은 추가하지 않았다.', '보충 그림');

s=baseSlide('집단 합산 SBS 민감도 분석','보충 그림');
await addImage(s,'Supplementary_Figure_Platinum_HRD_MMR_signatures.png',80,110,1120,510);
inputLine(s,['Input_Preop_platinum_pooled_signature_sensitivity.tsv','Input_Signature_analysis_definitions.tsv']);
explanationSlide('집단 합산 SBS 결과 해석',
  '환자별 추론이 불가능해 Yes/No군 SBS96을 합산하고 target SBS를 포함한 constrained refit을 시행했다. 환자 bootstrap 1,000회로 CI, label permutation 1,000회로 p-value를 구하고 6개 검정을 Holm 보정했다.',
  `Platinum SBS proportion은 No ${Number(platinumPooled.no_value).toFixed(3)}, Yes ${Number(platinumPooled.yes_value).toFixed(3)}, Holm p ${p(platinumPooled.p_holm_6_tests)}였다. HRD SBS proportion은 No ${Number(hrdPooled.no_value).toFixed(3)}, Yes ${Number(hrdPooled.yes_value).toFixed(3)}, Holm p ${p(hrdPooled.p_holm_6_tests)}였다. MMR process도 차이가 없었다.`,
  '어떤 process도 platinum Yes/No 차이를 지지하지 않았다. 그러나 이는 patient-level positive rate가 아니라 pooled spectrum의 방법론적 민감도 분석이다.',
  'Standard pooled SigProfilerAssignment에서는 target process count가 모두 0이었고 그림 값은 constrained refit에서만 나타났다. 검체가 모두 치료 전이므로 SBS31/35를 현재 platinum 치료에 의해 유발된 신호로 해석할 수 없다.', '보충 그림');

s=baseSlide('개별 COSMIC SBS pooled 민감도 분석','보충 그림');
await addImage(s,'Supplementary_Figure_Individual_SBS_Pooled_Sensitivity.png',45,105,1190,520);
inputLine(s,['Input_Preop_platinum_pooled_signature_by_SBS.tsv','Input_Signature_analysis_definitions.tsv']);
explanationSlide('개별 SBS 결과 해석',
  'SBS31·35, SBS3, SBS6·14·15·20·21·26·44 각각의 pooled attributed SNVs/patient와 proportion을 Yes-No 차이로 계산했다. Patient bootstrap CI와 label permutation p를 사용하고 20개 검정을 Holm 보정했다.',
  `Holm 보정 p<0.05인 개별 signature-metric 결과는 ${significantIndividual}개였다. Absolute proportion 차이가 가장 큰 것은 ${largestIndividual.signature}로 Yes-No ${Number(largestIndividual.difference_yes_minus_no).toFixed(3)}, raw p ${p(largestIndividual.permutation_p)}, Holm p ${p(largestIndividual.p_holm_20_tests)}였다.`,
  '개별 SBS까지 분리해도 노출군 차이를 지지하는 견고한 결과가 없다는 점이 핵심이다. 통합 process 결과에 의해 특정 SBS의 상반된 방향이 가려지는지도 함께 확인했다.',
  '이 결과는 pooled constrained refit 민감도 분석이다. 환자별 positive rate가 아니며 treatment-induced signature의 증거도 아니다.', '보충 그림');

s=baseSlide('Standard assignment와 constrained refit 감사','보충 그림');
await addImage(s,'Supplementary_Figure_Signature_Assignment_Audit.png',135,110,1010,505);
inputLine(s,['Input_Signature_target_assignment_audit.tsv','Input_Preop_platinum_pooled_signature_by_SBS.tsv']);
explanationSlide('Signature assignment 의존성 해석',
  '동일한 pooled SBS96에서 자료 주도 standard SigProfilerAssignment와 사전지정 target SBS를 강제로 후보에 포함한 constrained refit의 attributed count를 비교했다.',
  `Target SBS의 standard assigned count 합은 ${standardTargetTotal.toFixed(1)}, constrained-refit 합은 ${constrainedTargetTotal.toFixed(1)}였다. Target 신호는 standard 선택에서는 채택되지 않았고 forced inclusion에서만 분배되었다.`,
  '이 차이는 constrained 결과의 불확실성과 model dependence를 보여준다. 따라서 0보다 큰 constrained coefficient 자체를 signature 존재 증거로 해석하지 않는다.',
  '현재 자료에서 방어 가능한 결론은 patient-level 분석은 불가능하며 pooled sensitivity에서도 노출군 차이를 지지하지 않는다는 것이다. 원시 call과 BED 또는 WES/WGS가 필요하다.', '보충 그림');

s=baseSlide('임상적 novelty는 규모와 통합된 경과 분석에 있다','고찰');
const claims=[['대규모 통합 코호트','1,011명에서 sequencing, 수술 병리, 치료 맥락, 재발 및 생존을 함께 분석했다.'],['Subtype 수준의 예후','KRAS G12D는 G12V 및 G12R보다 높은 보정 사망 hazard와 연관되었다.'],['재발부위의 이질성','Local + distant의 불량한 경과와 함께 Distant only 내부의 Lung-only 양호 경과를 확인했다.'],['음성·불가 결과의 투명성','Platinum 관련 결과의 비유의성과 patient-level SBS 분석 불가능성을 품질기준과 함께 제시했다.']];
claims.forEach(([a,b],i)=>{const y=125+i*125;s.shapes.add({geometry:'rect',position:{left:70,top:y,width:1140,height:100},fill:i===2?'#F5E8DC':pale,line:{fill:'#D0DADD',width:1}});textbox(s,a,95,y+18,300,35,22,true,i===2?orange:blue);textbox(s,b,400,y+17,770,62,23,false,ink);});

s=baseSlide('연구의 제한점과 해석 범위','고찰');
const limits=['후향적 단일기관 연구이며 수술 환자 선택편향이 있음','치료 적응증, regimen, dose intensity 및 response 정보가 없음','Reported TMB의 단위, 검사법, 버전 및 결측기전이 확인되지 않음','재발부위 분석은 post-baseline이며 구제치료·전이부담의 잔여교란이 있음','Targeted-panel MAF의 낮은 SNV 수와 BED 부재로 환자별 SBS 추론 불가','CNV/LOH, germline 및 functional HRD/dMMR 미분석; 외부검증 필요'];
limits.forEach((t,i)=>{const col=i%2,row=Math.floor(i/2);s.shapes.add({geometry:'rect',position:{left:65+col*585,top:135+row*160,width:550,height:128},fill:'#F7F7F5',line:{fill:'#DADADA',width:1}});textbox(s,String(i+1).padStart(2,'0'),85+col*585,154+row*160,58,38,24,true,orange);textbox(s,t,150+col*585,150+row*160,435,75,23,false,ink);});

s=baseSlide('Summary. 임상적으로 해석 가능한 PDAC 경과 지도','요약','대규모 코호트, 분자-병리 연관성, KRAS subtype 예후, 재발부위 이질성, platinum 결과와 SBS 분석가능성.');
textbox(s,'자료에서 확인된 결과',65,125,500,35,23,true,blue);
textbox(s,'1. TP53 검출률은 분화도와 연관되었다\n2. CDKN2A 검출률은 림프절 병기와 연관되었다\n3. KRAS G12D는 subtype 수준에서 불량한 예후를 보였다\n4. Distant only 내부에서 Lung only가 상대적으로 양호했다',65,175,550,240,28,false,ink);
textbox(s,'임상적으로 중요한 부분',690,125,500,35,23,true,blue);
textbox(s,'가장 중요한 novelty는 1,011명에서 분자 아형과 수술 병리를 재발부위 및 재발 후 경과까지 연결해, 원격재발 내부의 임상적 이질성까지 정량화했다는 점이다.',690,175,500,170,30,true,ink);
s.shapes.add({geometry:'rect',position:{left:65,top:500,width:1125,height:100},fill:'#F5E8DC',line:{fill:'none',width:0}});
textbox(s,'중심 주장',90,522,180,34,22,true,orange);textbox(s,'대규모 PDAC 임상-유전체 코호트에서 subtype별 예후 차이와 재발부위별 경과 이질성을 확인했다. 관찰자료이므로 치료 인과효과와 mutational signature 기전은 주장하지 않는다.',265,515,890,65,26,true,ink);

await fs.writeFile(path.join(build,'slide_manifest.json'),JSON.stringify({slides:deck.slides.items.length,nativeTableSlides},null,2));
const candidate=path.join(build,'candidate.pptx');await(await PresentationFile.exportPptx(deck)).save(candidate);
await fs.mkdir(path.dirname(final),{recursive:true});
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:final,pythonExecutable:python,
  integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',...nativeTableSlides.flatMap(n=>['--require-native-table-slide',String(n)])],
  requiredNativeTableOwnerSlides:nativeTableSlides,fontPolicy:{basis:'user_request',families:['Pretendard']},verifyArtifactToolImport:true,
  receiptPath:path.join(build,'validation.json')});
for(let i=0;i<deck.slides.items.length;i++){
  const blob=await deck.export({slide:deck.slides.items[i],format:'png',scale:1.25});
  await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
console.log(final);console.log(build);
