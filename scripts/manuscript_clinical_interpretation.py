"""Add clinical-value results to the shared report/deck narrative.

Reads aggregate tables only. Findings are generated from actual estimates,
including negative or small incremental results; no success threshold is tuned.
"""
import csv
import math


def extend_content(content, root):
    def read(name):
        with (root / "tables" / name).open(encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    def number(v):
        try: return float(v)
        except (ValueError, TypeError): return float("nan")
    def fmt(v, digits=2): return f"{number(v):.{digits}f}" if math.isfinite(number(v)) else "NE"
    def pv(v): return "<0.001" if number(v) < .001 else fmt(v, 3)
    def get(rows, **keys): return next(r for r in rows if all(r[k] == v for k, v in keys.items()))
    def direction(hr): return "높게" if number(hr) >= 1 else "낮게"
    def hr_sentence(r, a=None, b=None):
        a = a or r.get("group1", r.get("comparison_group"))
        b = b or r.get("group2", r.get("reference_group"))
        adjusted = r.get("p_holm", r.get("q_BH", "NA"))
        evidence = "다중비교 보정 후에도 차이를 지지한다" if number(adjusted) < .05 else "다중비교 보정 후에는 차이를 확정하지 못했다"
        return (f"{b}군을 기준으로 {a}군의 사망 순간위험률은 {abs(number(r['HR'])-1)*100:.1f}% {direction(r['HR'])} 추정되었다. "
                f"HR {fmt(r['HR'])}, 95% CI {fmt(r['lower95'])}–{fmt(r['upper95'])}, 보정 p {pv(adjusted)}로 {evidence}.")
    def tab(identifier, title, headers, rows, inputs, note):
        return dict(id=identifier, title=title, headers=headers, rows=rows, inputs=inputs, note=note)
    figs = {f["id"]: f for f in content["figures"]}
    pairs = [r for r in read("Supplementary_Table4_Cox_pairs.tsv") if r["family"] == "Extended genomic" and r["exposure"] == "KRAS_group" and r["endpoint"] == "OS"]
    performance = read("Supplementary_Table9_Prediction_performance.tsv")
    differences = read("Supplementary_Table9_Prediction_differences.tsv")
    path = read("Supplementary_Table8_Adjusted_pathology.tsv")
    survival = read("Supplementary_Table10_Adjusted_survival.tsv")
    absolute = read("Supplementary_Table10_Survival_differences.tsv")
    interaction = read("Supplementary_Table11_Interaction.tsv")[0]
    models = read("Input_Prediction_models.tsv")
    tv = read("Supplementary_Table12_Time_varying_HR.tsv")
    missing = read("Supplementary_Table13_TMB_missingness.tsv")
    labels = {"Clinical":"임상병리정보", "Clinical_KRAS":"임상병리 + KRAS 아형", "Clinical_driver_count":"임상병리 + driver count", "Clinical_genes":"임상병리 + KRAS 아형 및 3개 유전자"}
    kr_delta = get(differences, model="Clinical_KRAS", metric="C_index")
    gene_delta = get(differences, model="Clinical_genes", metric="C_index")
    c_base = get(performance, model="Clinical", metric="C_index")
    c_kras = get(performance, model="Clinical_KRAS", metric="C_index")
    c_count = get(performance, model="Clinical_driver_count", metric="C_index")
    c_genes = get(performance, model="Clinical_genes", metric="C_index")
    primary_ph = read("Supplementary_Table3_PH_tests.tsv")
    global_ph = get(primary_ph, family="Extended genomic", exposure="KRAS_group", term="GLOBAL")
    timing = ("임상 검사 및 분자정보의 측정 시점은 연구자 확인에 따라 진단 후 선행치료 전으로 정의하였다. "
              "수술 병리와 절제연, 실제 치료 수진 정보는 해당 진료 과정에서 얻는 정보로 구분하였다. OS 시작점은 원자료 Dictionary에 정의된 수술 시점을 유지하였다. "
              "따라서 병리정보를 포함한 모형은 수술 환자에서 수술 후 예후를 평가하는 모형이며, 진단 시점에 모든 PDAC 환자에게 적용하는 모형이 아니다.")
    content["methods"].insert(1, timing)
    content["methods"] += [
        "분화도를 고정 HR로 넣은 초기 확장 모형에서 비례위험 가정 위반을 확인하여, 분화도별 기저위험함수를 허용하는 층화 Cox 모형으로 보완하였다. 연령은 자연 스플라인 3 자유도로 반영하였다. 분화도 자체의 단일 HR을 계산하지 않으며 다른 변수의 효과와 분화도별 기저위험을 결합해 생존확률을 계산하였다.",
        f"유전체의 추가 예후 가치는 동일한 complete-case {models[0]['n']}명, 사망 {models[0]['events']}건에서 비교하였다. 임상병리 모형, KRAS 아형 추가, driver count 추가, KRAS 아형과 TP53·SMAD4·CDKN2A 추가의 4개 모형을 고정하였다. 공변량 결측 {models[0]['excluded_missing']}명은 해당 확장 Cox 및 예측모형에서 제외하며, 전체 코호트의 기술통계에는 1,011명을 유지하였다.",
        f"환자 단위 bootstrap {c_base['requested']}회마다 모형과 연령 스플라인을 다시 적합하였다. 36개월 사망확률 점수와 36개월 행정 검열을 이용한 Harrell C-index, 12·24·36개월 IPCW Brier score, 분화도 층화 calibration slope의 낙관성을 보정하였다. 동일 bootstrap 표본을 모든 모형에 사용하였다. OOB는 해당 환자가 학습 표본에 포함되지 않은 반복을 뜻한다. 보정도는 OOB 생존예측을 평균한 뒤 5분위별 KM 관찰 생존과 비교하였다.",
        "OOB 모형 성능 차이의 2.5–97.5백분위수는 재표집에 따른 변동성 범위이며, 낙관성 보정 차이의 95% 신뢰구간이 아니다. 단순 무작위 분할 결과를 외부 검증으로 부르지 않는다. 임상적 의사결정과 임계값이 확정되지 않아 decision curve나 새로운 치료 기준을 만들지 않았다.",
        "KRAS별 12·24·36개월 생존확률은 같은 임상병리 공변량 분포에 표준화하였다. 환자 bootstrap으로 점별 CI를 계산하였다. 모든 10쌍과 3시점의 총 30개 절대 차이에 대해 bootstrap 표준오차 기반 근사 검정과 Holm 보정을 제시하였다. 이는 유전자 아형을 바꾸었을 때의 인과효과가 아니다.",
        "병리 확장은 CDKN2A–림프절 양성(N1/2 대 N0), TP53–저분화(PD 대 WD/MD)의 두 로지스틱 모형으로 한정하였다. 연령 스플라인, 성별, T/M category와 선행치료를 보정하고 두 검정에 Holm 보정을 적용하였다. 표준화 확률 차이의 CI는 bootstrap으로 구하였다. 전체 범주의 기존 CMH 검정은 함께 유지하였다.",
        "선행치료별 KRAS 예후 연관성은 전체 상호작용 검정으로 평가하였다. M0만의 분석과 선행치료별 분석, 선행치료 전 CA19-9 및 ASA 추가 보정은 민감도 분석이다. 보조항암치료 yes/no는 수술 시점 고정 공변량으로 넣지 않았다. 재발 및 TMB 모형에는 노출과 log(time/12)의 상호작용을 추가하여 시간별 HR을 제시하였다.",
        "이번 확장 분석은 이미 관찰한 결과와 모형 진단을 바탕으로 계획한 탐색적 분석이다. 사전등록된 확증 분석으로 표현하지 않는다. 유의한 결과뿐 아니라 비유의 결과, 작은 예측 개선, 검증 실패 횟수와 한계를 함께 제시한다.",
        "HR은 비교군의 사건 순간위험률을 기준군의 순간위험률로 나눈 값이다. 예를 들어 A군 / B군 HR 1.43은 B군을 기준으로 A군의 순간위험률이 43% 높게 추정된다는 뜻이다. 사망확률이 43%p 증가하거나 PDAC가 발병할 확률이 43% 증가한다는 뜻은 아니다. OR 역시 확률비가 아니라 odds의 비이다."
    ]
    figs["Figure 2"]["assets"].append("Figure2D_Adjusted_pathology.png")
    figs["Figure 2"]["inputs"] += ["Supplementary_Table8_Adjusted_pathology.tsv", "Supplementary_Table8_Pathology_probabilities.tsv"]
    path_sentences=[]
    for r in path:
        target = "림프절 양성" if r["gene"] == "CDKN2A" else "저분화"
        path_sentences.append(f"{r['gene']} 변이 미검출군을 기준으로 검출군의 {target} OR은 {fmt(r['OR'])} (95% CI {fmt(r['lower95'])}–{fmt(r['upper95'])}), Holm p {pv(r['p_holm'])}이다. 표준화 확률 차이는 검출군에서 {number(r['probability_difference'])*100:+.1f}%p (95% CI {number(r['difference_lower95'])*100:+.1f}–{number(r['difference_upper95'])*100:+.1f})였다.")
    figs["Figure 2"]["finding"] = " ".join(path_sentences)
    figs["Figure 2"]["method"] += " 추가 패널은 두 병리 결과의 조정 OR와 표준화 확률을 보여준다. OR의 비교군은 변이 검출군이며 기준군은 미검출군이다."
    figs["Figure 2"]["interpretation"] = "같은 임상적 조건을 고려한 뒤에도 남는 연관성이 있는지를 OR과 절대 확률 차이로 평가하였다. OR이 1보다 크면 변이 검출군의 해당 병리 소견 odds가 더 높다. 확률 차이는 임상적으로 차이의 크기를 읽기 쉽게 한다. 이 결과와 생존모형의 추가 가치는 서로 다른 질문이므로 연결하되 동일시하지 않는다."
    figs["Figure 2"]["interpretation"] += " TP53의 전체 WD·MD·PD 분포 차이는 유의했지만, PD 대 WD/MD로 이분화한 보정 비교에서는 유의하지 않았다. 범주 정의와 보정 수준이 달라 같은 가설의 반복 검정이 아니며, TP53가 독립적으로 저분화를 예측한다고 단정하지 않는다."
    figs["Figure 3"]["method"] += " 분화도는 비례위험 위반을 고려해 층화 변수로 반영하였다."
    figs["Figure 3"]["interpretation"] += " 각 유전자 HR의 기준은 변이 미검출군이다. HR>1은 검출군의 사망 순간위험률 증가, HR<1은 감소 추정이며, q와 CI를 함께 보아야 한다."
    figs["Supplementary Figure 5"]["method"] += " 확장 모형에는 분화도 층화도 적용하였다. 모든 유전자 HR은 변이 검출군 / 미검출군의 사망 순간위험률 비이다."
    dr=get(pairs,group1="G12D",group2="G12R");dv=get(pairs,group1="G12D",group2="G12V")
    ad=get(absolute,group1="G12D",reference_group="G12R",month="36")
    figs["Figure 4"]["finding"] = hr_sentence(dr)+" "+hr_sentence(dv)+f" G12D−G12R의 표준화 36개월 생존확률 차이는 {number(ad['survival_difference'])*100:+.1f}%p (점별 95% CI {number(ad['lower95'])*100:+.1f}–{number(ad['upper95'])*100:+.1f})였다."
    figs["Figure 4"]["method"] += " 본문 Cox는 분화도 층화 및 확장 병리 보정 모형이다. 추가 패널은 표준화 생존곡선과 36개월 생존확률 차이를 제시한다."
    figs["Figure 4"]["interpretation"] = "KRAS 아형의 차이는 변이 보유 여부만으로는 표현되지 않는 예후의 이질성을 보여준다. HR에서 분자는 앞에 표시한 아형, 분모는 Reference로 표시한 아형이다. 반면 생존확률 차이는 앞 아형의 생존확률에서 기준 아형의 생존확률을 뺀 값으로, 양수이면 앞 아형의 생존이 더 좋다는 뜻이다. 두 지표의 방향을 혼동하지 않아야 한다."
    figs["Figure 4"]["limitation"] += f" 분화도 층화 후 전체 모형 PH 검정 p는 {pv(global_ph['p'])}이다. 유의하지 않은 PH 검정이 가정의 완전한 충족을 보증하지는 않는다."
    figs["Figure 4"]["assets"] += ["Figure4C_Standardized_survival.png", "Figure4D_Absolute_survival_difference.png"]
    figs["Figure 4"]["inputs"] += ["Supplementary_Table10_Adjusted_survival.tsv", "Supplementary_Table10_Survival_differences.tsv"]
    perftext = (f"낙관성 보정 36개월 C-index는 임상병리 단독 {fmt(c_base['corrected'],3)}, KRAS 아형 추가 {fmt(c_kras['corrected'],3)}, "
                f"driver count 추가 {fmt(c_count['corrected'],3)}, KRAS 아형 및 3개 유전자 추가 {fmt(c_genes['corrected'],3)}였다. "
                f"KRAS 추가의 C-index 변화는 {number(kr_delta['corrected_delta']):+.4f}였고, paired OOB 변화의 중앙 95% 범위는 {number(kr_delta['oob_p025']):+.4f}–{number(kr_delta['oob_p975']):+.4f}였다.")
    prediction_fig=dict(id="Figure 5", title="유전체 정보의 추가 예후 가치와 내부 검증",image="Figure5_Added_prognostic_value.png",
        assets=["Figure5A_Incremental_discrimination.png","Figure5B_Prediction_error.png","Figure5C_OOB_calibration.png"],
        method=f"동일한 {models[0]['n']}명에서 고정된 네 모형을 비교하고 {c_base['requested']}회 환자 bootstrap으로 낙관성을 보정하였다. C-index는 클수록, IPCW Brier score는 작을수록 좋다. 보정도 그림은 환자별 OOB 예측의 평균을 이용하며 대각선에 가까울수록 예측확률과 관찰 생존확률이 잘 맞는다.",
        finding=perftext,
        interpretation="유전체와 생존의 통계적 연관성이 확인되더라도 이미 강한 병리정보를 알고 있을 때 얻는 예측 개선은 작을 수 있다. 본 분석은 개선의 크기와 예측 오차, 보정도를 함께 보여준다. KRAS와 driver count의 우열은 각각의 p값으로 판단하지 않고 같은 환자에서의 성능 차이로 평가한다. 개선이 양수라는 이유만으로 임상적으로 유용하다고 확정하지 않는다.",
        limitation="이 검증은 동일 기관 자료에서 재표집한 내부 검증이다. OOB 변동성 범위는 보정 성능 차이의 CI가 아니다. 새 환자·다른 기관에서의 외부 검증이 필요하며, 이 모형만으로 항암제나 추적검사 간격을 변경하지 않는다. 결측 병리정보가 있는 환자를 제외한 결과이므로 선택 영향이 남는다.",
        inputs=["Supplementary_Table9_Prediction_performance.tsv","Supplementary_Table9_Prediction_differences.tsv","Supplementary_Table9_OOB_calibration.tsv","Supplementary_Table9_Model_global_tests.tsv","Input_Prediction_models.tsv"])
    content["figures"].insert(4,prediction_fig)
    content["figures"] += [dict(id="Supplementary Figure 6",title="임상 상황별 KRAS 예후 연관성",image="Supplementary_Figure6_Clinical_context.png",assets=["Supplementary_Figure6_Clinical_context.png"],
        method="전체 병리 보정, M0 민감도, 선행치료별 모형과 CA19-9·ASA 추가 모형에서 G12D 대 G12V/G12R의 HR을 표시하였다. 전체 아형 상호작용과 모형별 10개 쌍의 Holm p는 결과표에 보존하였다.",
        finding=f"KRAS 아형과 선행치료의 전체 상호작용 p는 {pv(interaction['p'])}이다. 그림의 모든 HR은 G12D의 사망 hazard를 범례에 표시한 기준 아형의 hazard로 나눈 값이다. HR>1이면 G12D의 위험이 더 높다.",
        interpretation="서로 다른 치료 맥락과 보정 조건에서 아형의 예후 연관성이 얼마나 일관되는지 평가한다. 한 군에서만 p<0.05라는 사실을 군 간 차이로 해석하지 않는다. 전체 상호작용의 크기와 불확실성을 고려하며, 비유의한 상호작용은 동등성의 증거가 아니다.",
        limitation="선행치료 수진은 무작위 배정이 아니며, 분석 대상은 수술에 도달한 환자들이다. 상호작용이 있더라도 약제 반응 예측이나 선행치료의 인과적 이득을 의미하지 않는다. M1은 전체 분석에 유지했고 M0 제한은 민감도 분석에만 적용하였다.",inputs=["Supplementary_Table11_Interaction.tsv","Supplementary_Table11_Sensitivity_pairs.tsv","Supplementary_Table11_Sensitivity_PH.tsv"]),
        dict(id="Supplementary Figure 7",title="재발과 TMB의 시간별 위험도",image="Supplementary_Figure7_Time_varying_effects.png",assets=["Supplementary_Figure7_Time_varying_effects.png"],
        method="기록된 재발과 Reported TMB의 모형에 노출×log(time/12) 항을 넣고 12·24·36개월의 HR과 점별 CI를 계산하였다. 선행치료의 시간가변 효과도 허용했으며 TMB 모형은 T category의 시간가변 효과를 추가하였다.",
        finding="시간별 HR과 95% CI를 표시하였다. 재발 모형의 비교군은 각 KRAS 아형, 기준군은 Not detected이다. TMB 모형은 1 SD 높은 값 대 낮은 값의 사망 HR이다. HR>1은 해당 시점에서 비교군의 사건 순간위험률이 더 높다는 뜻이다.",
        interpretation="단일 평균 HR로 설명하기 어려운 시간별 차이를 탐색한다. 시간에 따른 변화와 CI를 함께 읽으며 특정 시점의 p값만 선택하여 결론 내리지 않는다. 동일 시점의 순간위험률과 그 시점까지의 누적 재발·사망확률은 서로 다르다.",
        limitation="시간가변 효과의 log-time 형태 역시 모형 가정이다. 재발 전 사망과 재발 평가 종료 시점이 완전히 정리되지 않아 경쟁위험 누적발생 분석으로 확대하지 않았다. TMB 결측 원인과 측정 세부사항의 불확실성도 남아 있다.",inputs=["Supplementary_Table12_Time_varying_HR.tsv","Supplementary_Table13_TMB_missingness.tsv"])]
    # Individual panel pages keep each figure legible and unambiguous in Word.
    panel_titles={"Figure 2":["유전자와 병리 항목의 전체 연관성","TP53와 분화도 분포","CDKN2A와 림프절 병기 분포","보정된 병리 소견 확률"],
        "Figure 3":["주요 유전자와 전체생존","변이 유전자 수별 전체생존","변이 유전자 수의 보정 위험도"],
        "Figure 4":["KRAS 아형별 관찰 전체생존","KRAS 아형의 보정 쌍별 위험도","KRAS 아형의 표준화 생존곡선","36개월 생존확률의 절대 차이"],
        "Figure 5":["유전체 정보 추가에 따른 판별력","36개월 생존 예측의 오차","OOB 예측의 보정도"]}
    for f in content["figures"]:
        if f["id"] in panel_titles:
            f["panels"]=[dict(id=f"{f['id']}{chr(97+i)}",title=panel_titles[f['id']][i],image=name,
                caption=f["method"],inputs=f["inputs"][:3]) for i,name in enumerate(f["assets"])]
    panel_details={
        "Figure 2a":("전체 24개 유전자–병리 연관성을 선행치료 여부로 층화한 generalized CMH 검정으로 평가하였다. 각 칸은 원 p와 BH q를 표시한다.",["Supplementary_Table2_Pathology_CMH.tsv"]),
        "Figure 2b":("TP53 변이 검출군과 미검출군의 WD·MD·PD 분포이다. 백분율은 해당 유전자 상태 집단 내 비율이며, 그림 아래에는 Holm 보정 후 유의한 범주 쌍을 표시하였다.",["Input_Pathology_counts.tsv","Supplementary_Table2_Pathology_pairs.tsv"]),
        "Figure 2c":("CDKN2A 변이 검출군과 미검출군의 N0·N1·N2 분포이다. 해당 상태 집단 내 비율과 Holm 보정 후 유의한 범주 쌍을 표시하였다.",["Input_Pathology_counts.tsv","Supplementary_Table2_Pathology_pairs.tsv"]),
        "Figure 2d":("임상 변수 보정 후 표준화한 림프절 양성 및 저분화 확률과 점별 bootstrap 95% CI이다. OR은 변이 검출군 / 미검출군이며 두 검정에 Holm 보정을 적용하였다.",["Supplementary_Table8_Adjusted_pathology.tsv","Supplementary_Table8_Pathology_probabilities.tsv"]),
        "Figure 3a":("각 유전자의 검출군 / 미검출군 사망 HR이다. 분화도 층화 및 확장 임상병리 보정을 적용하고 원 p와 BH q를 표시하였다.",["Supplementary_Table3_Cox_coefficients.tsv","Supplementary_Table3_PH_tests.tsv"]),
        "Figure 3b":("KRAS·TP53·SMAD4·CDKN2A 중 변이 유전자 수에 따른 관찰 KM 곡선이다. 전체 및 쌍별 log-rank와 number at risk를 함께 제시하였다.",["Supplementary_Table4_Logrank.tsv","Supplementary_Table4_Logrank_pairs.tsv","Supplementary_Table4_Number_at_risk.tsv"]),
        "Figure 3c":("변이 유전자 수 1개를 기준으로 0개·2개·3개 이상의 사망 HR을 비교하였다. 분화도 층화와 확장 임상병리 보정 모형이며 HR>1은 표시된 개수 집단의 사망 hazard가 높다는 뜻이다.",["Supplementary_Table3_Cox_coefficients.tsv"]),
        "Figure 4a":("임상정보 KRAS 아형별 관찰 KM 전체생존곡선이다. OS 시작점은 수술이며 전체 log-rank와 유의한 Holm 보정 쌍별 비교, number at risk를 표시하였다.",["Supplementary_Table4_Logrank.tsv","Supplementary_Table4_Logrank_pairs.tsv","Supplementary_Table4_Number_at_risk.tsv"]),
        "Figure 4b":("분화도 층화 및 확장 임상병리 보정 Cox의 10쌍 비교이다. HR은 앞 아형 / Reference 아형의 사망 hazard이며 95% CI와 Holm p를 표시한다.",["Supplementary_Table4_Cox_pairs.tsv","Supplementary_Table3_PH_tests.tsv"]),
        "Figure 4c":("같은 1,000명 공변량 분포로 표준화한 아형별 모형 기반 생존곡선이다. 관찰 KM 또는 인과효과가 아니며 12·24·36개월 점별 bootstrap CI는 결과표에 제시하였다.",["Input_Standardized_survival_curves.tsv","Supplementary_Table10_Adjusted_survival.tsv"]),
        "Figure 4d":("36개월 표준화 생존확률의 앞 아형−기준 아형 차이(%p)와 점별 bootstrap CI이다. 양수는 앞 아형의 생존이 높음을 뜻한다. Holm p는 10쌍×3시점의 30개 비교 전체를 보정하였다.",["Supplementary_Table10_Survival_differences.tsv"]),
        "Figure 5a":("동일 환자의 4개 모형에서 500회 bootstrap 낙관성 보정 36개월 Harrell C-index를 비교하였다. 작은 차이를 과장하지 않도록 0.50–0.80 척도로 표시한다.",["Supplementary_Table9_Prediction_performance.tsv","Supplementary_Table9_Prediction_differences.tsv"]),
        "Figure 5b":("동일 환자의 낙관성 보정 36개월 IPCW Brier score이다. 낮을수록 생존확률 예측의 평균 오차가 작다. 점추정이며 외부 검증 성능이 아니다.",["Supplementary_Table9_Prediction_performance.tsv"]),
        "Figure 5c":("학습에 포함되지 않은 환자의 OOB 예측을 평균하여 5분위별 관찰 KM 생존확률과 비교하였다. 대각선은 완전한 일치이며 막대는 관찰 생존의 점별 95% CI이다.",["Supplementary_Table9_OOB_calibration.tsv","Input_Validation_coverage.tsv"]),
    }
    for f in content["figures"]:
        for panel in f.get("panels",[]):
            if panel["id"] in panel_details:
                panel["caption"],panel["inputs"]=panel_details[panel["id"]]
    # Explicit comparison/reference table replaces ambiguous 'A vs B'.
    st4=next(t for t in content["tables"] if t["id"]=="Supplementary Table 4")
    st4.update(headers=["Comparison group","Reference group","Death HR [95% CI]","Holm p","Estimate direction"],
        rows=[[r["group1"],r["group2"],f"{fmt(r['HR'])} [{fmt(r['lower95'])}, {fmt(r['upper95'])}]",pv(r["p_holm"]),
               f"비교군에서 {abs(number(r['HR'])-1)*100:.1f}% {'높음' if number(r['HR'])>=1 else '낮음'}"] for r in pairs],
        note="HR=비교군 사망 hazard / 기준군 사망 hazard. 방향 열은 점추정의 방향이며 유의성 주장이 아니다. 10쌍 Holm 보정, 점별 95% CI. 누적 사망확률 또는 PDAC 발병위험의 비율과 다르다.")
    st3=next(t for t in content["tables"] if t["id"]=="Supplementary Table 3")
    st3["note"]="분화도 층화와 확장 임상병리 보정 모형이다. 각 유전자는 검출군 / 미검출군, driver count는 표시 개수 / 1개이다. 전체 계수표에 비교군·기준군·사건 정의·해석을 기록하였다. q는 동일 보정 수준의 유전체 노출 계수에 대한 BH 보정이다."
    content["tables"] += [
        tab("Supplementary Table 8","병리 소견의 보정 연관성",["Gene and outcome","Detected / not detected OR [95% CI]","Holm p","Probability difference (%p)"],
            [[r['gene']+" / "+r['outcome'],f"{fmt(r['OR'])} [{fmt(r['lower95'])}, {fmt(r['upper95'])}]",pv(r['p_holm']),f"{number(r['probability_difference'])*100:+.1f} [{number(r['difference_lower95'])*100:+.1f}, {number(r['difference_upper95'])*100:+.1f}]"] for r in path],
            ["Supplementary_Table8_Adjusted_pathology.tsv","Supplementary_Table8_Pathology_probabilities.tsv"],"OR의 비교군은 변이 검출, 기준군은 미검출이다. 확률 차이는 검출−미검출이다. 두 모형에 Holm 보정, 확률 차이에 점별 bootstrap CI를 적용하였다."),
        tab("Supplementary Table 9","유전체 추가 모형의 내부 검증",["Model","N","Corrected C at 36 months","Corrected Brier at 36 months","Corrected calibration slope"],
            [[labels[m],models[0]['n'],fmt(get(performance,model=m,metric="C_index")['corrected'],3),fmt(get(performance,model=m,metric="Brier_36")['corrected'],3),fmt(get(performance,model=m,metric="Calibration_slope")['corrected'],3)] for m in labels],
            prediction_fig['inputs'],"C는 클수록, Brier는 작을수록 좋고 slope의 이상값은 1이다. 500회 bootstrap 낙관성 보정. 분화도별 기저위험을 허용하며 외부 검증은 아니다. 전체 시점·OOB 변동성·실패 수는 입력표에 보존하였다."),
        tab("Supplementary Table 10","KRAS 아형별 표준화 생존확률",["KRAS subtype","Month","Survival % [95% CI]"],
            [[r['group'],r['month'],f"{number(r['survival'])*100:.1f} [{number(r['lower95'])*100:.1f}, {number(r['upper95'])*100:.1f}]"] for r in survival],
            ["Supplementary_Table10_Adjusted_survival.tsv","Supplementary_Table10_Survival_differences.tsv"],"같은 임상병리 분포로 표준화한 예측 생존확률이다. 점별 bootstrap CI이며 개별 환자 예측구간이 아니다. 전체 30개 확률 차이와 Holm p는 Input 2에 제시하였다."),
        tab("Supplementary Table 11","선행치료와 KRAS 아형의 상호작용",["Interaction","N","Global p"],[["KRAS subtype × neoadjuvant",interaction['n'],pv(interaction['p'])]],
            ["Supplementary_Table11_Interaction.tsv","Supplementary_Table11_Sensitivity_pairs.tsv","Supplementary_Table11_Sensitivity_PH.tsv"],"전체 아형의 상호작용 검정이다. 수술 환자에서의 예후 이질성을 평가하며 선행치료의 인과적 효과 또는 치료반응 예측을 의미하지 않는다."),
        tab("Supplementary Table 12","시간별 사건 위험도",["Endpoint and comparison / reference","Month","HR [95% CI]","Holm p"],
            [[r['endpoint']+" / "+r['comparison']+" / ref: "+r['reference'],r['month'],f"{fmt(r['HR'])} [{fmt(r['lower95'])}, {fmt(r['upper95'])}]",pv(r['p_holm'])] for r in tv],
            ["Supplementary_Table12_Time_varying_HR.tsv"],"HR>1이면 표시된 비교군의 해당 시점 사건 hazard가 기준군보다 높다. 사건은 기록된 재발 또는 전체 사망으로 구분한다. 점별 Wald CI와 endpoint 내 Holm p를 제시하였다."),
        tab("Supplementary Table 13","TMB 결측에 따른 환자 구성",["TMB status","N","Deaths n","Median reverse KM follow-up months"],
            [[r['TMB_status'],r['n'],r['deaths'],fmt(r['median_reverse_KM_followup'],1)] for r in missing],
            ["Supplementary_Table13_TMB_missingness.tsv"],"사망 건수는 동일 추적기간의 위험 비교가 아니다. 결측과 관측 집단의 추적기간 차이를 함께 보고하며 TMB 결측을 0이나 Variant count로 대체하지 않았다.")]
    content["references"] += [
        {"label":"McIntyre CA et al. Cancer Cell. 2024;42:1614–1629.e5.","url":"https://doi.org/10.1016/j.ccell.2024.08.002"},
        {"label":"Prognostic Implications of Codon-Specific KRAS Mutations. JCO Precis Oncol. 2026.","url":"https://doi.org/10.1200/PO-25-00115"},
        {"label":"Collins GS et al. TRIPOD+AI statement. BMJ. 2024;385:e078378.","url":"https://doi.org/10.1136/bmj-2023-078378"}]
    novelty = ("KRAS 아형의 예후 차이 자체는 기존 연구에서도 보고되었다 [7,8]. 본 연구의 차별화된 질문은 상세한 수술 병리정보를 이미 알고 있을 때 유전체 정보가 추가로 제공하는 예후적 가치의 크기이다. "
               f"KRAS 아형 추가의 낙관성 보정 C-index 변화는 {number(kr_delta['corrected_delta']):+.4f}, 아형과 3개 유전자를 함께 추가한 변화는 {number(gene_delta['corrected_delta']):+.4f}였다. "
               "통계적으로 유의한 HR과 실제 예측 개선을 분리해 제시한 점이 임상적 기여이다. 이 결과만으로 새로운 진료 기준이나 외부 임상 유용성을 확정하지 않는다.")
    content["discussion"] = [
        dict(title="병리 표현형과 유전체의 연결",text=" ".join(path_sentences)+" 병리 소견과의 연관성이 곧 독립적인 생존 영향이나 전이의 인과기전을 의미하지는 않는다."),
        dict(title="KRAS 아형의 예후 차이",text=hr_sentence(dr)+" "+hr_sentence(dv)+" 모든 HR은 표시된 비교군을 기준군과 비교한다. 이 수치는 누적 사망확률의 상대 증가가 아니다."),
        dict(title="유전체의 추가 예후 가치",text=perftext+" 판별력뿐 아니라 Brier score와 OOB 보정도를 함께 검토해야 한다. 개선이 작거나 불확실하면 그 제한을 연구 결과로 보고한다."),
        dict(title="기존 연구와 연구의 기여",text=novelty+" 1,011명은 McIntyre 2020의 283명과 Campbell 2025의 508명보다 많지만, 전체 2,336명의 Varghese 연구도 있어 최대 규모를 주장하지 않는다 [1–3]."),
        dict(title="적용 범위와 남은 한계",text=timing+" NGS 대상 선택과 잔여 교란, 유전자 검사 범위 및 변이 병원성, TMB 결측 원인은 여전히 제한이다. 이번 결과 기반 확장과 모형 진단 후 수정은 탐색적 성격을 가지며 외부 검증이 필요하다. 재발은 사망 포함 DFS나 경쟁위험 누적발생률이 아니다."),
        dict(title="Summary",text="1,011명의 PDAC 수술 환자에서 유전체와 병리 소견, 생존을 연결하였다. "+hr_sentence(dr)+" "+
             f"임상병리 단독과 KRAS 추가의 내부 검증 C-index는 각각 {fmt(c_base['corrected'],3)}, {fmt(c_kras['corrected'],3)}였다. "
             "연구의 임상적 의미는 아형별 예후 연관성과 실제 추가 예측 가치의 크기를 구분해 제시하는 데 있다. 환자 수와 사건 수, 상세 병리 보정, 일관된 비교군 정의와 내부 검증이 강점이다. 치료 선택이나 추적 전략 변경을 뒷받침하려면 별도 외부 검증과 임상적 유용성 평가가 필요하다.")]
    content["subtitle"]="1,011명 코호트의 병리 연관성과 유전체의 추가 예후 가치"
    content["main_figure_count"]=5
    content.pop("review_note",None)
    return content
