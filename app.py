import re
import math
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="Lab-ready His-tag Toolkit", layout="wide")

# ==============================
# 0) Global constants & helpers
# ==============================
# Average residue masses (Da) for polypeptides (residue form)
AA_RESIDUE_MW = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09,
    "C": 103.15, "E": 129.12, "Q": 128.13, "G": 57.05,
    "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
    "M": 131.19, "F": 147.18, "P": 97.12,  "S": 87.08,
    "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13
}
H2O = 18.015  # Da


def clean_sequence(raw: str) -> str:
    """Accept FASTA or raw. Remove headers and keep only valid AA letters."""
    if not raw:
        return ""
    lines = []
    for line in raw.splitlines():
        if line.strip().startswith(">"):
            continue
        lines.append(line)
    seq = "".join(lines).upper()
    seq = re.sub(r"[^A-Z]", "", seq)
    seq = "".join([aa for aa in seq if aa in AA_RESIDUE_MW])
    return seq


def calc_mw_da(seq: str) -> float:
    """Residue sum + H2O(termini)"""
    if not seq:
        return 0.0
    mw = sum(AA_RESIDUE_MW[aa] for aa in seq) + H2O
    return round(mw, 2)


def find_his_runs(seq: str, min_len: int = 6):
    return [(m.start(), m.end()) for m in re.finditer(rf"H{{{min_len},}}", seq)]


def c1v1_ml(final_vol_ml: float, final_conc, stock_conc) -> float:
    """C1V1=C2V2; units must match. Returns mL."""
    if stock_conc == 0:
        return 0.0
    return (final_conc * final_vol_ml) / stock_conc


def pct_vv_ml(final_vol_ml: float, final_pct: float, stock_pct: float) -> float:
    """v/v percent from percent stock. Returns mL."""
    if stock_pct == 0:
        return 0.0
    return (final_pct / stock_pct) * final_vol_ml


def recipe_df(final_vol_ml: float, rows: list[dict]) -> pd.DataFrame:
    used = sum(float(r["vol_ml"]) for r in rows)
    dw = max(0.0, float(final_vol_ml) - used)
    out_rows = []
    for r in rows:
        out_rows.append({
            "Component": r["Component"],
            "Stock": r.get("Stock", "-"),
            "Target (final)": r.get("Target (final)", "-"),
            "Volume (mL)": float(r["vol_ml"])
        })
    out_rows.append({"Component": "DW", "Stock": "-", "Target (final)": "to volume", "Volume (mL)": dw})
    df = pd.DataFrame(out_rows)
    df["Volume (mL)"] = df["Volume (mL)"].map(lambda x: round(x, 4))
    return df


def download_df_button(df: pd.DataFrame, filename: str, label: str):
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(label=label, data=csv, file_name=filename, mime="text/csv")


# ==============================
# Sidebar: shared inputs
# ==============================
st.sidebar.title("⚙️ Experiment Inputs")

pellet_g = st.sidebar.number_input("E. coli pellet (g)", min_value=0.1, max_value=20.0, value=3.0, step=0.1)
resin_ul = st.sidebar.number_input("Ni-NTA resin (µL) — 50% slurry", min_value=10, max_value=2000, value=100, step=10)
st.sidebar.caption("50% slurry = resin:20% EtOH = 1:1 (실제 resin volume은 slurry의 50%)")

page = st.sidebar.radio(
    "📄 Pages",
    [
        "Protein MW",
        "Ni-NTA Buffer (stock-based)",
        "SOP & Checklist",
        "Virtual SDS-PAGE",
    ],
)


# ==============================
# Page 1: Protein MW
# ==============================

def page_protein_mw():
    st.title("🧬 Protein MW Calculator")

    colA, colB = st.columns([1.15, 0.85], gap="large")

    with colA:
        raw = st.text_area(
    "Protein sequence (FASTA 가능)",
    height=240,
    placeholder=">protein_name\nMHHHHHH..."
)


        c1, c2, c3 = st.columns(3)
        with c1:
            add_his = st.checkbox("N-말단 6×His 추가해서 계산", value=False)
        with c2:
            detect_his = st.checkbox("His-run(≥6) 위치 표시", value=True)
        with c3:
            show_len = st.checkbox("Length 표시", value=True)

        seq = clean_sequence(raw)
        seq_calc = ("HHHHHH" + seq) if (add_his and seq) else seq

        if not seq_calc:
            st.info("서열을 입력하면 MW가 계산됩니다.")
            return

        mw_da = calc_mw_da(seq_calc)
        mw_kda = mw_da / 1000.0

        st.metric("Molecular Weight", f"{mw_da:,.2f} Da", f"{mw_kda:,.3f} kDa")
        if show_len:
            st.write(f"Length: **{len(seq_calc)} aa**")

        if detect_his:
            runs = find_his_runs(seq_calc, 6)
            if runs:
                st.success(f"His-run(≥6) 탐지: {len(runs)}개 | 위치(0-index): {runs}")
            else:
                st.info("His-run(≥6) 없음")

        st.divider()
        st.subheader("실험 팁")
        if mw_kda < 10:
            st.warning("10 kDa 미만이면 Tricine-SDS-PAGE가 밴드 분리/식별에 유리합니다.")
        else:
            st.info("일반 Tris-Glycine gel에서도 확인 가능할 확률이 큽니다.")

    with colB:
        st.subheader("Quick Notes")
        st.markdown(
            f"""
- 입력은 FASTA/일반 텍스트 모두 OK (헤더 자동 제거)
- 알파벳 외 문자 자동 제거
- 6×His는 옵션(자동 추가)

**현재 입력 조건(고정 컨텍스트)**
- Pellet: **{pellet_g} g**
- Resin: **{resin_ul} µL (50% slurry)**
"""
        )


# ==============================
# Page 2: Buffer designer
# ==============================

def page_buffer_designer():
    st.title("🧪 Ni-NTA Buffer Designer (Stock-based)")
    st.caption("Stock 기준(C1V1)으로 buffer 조성(물질/vol/final vol)을 자동 계산합니다.")

    base_options = {
        "0.5 M Phosphate buffer (pH 7.0)": 500.0,
        "1 M HEPES (pH 7.4)": 1000.0,
        "1 M Tris-HCl (pH 8.0)": 1000.0,
        "1 M Tris-HCl (pH 8.5)": 1000.0,
    }

    with st.sidebar.expander("🧴 Stock 선택", expanded=True):
        base_stock_name = st.selectbox("Base buffer stock", list(base_options.keys()), index=2)
        nacl_stock_mM = st.number_input("NaCl stock (mM) [5M=5000mM]", value=5000.0, step=100.0)
        tx_stock_pct = st.number_input("Triton X-100 stock (%) [10%]", value=10.0, step=1.0)
        pi_stock_x = st.number_input("Protease inhibitor stock (×) [100×]", value=100.0, step=10.0)
        imid_stock_mM = st.number_input("Imidazole stock (mM) [1M=1000mM]", value=1000.0, step=50.0)
        dtt_stock_mM = st.number_input("DTT stock (mM) [1M=1000mM]", value=1000.0, step=50.0)
        bme_stock_mM = st.number_input("β-ME stock (mM) [14.3M=14300mM]", value=14300.0, step=100.0)

    st.subheader("1) Volume 설정")
    strategy = st.radio(
        "Volume strategy",
        [
            "Use SOP defaults (Lysis 30 mL / Wash 20 mL / Elution 10 mL)",
            "Auto-scale from pellet (Lysis = pellet × factor)",
            "Manual",
        ],
        horizontal=False,
    )

    if strategy == "Auto-scale from pellet (Lysis = pellet × factor)":
        factor = st.number_input("Lysis factor (mL per g pellet)", min_value=1.0, max_value=50.0, value=10.0, step=1.0)
        lysis_vol = float(pellet_g) * float(factor)
        wash_vol = st.number_input("Wash buffer (mL)", min_value=2.0, max_value=500.0, value=20.0, step=5.0)
        elution_vol = st.number_input("Elution buffer (mL)", min_value=2.0, max_value=200.0, value=10.0, step=2.0)
    elif strategy == "Manual":
        lysis_vol = st.number_input("Lysis buffer (mL)", min_value=5.0, max_value=500.0, value=30.0, step=5.0)
        wash_vol = st.number_input("Wash buffer (mL)", min_value=2.0, max_value=500.0, value=20.0, step=5.0)
        elution_vol = st.number_input("Elution buffer (mL)", min_value=2.0, max_value=200.0, value=10.0, step=2.0)
    else:
        lysis_vol, wash_vol, elution_vol = 30.0, 20.0, 10.0
        st.info("SOP 기본값: Lysis 30 mL / Wash 20 mL / Elution 10 mL")

    st.divider()
    st.subheader("2) Target 조성 (final)")

    cA, cB, cC = st.columns(3)
    with cA:
        base_final_mM = st.number_input("Base buffer final (mM)", value=50.0, step=5.0)
        nacl_final_mM = st.number_input("NaCl final (mM)", value=300.0, step=25.0)
    with cB:
        tx_final_pct = st.number_input("Triton X-100 final (%)", value=1.0, step=0.1)
        pi_final_x = st.number_input("Protease inhibitor final (×)", value=1.0, step=1.0)
    with cC:
        imid_lysis_mM = st.number_input("Imidazole in Lysis (mM)", value=0.0, step=5.0)
        imid_wash_mM = st.number_input("Imidazole in Wash (mM)", value=20.0, step=5.0)
        imid_elution_mM = st.number_input("Imidazole in Elution (mM)", value=300.0, step=25.0)

    st.divider()
    st.subheader("3) Reducing agent 옵션")
    r1, r2 = st.columns(2)
    with r1:
        add_dtt = st.checkbox("Add DTT", value=False)
        dtt_final_mM = st.number_input("DTT final (mM)", value=1.0, step=0.5, disabled=not add_dtt)
    with r2:
        add_bme = st.checkbox("Add β-ME", value=False)
        bme_final_mM = st.number_input("β-ME final (mM)", value=5.0, step=1.0, disabled=not add_bme)

    st.divider()

    if st.button("✅ Calculate buffer recipes", type="primary"):
        base_stock_mM = base_options[base_stock_name]

        # Lysis
        l_rows = []
        l_rows.append({"Component": "Base buffer", "Stock": base_stock_name,
                       "Target (final)": f"{base_final_mM} mM",
                       "vol_ml": c1v1_ml(lysis_vol, base_final_mM, base_stock_mM)})
        l_rows.append({"Component": "NaCl", "Stock": f"{nacl_stock_mM/1000:.1f} M",
                       "Target (final)": f"{nacl_final_mM} mM",
                       "vol_ml": c1v1_ml(lysis_vol, nacl_final_mM, nacl_stock_mM)})
        if tx_final_pct > 0:
            l_rows.append({"Component": "Triton X-100", "Stock": f"{tx_stock_pct}%",
                           "Target (final)": f"{tx_final_pct}%",
                           "vol_ml": pct_vv_ml(lysis_vol, tx_final_pct, tx_stock_pct)})
        if pi_final_x > 0:
            l_rows.append({"Component": "Protease inhibitor", "Stock": f"{pi_stock_x}×",
                           "Target (final)": f"{pi_final_x}×",
                           "vol_ml": (pi_final_x / pi_stock_x) * lysis_vol})
        if imid_lysis_mM > 0:
            l_rows.append({"Component": "Imidazole", "Stock": f"{imid_stock_mM/1000:.1f} M",
                           "Target (final)": f"{imid_lysis_mM} mM",
                           "vol_ml": c1v1_ml(lysis_vol, imid_lysis_mM, imid_stock_mM)})
        if add_dtt and dtt_final_mM > 0:
            l_rows.append({"Component": "DTT", "Stock": f"{dtt_stock_mM/1000:.1f} M",
                           "Target (final)": f"{dtt_final_mM} mM",
                           "vol_ml": c1v1_ml(lysis_vol, dtt_final_mM, dtt_stock_mM)})
        if add_bme and bme_final_mM > 0:
            l_rows.append({"Component": "β-ME", "Stock": "14.3 M",
                           "Target (final)": f"{bme_final_mM} mM",
                           "vol_ml": c1v1_ml(lysis_vol, bme_final_mM, bme_stock_mM)})

        lysis_df = recipe_df(lysis_vol, l_rows)

        # Wash
        w_rows = []
        w_rows.append({"Component": "Base buffer", "Stock": base_stock_name,
                       "Target (final)": f"{base_final_mM} mM",
                       "vol_ml": c1v1_ml(wash_vol, base_final_mM, base_stock_mM)})
        w_rows.append({"Component": "NaCl", "Stock": f"{nacl_stock_mM/1000:.1f} M",
                       "Target (final)": f"{nacl_final_mM} mM",
                       "vol_ml": c1v1_ml(wash_vol, nacl_final_mM, nacl_stock_mM)})
        w_rows.append({"Component": "Imidazole", "Stock": f"{imid_stock_mM/1000:.1f} M",
                       "Target (final)": f"{imid_wash_mM} mM",
                       "vol_ml": c1v1_ml(wash_vol, imid_wash_mM, imid_stock_mM)})
        wash_df = recipe_df(wash_vol, w_rows)

        # Elution
        e_rows = []
        e_rows.append({"Component": "Base buffer", "Stock": base_stock_name,
                       "Target (final)": f"{base_final_mM} mM",
                       "vol_ml": c1v1_ml(elution_vol, base_final_mM, base_stock_mM)})
        e_rows.append({"Component": "NaCl", "Stock": f"{nacl_stock_mM/1000:.1f} M",
                       "Target (final)": f"{nacl_final_mM} mM",
                       "vol_ml": c1v1_ml(elution_vol, nacl_final_mM, nacl_stock_mM)})
        e_rows.append({"Component": "Imidazole", "Stock": f"{imid_stock_mM/1000:.1f} M",
                       "Target (final)": f"{imid_elution_mM} mM",
                       "vol_ml": c1v1_ml(elution_vol, imid_elution_mM, imid_stock_mM)})
        elution_df = recipe_df(elution_vol, e_rows)

        st.success("Buffer recipes calculated.")

        t1, t2, t3 = st.tabs([f"Lysis ({lysis_vol:.1f} mL)", f"Wash ({wash_vol:.1f} mL)", f"Elution ({elution_vol:.1f} mL)"])
        with t1:
            st.dataframe(lysis_df, use_container_width=True)
            download_df_button(lysis_df, "lysis_buffer.csv", "⬇️ Download Lysis CSV")
        with t2:
            st.dataframe(wash_df, use_container_width=True)
            download_df_button(wash_df, "wash_buffer.csv", "⬇️ Download Wash CSV")
        with t3:
            st.dataframe(elution_df, use_container_width=True)
            download_df_button(elution_df, "elution_buffer.csv", "⬇️ Download Elution CSV")

        st.divider()
        st.subheader("현장용 요약")
        st.markdown(
            f"""
- Pellet: **{pellet_g} g** | Resin slurry: **{resin_ul} µL** (실제 resin ~{resin_ul/2:.0f} µL)
- Base: **{base_stock_name}** → final **{base_final_mM} mM**
- NaCl final **{nacl_final_mM} mM**
- Imidazole: Wash **{imid_wash_mM} mM** / Elution **{imid_elution_mM} mM**
"""
        )


# ==============================
# Page 3: SOP & Checklist
# ==============================

def page_sop():
    st.title("📋 SOP & Checklist")
    st.caption("실험할 때 그대로 따라가기 쉬운 형태로 정리했습니다.")

    st.subheader("1) Buffer Recipes (SOP 예시)")
    st.markdown(
        """
**1.1. Lysis Buffer (Final Vol: 30 mL)**
- 1 M Tris-HCl (pH 8.0): 1.5 mL
- 5 M NaCl: 1.8 mL
- 50% Glycerol: 6.0 mL
- 10% Triton X-100: 3.0 mL
- 100× Protease Inhibitor: 0.3 mL
- DW: 17.4 mL
- Note: 50mM Tris, 300mM NaCl, 10% Glycerol, 1% Tx-100

**1.2. Wash Buffer (Final Vol: 20 mL)**
- 1 M Tris-HCl (pH 8.0): 1.0 mL
- 5 M NaCl: 1.2 mL
- 1 M Imidazole: 0.4 mL (Final 20mM)
- DW: 17.4 mL
- Note: 50mM Tris, 300mM NaCl, 20mM Imidazole

**1.3. Elution Buffer (Final Vol: 10 mL)**
- 1 M Tris-HCl (pH 8.0): 0.5 mL
- 5 M NaCl: 0.6 mL
- 1 M Imidazole: 3.0 mL (Final 300mM)
- DW: 5.9 mL
- Note: 50mM Tris, 300mM NaCl, 300mM Imidazole
"""
    )

    st.divider()
    st.subheader("2) Detailed Protocol")

    with st.expander("2.1 Lysis & Clarification", expanded=True):
        st.markdown(
            """
- Cell Harvesting: 4,000 rpm, 4°C, 20 min → pellet 3 g
- Resuspension: Lysis Buffer 30 mL로 완전 현탁
- Sonication: Ice 위, 30% Amp, 5s On / 10s Off, total 10 min
- Clarification: 13,000 rpm, 4°C, 30 min
- Fraction: 상층액(soluble) 회수 + SDS-PAGE용 Total/Soluble 샘플 보관
"""
        )

    with st.expander("2.2 Ni-NTA Purification", expanded=True):
        st.markdown(
            f"""
- Resin Preparation: Ni-NTA slurry {resin_ul} µL (실제 resin 약 {resin_ul/2:.0f} µL)
- Equilibration: Lysis Buffer 2 mL
- Loading: 상층액 전체 로딩 (gravity flow) → flow-through 보관
- Washing: Wash Buffer 20 mL (분할 로딩) → 마지막 wash 일부 보관
- Elution: Elution Buffer 10 mL, 1 mL × 10 fraction 권장
"""
        )

    with st.expander("2.3 Analysis", expanded=True):
        st.markdown(
            """
- SDS-PAGE: Total / Soluble / Flow-through / Wash / Elution 로딩
- 저분자(예: IGF-1 ~7.6 kDa)면 Tricine-SDS-PAGE 추천
"""
        )

    st.divider()
    st.subheader("3) Checklist")
    st.markdown(
        """
**실험 전**
- [ ] Buffer 라벨(L/W/E) + 충분한 여유분
- [ ] 샘플링 튜브(총/soluble/FT/wash/elution) 라벨
- [ ] 얼음/소닉/원심 조건 확인

**실험 중**
- [ ] 소닉 과열 방지(ice 유지)
- [ ] 원심 후 pellet disturbance 방지
- [ ] Elution fraction 분취(해석 용이)

**실험 후**
- [ ] SDS-PAGE 로딩 조건 정리
- [ ] Elution pool 기준(순도/수율) 기록
"""
    )


# ==============================
# Page 4: Virtual SDS-PAGE
# ==============================

def _mw_to_migration(mw_kda: float, gel_system: str, gel_pct: float) -> float:
    """Simple log mapping + gel% effect (higher % -> smaller proteins run further). Returns 0~1."""
    # base ranges
    if gel_system == "Tris-Glycine":
        hi, lo = 250.0, 10.0
    else:  # Tricine
        hi, lo = 100.0, 1.0

    mw_kda = min(max(mw_kda, lo), hi)

    # log distance (0 top, 1 bottom)
    base = (math.log10(hi) - math.log10(mw_kda)) / (math.log10(hi) - math.log10(lo))

    # gel% tweak: relative to 12%
    # higher %: pushes small proteins further down (increase base a bit)
    pct_effect = (gel_pct - 12.0) * 0.015

    # keep within 0~1
    y = base + pct_effect
    return min(max(y, 0.02), 0.98)


def _default_ladders():
    return {
        "10–250 kDa ladder": [250, 150, 100, 75, 50, 37, 25, 20, 15, 10],
        "5–100 kDa ladder": [100, 75, 50, 37, 25, 20, 15, 10, 5],
        "2–40 kDa ladder": [40, 30, 25, 20, 15, 10, 5, 2],
    }


def _recommend_gel(mw_kda: float):
    # practical defaults
    if mw_kda <= 10:
        return "Tricine", 16.5
    if mw_kda <= 25:
        return "Tris-Glycine", 15.0
    if mw_kda <= 60:
        return "Tris-Glycine", 12.0
    return "Tris-Glycine", 10.0


def page_virtual_sds():
    st.title("🧫 Virtual SDS-PAGE (Band Predictor)")
    st.caption("예상 MW 기반으로 밴드 위치를 사전 예측하고, gel%/marker 선택 근거를 확보합니다.")

    left, right = st.columns([1.05, 0.95], gap="large")

    with left:
        raw = st.text_area("Protein sequence (FASTA 가능)", height=180)
        seq = clean_sequence(raw)

        manual_mw = st.checkbox("서열 없이 MW만으로 시뮬레이션", value=False)
        add_his = st.checkbox("N-말단 6×His 포함", value=False, disabled=manual_mw)

        if manual_mw:
            mw_kda = st.number_input("Target MW (kDa)", min_value=0.5, max_value=300.0, value=7.6, step=0.1)
        else:
            if not seq:
                st.info("서열을 넣거나 'MW만으로 시뮬레이션'을 켜세요.")
                return
            seq_calc = ("HHHHHH" + seq) if add_his else seq
            mw_kda = calc_mw_da(seq_calc) / 1000.0

        st.metric("Predicted MW", f"{mw_kda:.3f} kDa")

    with right:
        ladders = _default_ladders()
        ladder_name = st.selectbox("Marker ladder", list(ladders.keys()), index=1)

        # gel system + gel%
        auto_gel = st.checkbox("MW 기반 gel 추천값 자동 적용", value=True)
        if auto_gel:
            rec_system, rec_pct = _recommend_gel(mw_kda)
            gel_system = st.selectbox("Gel system", ["Tris-Glycine", "Tricine"], index=0 if rec_system == "Tris-Glycine" else 1)
            gel_pct = st.selectbox("Gel %", [8.0, 10.0, 12.0, 15.0, 16.5], index=[8.0,10.0,12.0,15.0,16.5].index(rec_pct))
            st.caption(f"추천: **{rec_system} / {rec_pct}%**")
        else:
            gel_system = st.selectbox("Gel system", ["Tris-Glycine", "Tricine"], index=0)
            gel_pct = st.selectbox("Gel %", [8.0, 10.0, 12.0, 15.0, 16.5], index=2)

        lane_count = st.slider("Lanes (visual)", min_value=2, max_value=6, value=3)

    # Practical warnings
    if mw_kda < 10 and gel_system != "Tricine":
        st.warning("10 kDa 미만이면 Tricine gel이 밴드 분리/식별에 유리합니다.")

    # Plot
    fig, ax = plt.subplots(figsize=(5.6, 7.6))
    ax.set_xlim(0, lane_count)
    ax.set_ylim(1, 0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"Virtual SDS-PAGE ({gel_system}, {gel_pct}%)")

    # lanes background
    for i in range(lane_count):
        ax.plot([i + 0.5, i + 0.5], [0.02, 0.98], linewidth=6, alpha=0.06)

    # ladder lane (0)
    ladder = ladders[ladder_name]
    for mw in ladder:
        y = _mw_to_migration(float(mw), gel_system, float(gel_pct))
        ax.plot([0.35, 0.65], [y, y], linewidth=2)
        ax.text(0.68, y, f"{mw}k", fontsize=8, va="center")

    # target lane (1)
    y_target = _mw_to_migration(float(mw_kda), gel_system, float(gel_pct))
    ax.plot([1.35, 1.65], [y_target, y_target], linewidth=3)
    ax.text(1.7, y_target, f"Target ~{mw_kda:.2f}k", fontsize=9, va="center")

    st.pyplot(fig)

    st.divider()
    st.subheader("로딩/마커 선택 근거 (실험용)")

    if mw_kda < 10:
        marker_tip = "저분자(≤10 kDa)는 Tricine gel + low range ladder(2–40k 또는 5–100k) 추천"
        loading_tip = "초기 확인은 2–5 µg 로딩부터(샘플 농도/염/이미다졸 상태 고려)"
    else:
        marker_tip = "타깃 주변 눈금이 촘촘한 ladder를 선택(5–100k 또는 10–250k 중)"
        loading_tip = "초기 확인은 1–3 µg 로딩부터(과로딩 시 smear 주의)"

    st.markdown(
        f"""
- **Gel 선택 근거:** MW {mw_kda:.2f} kDa 기준 {gel_system} / {gel_pct}%
- **Marker 추천:** {marker_tip}
- **로딩 추천:** {loading_tip}

**해석 포인트(빠른 CAPA)**
- FT에서 타깃이 강함 → binding 부족(레진량/이미다졸/pH/염/접촉시간)
- Wash에서 타깃이 강함 → wash imidazole 과다 또는 결합 약함
- Elution에서 없음 → imidazole/pH/레진 상태/His-tag 노출 확인
"""
    )


# ==============================
# Router
# ==============================
if page == "Protein MW":
    page_protein_mw()
elif page == "Ni-NTA Buffer (stock-based)":
    page_buffer_designer()
elif page == "SOP & Checklist":
    page_sop()
elif page == "Virtual SDS-PAGE":

    page_virtual_sds()
