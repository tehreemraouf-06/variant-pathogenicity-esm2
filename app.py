import json
import joblib
import numpy as np
import pandas as pd
import streamlit as st
import torch
from Bio.Align import substitution_matrices
from transformers import AutoTokenizer, EsmForMaskedLM

st.set_page_config(page_title="Variant Pathogenicity (ESM-2)", page_icon="🧬")

LETTERS = list("ACDEFGHIKLMNPQRSTVWY")
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2,
           I=166.7, L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8,
           Y=193.6, V=140.0)
CHG = {a: 0.0 for a in LETTERS}
CHG.update(K=1, R=1, H=0.1, D=-1, E=-1)
FEATURES = ["blosum", "d_hydro", "d_vol", "d_charge", "rel_pos", "esm_score", "esm_wt_logp"]


@st.cache_resource
def load_all():
    name = "facebook/esm2_t12_35M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    esm = EsmForMaskedLM.from_pretrained(name).eval()
    clf = joblib.load("model.joblib")
    seqs = json.load(open("sequences.json"))
    return tok, esm, clf, seqs


def position_logp(tok, esm, seq, pos, window=1000):
    if len(seq) > window:
        start = max(0, min(pos - 1 - window // 2, len(seq) - window))
        seq = seq[start:start + window]
        pos = pos - start
    enc = tok(seq, return_tensors="pt")
    ids = enc["input_ids"].clone()
    ids[0, pos] = tok.mask_token_id
    with torch.no_grad():
        logits = esm(input_ids=ids, attention_mask=enc["attention_mask"]).logits
    logp = torch.log_softmax(logits[0, pos], dim=-1)
    return {a: logp[tok.convert_tokens_to_ids(a)].item() for a in LETTERS}


st.title("🧬 Missense Variant Pathogenicity Predictor")
st.caption("ESM-2 (35M) + LightGBM, trained on ClinVar with gene-held-out validation. "
           "**Research prototype. Not for clinical use.**")

with st.spinner("Loading models (the first load takes a minute)..."):
    tok, esm, clf, seqs = load_all()

gene = st.selectbox("Gene", sorted(seqs))
seq = seqs[gene]
pos = st.number_input(f"Position (1 to {len(seq)})", min_value=1, max_value=len(seq), value=1, step=1)
wt = seq[pos - 1]
mut = st.selectbox("New amino acid", [a for a in LETTERS if a != wt])
st.write(f"Original amino acid at this position: **{wt}**")

if st.button("Predict"):
    lp = position_logp(tok, esm, seq, int(pos))
    B = substitution_matrices.load("BLOSUM62")
    row = pd.DataFrame([{
        "blosum": B[wt][mut],
        "d_hydro": KD[mut] - KD[wt],
        "d_vol": VOL[mut] - VOL[wt],
        "d_charge": CHG[mut] - CHG[wt],
        "rel_pos": pos / len(seq),
        "esm_score": lp[mut] - lp[wt],
        "esm_wt_logp": lp[wt],
    }])[FEATURES]
    prob = clf.predict_proba(row)[0, 1]

    st.metric("Predicted probability of pathogenic", f"{prob:.2f}")
    st.progress(float(prob))
    st.write(f"ESM-2 score: **{row['esm_score'].iloc[0]:.2f}** "
             "(more negative means a more surprising change)")
    st.info("This is a computational estimate for research and prioritization only. "
            "It is not a diagnosis.")

st.markdown("---")
st.caption("Limits: only genes seen in training are available, only missense variants are "
           "supported, and ClinVar labels are imperfect.")
