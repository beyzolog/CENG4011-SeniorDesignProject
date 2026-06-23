import io
import json
import math
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

# ── Page config ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Drug Discovery Portal | MYC & CTNNB1",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ─────────────────────────────────────────────────────────────────
st.markdown(
    """
    <style>
    /* ---- Base & fonts ---- */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', sans-serif;
    }

    /* ---- Sidebar ---- */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #0d1117 0%, #161b22 100%);
        border-right: 1px solid #30363d;
    }
    [data-testid="stSidebar"] * {
        color: #e6edf3 !important;
    }

    /* ---- Sidebar collapse/expand toggle button — all known selectors ---- */
    [data-testid="collapsedControl"],
    [data-testid="stSidebarCollapsedControl"],
    section[data-testid="stSidebarCollapsedControl"] {
        background-color: #21262d !important;
        border: 2px solid #388bfd !important;
        border-radius: 0 10px 10px 0 !important;
        visibility: visible !important;
        opacity: 1 !important;
        display: flex !important;
        z-index: 9999 !important;
    }
    [data-testid="collapsedControl"] button,
    [data-testid="stSidebarCollapsedControl"] button {
        background-color: #21262d !important;
        color: #e6edf3 !important;
        visibility: visible !important;
        opacity: 1 !important;
    }
    [data-testid="collapsedControl"] svg,
    [data-testid="stSidebarCollapsedControl"] svg,
    [data-testid="collapsedControl"] svg path,
    [data-testid="stSidebarCollapsedControl"] svg path {
        fill: #e6edf3 !important;
        stroke: #e6edf3 !important;
        visibility: visible !important;
        opacity: 1 !important;
    }
    [data-testid="collapsedControl"]:hover,
    [data-testid="stSidebarCollapsedControl"]:hover {
        background-color: #1f6feb !important;
        border-color: #58a6ff !important;
    }

    /* ---- Main background ---- */
    .stApp {
        background-color: #0d1117;
        color: #e6edf3;
    }

    /* ---- Metric cards ---- */
    .metric-card {
        background: linear-gradient(135deg, #161b22 0%, #21262d 100%);
        border: 1px solid #30363d;
        border-radius: 12px;
        padding: 20px 24px;
        text-align: center;
        box-shadow: 0 4px 16px rgba(0,0,0,0.4);
        transition: border-color 0.2s;
    }
    .metric-card:hover { border-color: #58a6ff; }
    .metric-label {
        font-size: 0.78rem;
        font-weight: 600;
        letter-spacing: 0.08em;
        text-transform: uppercase;
        color: #8b949e;
        margin-bottom: 6px;
    }
    .metric-value {
        font-size: 2rem;
        font-weight: 700;
        color: #58a6ff;
        line-height: 1.1;
    }
    .metric-sub {
        font-size: 0.75rem;
        color: #8b949e;
        margin-top: 4px;
    }

    /* ---- Section headers ---- */
    .section-header {
        font-size: 1.1rem;
        font-weight: 600;
        color: #58a6ff;
        border-left: 3px solid #1f6feb;
        padding-left: 10px;
        margin: 24px 0 12px 0;
        letter-spacing: 0.02em;
    }

    /* ---- Lipinski property card ---- */
    .prop-card {
        background: #161b22;
        border: 1px solid #30363d;
        border-radius: 8px;
        padding: 12px 16px;
        margin-bottom: 8px;
        display: flex;
        justify-content: space-between;
        align-items: center;
    }
    .prop-name  { color: #8b949e; font-size: 0.85rem; }
    .prop-value { color: #e6edf3; font-weight: 600; font-family: 'JetBrains Mono', monospace; }
    .prop-pass  { color: #3fb950; font-size: 0.75rem; font-weight: 600; }
    .prop-fail  { color: #f85149; font-size: 0.75rem; font-weight: 600; }

    /* ---- TÜBİTAK placeholder ---- */
    .tubitak-box {
        background: linear-gradient(135deg, #1a2332 0%, #1f2d3d 100%);
        border: 2px dashed #30363d;
        border-radius: 10px;
        padding: 16px;
        text-align: center;
        margin: 12px 0;
    }
    .tubitak-title {
        font-size: 0.65rem;
        font-weight: 700;
        letter-spacing: 0.15em;
        color: #8b949e;
        text-transform: uppercase;
    }
    .tubitak-main {
        font-size: 1.1rem;
        font-weight: 700;
        color: #58a6ff;
        margin: 4px 0;
    }
    .tubitak-sub {
        font-size: 0.7rem;
        color: #8b949e;
    }

    /* ---- Footer ---- */
    .footer-bar {
        background: linear-gradient(90deg, #1a2332 0%, #1f2d3d 100%);
        border: 1px solid #30363d;
        border-radius: 10px;
        padding: 16px 24px;
        margin-top: 40px;
        text-align: center;
        color: #8b949e;
        font-size: 0.82rem;
    }
    .footer-icon { font-size: 1.1rem; margin-right: 6px; }

    /* ---- Dataframe tweaks ---- */
    [data-testid="stDataFrame"] { border-radius: 10px; overflow: hidden; }

    /* ---- Gene badge ---- */
    .gene-badge {
        display: inline-block;
        background: linear-gradient(90deg, #1f6feb, #388bfd);
        color: white;
        border-radius: 20px;
        padding: 4px 14px;
        font-size: 0.78rem;
        font-weight: 600;
        letter-spacing: 0.04em;
        margin-left: 10px;
        vertical-align: middle;
    }

    /* ---- Page title ---- */
    .portal-title {
        font-size: 1.9rem;
        font-weight: 700;
        color: #e6edf3;
        line-height: 1.2;
    }
    .portal-subtitle {
        color: #8b949e;
        font-size: 0.9rem;
        margin-top: 4px;
    }

    /* ---- PDB info card ---- */
    .pdb-card {
        background: #161b22;
        border: 1px solid #30363d;
        border-radius: 8px;
        padding: 10px 16px;
        font-size: 0.82rem;
        color: #8b949e;
    }
    .pdb-code {
        font-family: 'JetBrains Mono', monospace;
        color: #58a6ff;
        font-weight: 600;
    }

    /* ---- Hide Streamlit default elements ---- */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    /* header {visibility: hidden;}  <-- Bu satırı sildik veya yorum satırı yaptık */

    /* Sadece üst barın arka planını şeffaf yapıp butonu koruyalım */
    header {
    background-color: rgba(0,0,0,0) !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ── JavaScript: keep sidebar toggle button always visible ──────────────────────
# CSS alone is unreliable because Streamlit re-renders and wipes inline styles.
# A MutationObserver watches for DOM changes and re-applies styles immediately.
import streamlit.components.v1 as components

components.html(
    """
    <script>
    (function() {
        function styleToggle() {
            // Try every known selector across Streamlit versions
            var selectors = [
                '[data-testid="collapsedControl"]',
                '[data-testid="stSidebarCollapsedControl"]',
                'section[data-testid="stSidebarCollapsedControl"]'
            ];
            selectors.forEach(function(sel) {
                var els = window.parent.document.querySelectorAll(sel);
                els.forEach(function(el) {
                    el.style.setProperty('background-color', '#21262d', 'important');
                    el.style.setProperty('border', '2px solid #388bfd', 'important');
                    el.style.setProperty('border-radius', '0 10px 10px 0', 'important');
                    el.style.setProperty('visibility', 'visible', 'important');
                    el.style.setProperty('opacity', '1', 'important');
                    el.style.setProperty('display', 'flex', 'important');
                    el.style.setProperty('z-index', '9999', 'important');

                    // Inner button
                    var btns = el.querySelectorAll('button');
                    btns.forEach(function(btn) {
                        btn.style.setProperty('background-color', '#21262d', 'important');
                        btn.style.setProperty('color', '#e6edf3', 'important');
                        btn.style.setProperty('visibility', 'visible', 'important');
                        btn.style.setProperty('opacity', '1', 'important');
                    });

                    // SVG icons
                    var svgs = el.querySelectorAll('svg, svg path, svg polyline, svg line');
                    svgs.forEach(function(s) {
                        s.style.setProperty('fill', '#e6edf3', 'important');
                        s.style.setProperty('stroke', '#e6edf3', 'important');
                        s.style.setProperty('visibility', 'visible', 'important');
                        s.style.setProperty('opacity', '1', 'important');
                    });
                });
            });
        }

        // Run immediately
        styleToggle();

        // Run on every DOM mutation in the parent frame
        var observer = new MutationObserver(styleToggle);
        observer.observe(window.parent.document.body, {
            childList: true,
            subtree: true,
            attributes: true
        });

        // Fallback interval every 300ms
        setInterval(styleToggle, 300);
    })();
    </script>
    """,
    height=0,
    scrolling=False,
)

# ── Constants ──────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

GENE_CONFIG = {
    "MYC (Transcription Factor)": {
        "csv": os.path.join(BASE_DIR, "myc_exp06_final_1000.csv"),
        "pdb": "6G6K",
        "color": "#58a6ff",
        "pathway": "Transcription Factor",
        "model_type": "Extra Trees",
        "experiment": "exp_06",
        "auc_cv": 0.9084,
        "auc_train": 0.9443,
        "auc_test": 0.9082,
        "description": "MYC oncoprotein — master regulator of cell proliferation",
    },
    "CTNNB1 (Wnt Signaling)": {
        "csv": os.path.join(BASE_DIR, "ctnnb1_exp05_final_1000.csv"),
        "pdb": "1JPW",
        "color": "#3fb950",
        "pathway": "Wnt Signaling",
        "model_type": "Random Forest",
        "experiment": "exp_05",
        "auc_cv": 0.8723,
        "auc_train": 0.8968,
        "auc_test": 0.8989,
        "description": "β-catenin (CTNNB1) — key mediator of Wnt/β-catenin pathway",
    },
}

# ── Helper: load CSV ───────────────────────────────────────────────────────────
@st.cache_data(show_spinner=False)
def load_data(path: str) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(path)
        df = df.sort_values("prediction_score", ascending=False).reset_index(drop=True)
        return df
    except FileNotFoundError:
        return None
    except Exception as e:
        st.error(f"Veri yüklenirken beklenmeyen hata: {e}")
        return None


# ── Helper: RDKit 2D structure ─────────────────────────────────────────────────
def draw_molecule_2d(smiles: str):
    """Returns a PIL Image of the 2D structure or None on failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        img = Draw.MolToImage(mol, size=(400, 300))
        return img
    except ImportError:
        return None
    except Exception:
        return None


# ── Helper: Lipinski properties ────────────────────────────────────────────────
def compute_lipinski(smiles: str) -> dict | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        props = {
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "RotBonds": Descriptors.NumRotatableBonds(mol),
        }
        violations = sum([
            props["MW"] > 500,
            props["LogP"] > 5,
            props["HBD"] > 5,
            props["HBA"] > 10,
        ])
        props["violations"] = violations
        props["Lipinski_Pass"] = (violations == 0)
        props["bRo5_Pass"] = (violations <= 2)
        return props
    except ImportError:
        return None
    except Exception:
        return None


# ── Helper: render Lipinski card ───────────────────────────────────────────────
def render_lipinski_card(props: dict):
    rules = {
        "MW": ("<= 500 Da", props["MW"] <= 500),
        "LogP": ("<= 5", props["LogP"] <= 5),
        "HBD": ("<= 5", props["HBD"] <= 5),
        "HBA": ("<= 10", props["HBA"] <= 10),
        "TPSA": ("<= 140 Å²", props["TPSA"] <= 140),
        "RotBonds": ("<= 10", props["RotBonds"] <= 10),
    }
    labels = {
        "MW": "Molecular Weight",
        "LogP": "LogP (lipophilicity)",
        "HBD": "H-Bond Donors",
        "HBA": "H-Bond Acceptors",
        "TPSA": "Topological PSA",
        "RotBonds": "Rotatable Bonds",
    }
    units = {
        "MW": " Da",
        "LogP": "",
        "HBD": "",
        "HBA": "",
        "TPSA": " Å²",
        "RotBonds": "",
    }
    html = ""
    for key, (rule_text, passes) in rules.items():
        status_cls = "prop-pass" if passes else "prop-fail"
        status_txt = f"✓ {rule_text}" if passes else f"✗ {rule_text}"
        html += f"""
        <div class="prop-card">
            <span class="prop-name">{labels[key]}</span>
            <span>
                <span class="prop-value">{props[key]}{units[key]}</span>&nbsp;&nbsp;
                <span class="{status_cls}">{status_txt}</span>
            </span>
        </div>"""
    st.markdown(html, unsafe_allow_html=True)

    v = props["violations"]

    ro5_color = "#3fb950" if props["Lipinski_Pass"] else "#f85149"
    ro5_text  = "PASSES Lipinski Ro5 (0 violations)" if props["Lipinski_Pass"] \
                else f"FAILS Lipinski Ro5 ({v} violation{'s' if v != 1 else ''})"

    bro5_color = "#3fb950" if props["bRo5_Pass"] else "#f85149"
    bro5_text  = f"PASSES bRo5 (≤ 2 violations allowed)" if props["bRo5_Pass"] \
                 else f"FAILS bRo5 ({v} violations — exceeds relaxed threshold)"

    st.markdown(
        f'<div style="margin-top:10px; border-radius:8px; overflow:hidden; '
        f'border:1px solid #30363d;">'
        f'<div style="padding:7px 12px; background:#161b22; border-bottom:1px solid #30363d; '
        f'color:{ro5_color}; font-weight:700; font-size:0.82rem; text-align:center;">'
        f'{ro5_text}</div>'
        f'<div style="padding:7px 12px; background:#0d1117; '
        f'color:{bro5_color}; font-size:0.78rem; text-align:center;">'
        f'{bro5_text}</div>'
        f'</div>',
        unsafe_allow_html=True,
    )


# ══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown(
        """
        <div style="text-align:center; padding: 10px 0 6px 0;">
            <div style="font-size:1.15rem; font-weight:700; color:#e6edf3; line-height:1.3;">
                Non-Globular<br>
                <span style="color:#58a6ff;">Oncoprotein Drug Discovery</span>
            </div>
            <div style="font-size:0.72rem; color:#8b949e; margin-top:4px; letter-spacing: 0.05em;">
                VIRTUAL SCREENING PORTAL
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.divider()

    # # TÜBİTAK placeholder
    # st.markdown(
    #     """
    #     <div class="tubitak-box">
    #         <div class="tubitak-title">Supported by</div>
    #         <div class="tubitak-main">🇹🇷 TÜBİTAK</div>
    #         <div style="font-size:0.85rem; font-weight:700; color:#58a6ff; margin:2px 0;">2209-A</div>
    #         <div class="tubitak-sub">Undergraduate Research Grant</div>
    #     </div>
    #     """,
    #     unsafe_allow_html=True,
    # )

    # st.divider()

    # Gene selection
    st.markdown(
        '<div style="font-size:0.72rem; font-weight:600; letter-spacing:0.1em; '
        'color:#8b949e; text-transform:uppercase; margin-bottom:8px;">Target Selection</div>',
        unsafe_allow_html=True,
    )
    selected_gene = st.radio(
        label="Select Target Gene",
        options=list(GENE_CONFIG.keys()),
        label_visibility="collapsed",
    )

    cfg = GENE_CONFIG[selected_gene]

    st.divider()

    # Score filter
    st.markdown(
        '<div style="font-size:0.72rem; font-weight:600; letter-spacing:0.1em; '
        'color:#8b949e; text-transform:uppercase; margin-bottom:8px;">Score Filter</div>',
        unsafe_allow_html=True,
    )
    min_score = st.slider(
        "Minimum Prediction Score",
        min_value=0.80,
        max_value=1.00,
        value=0.85,
        step=0.005,
        format="%.3f",
        label_visibility="collapsed",
    )

    st.divider()

    # Methodology summary
    st.markdown(
        """
        <div style="font-size:0.72rem; font-weight:600; letter-spacing:0.1em;
            color:#8b949e; text-transform:uppercase; margin-bottom:10px;">
            Methodology
        </div>
        <div style="font-size:0.78rem; color:#c9d1d9; line-height:1.7;">
            <b style="color:#58a6ff;">1.</b> PubChem screening: <b>119M</b> compounds<br>
            <b style="color:#58a6ff;">2.</b> Morgan fingerprints (r=2, 2048-bit)<br>
            <b style="color:#58a6ff;">3.</b> MYC: <b>Extra Trees</b> (exp_06)<br>
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CTNNB1: <b>Random Forest</b> (exp_05)<br>
            <b style="color:#58a6ff;">4.</b> Lipinski + Toxicity filters applied<br>
            <b style="color:#58a6ff;">5.</b> Top 1000 candidates retained
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.divider()

    st.markdown(
    '<div style="font-size:0.68rem; color:#8b949e; text-align:center; font-family:monospace;">'
    "Muğla Sıtkı Koçman University<br>"
    "Süzek Lab Bioinformatics Infrastructure<br>"
    "</div>",
    unsafe_allow_html=True,
)

# ══════════════════════════════════════════════════════════════════════════════
# MAIN PANEL
# ══════════════════════════════════════════════════════════════════════════════

# ── Page title ─────────────────────────────────────────────────────────────────
gene_short = "MYC" if "MYC" in selected_gene else "CTNNB1"
st.markdown(
    f"""
    <div style="padding: 4px 0 20px 0;">
        <div class="portal-title">
            Drug Discovery Portal
            <span class="gene-badge">{gene_short}</span>
        </div>
        <div class="portal-subtitle">{cfg['description']}</div>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Gene Information Card ─────────────────────────────────────────────────────
GENE_INFO = {
    "MYC (Transcription Factor)": {
        "symbol": "MYC",
        "full_name": "MYC Proto-Oncogene, bHLH Transcription Factor",
        "ensembl_id": "ENSG00000136997",
        "description": "MYC proto-oncogene, bHLH transcription factor",
        "aliases": ["MRTL", "MYCC", "bHLHe39", "c-Myc"],
        "summary": (
            "This gene is a proto-oncogene and encodes a nuclear phosphoprotein that plays a "
            "role in cell cycle progression, apoptosis and cellular transformation. The encoded "
            "protein forms a heterodimer with the related transcription factor MAX. This complex "
            "binds to the E box DNA consensus sequence and regulates the transcription of specific "
            "target genes. Amplification of this gene is frequently observed in numerous human "
            "cancers. Translocations involving this gene are associated with Burkitt lymphoma and "
            "multiple myeloma in human patients. There is evidence to show that translation "
            "initiates both from an upstream, in-frame non-AUG (CUG) and a downstream AUG start "
            "site, resulting in the production of two isoforms with distinct N-termini. "
            "[provided by RefSeq, Aug 2017]"
        ),
        "links": {
            "NCBI": "https://www.ncbi.nlm.nih.gov/gene/4609",
            "Ensembl": "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000136997",
            "EBI": "https://www.ebi.ac.uk/geneontology/ENSG00000136997",
            "GeneCards": "https://www.genecards.org/cgi-bin/carddisp.pl?gene=MYC",
            "OMIM": "https://www.omim.org/entry/190080",
            "COSMIC": "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=MYC",
            "HPA": "https://www.proteinatlas.org/ENSG00000136997-MYC",
            "DrugBank": "https://go.drugbank.com/genes/BE0000533",
            "Xena": "https://xenabrowser.net/datapages/?gene=MYC",
            "cBioPortal": "https://www.cbioportal.org/results?gene_list=MYC",
        },
    },
    "CTNNB1 (Wnt Signaling)": {
        "symbol": "CTNNB1",
        "full_name": "Catenin Beta 1",
        "ensembl_id": "ENSG00000168036",
        "description": "catenin beta 1",
        "aliases": ["CTNNB", "EVR7", "MRD19", "NEDSDV", "armadillo"],
        "summary": (
            "The protein encoded by this gene is part of a complex of proteins that constitute "
            "adherens junctions (AJs). AJs are necessary for the creation and maintenance of "
            "epithelial cell layers by regulating cell growth and adhesion between cells. The "
            "encoded protein also anchors the actin cytoskeleton and may be responsible for "
            "transmitting the contact inhibition signal that causes cells to stop dividing once "
            "the epithelial sheet is complete. Finally, this protein binds to the product of the "
            "APC gene, which is mutated in adenomatous polyposis of the colon. Mutations in this "
            "gene are a cause of colorectal cancer (CRC), pilomatrixoma (PTR), medulloblastoma "
            "(MDB), and ovarian cancer. Alternative splicing results in multiple transcript "
            "variants. [provided by RefSeq, Aug 2016]"
        ),
        "links": {
            "NCBI": "https://www.ncbi.nlm.nih.gov/gene/1499",
            "Ensembl": "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000168036",
            "EBI": "https://www.ebi.ac.uk/geneontology/ENSG00000168036",
            "GeneCards": "https://www.genecards.org/cgi-bin/carddisp.pl?gene=CTNNB1",
            "OMIM": "https://www.omim.org/entry/116806",
            "COSMIC": "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=CTNNB1",
            "HPA": "https://www.proteinatlas.org/ENSG00000168036-CTNNB1",
            "DrugBank": "https://go.drugbank.com/genes/BE0000551",
            "Xena": "https://xenabrowser.net/datapages/?gene=CTNNB1",
            "cBioPortal": "https://www.cbioportal.org/results?gene_list=CTNNB1",
        },
    },
}

gi = GENE_INFO[selected_gene]

st.markdown('<div class="section-header">Gene Information</div>', unsafe_allow_html=True)

with st.container():
    gi_left, gi_right = st.columns([3, 2], gap="large")

    with gi_left:
        # Title row
        aliases_html = " &nbsp;·&nbsp; ".join(
            f'<span style="font-family:\'JetBrains Mono\',monospace; '
            f'font-size:0.72rem; color:#58a6ff; background:#1f3a5f; '
            f'border-radius:4px; padding:1px 6px;">{a}</span>'
            for a in gi["aliases"]
        )
        st.markdown(
            f"""
            <div style="background:#161b22; border:1px solid #30363d; border-radius:12px;
                padding:20px 24px 16px 24px;">
                <div style="font-size:1.5rem; font-weight:800; color:#e6edf3; line-height:1.1;">
                    {gi['symbol']}
                    <span style="font-size:0.9rem; font-weight:400; color:#8b949e; margin-left:8px;">
                        {gi['full_name']}
                    </span>
                </div>
                <div style="margin-top:8px; display:flex; flex-wrap:wrap; align-items:center; gap:8px;">
                    <span style="font-size:0.75rem; color:#8b949e;">Ensembl:</span>
                    <span style="font-family:'JetBrains Mono',monospace; font-size:0.75rem;
                        color:#3fb950;">{gi['ensembl_id']}</span>
                    <span style="color:#30363d;">|</span>
                    <span style="font-size:0.75rem; color:#8b949e;">Aliases:</span>
                    {aliases_html}
                </div>
                <div style="margin-top:14px; font-size:0.8rem; color:#8b949e;
                    font-style:italic; letter-spacing:0.02em;">
                    {gi['description']}
                </div>
                <div style="margin-top:12px; font-size:0.82rem; color:#c9d1d9;
                    line-height:1.75; text-align:justify;">
                    {gi['summary']}
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

    with gi_right:
        link_buttons = "".join(
            f'<a href="{url}" target="_blank" style="'
            f"display:inline-block; margin:4px 4px 0 0; padding:5px 12px; "
            f"background:#21262d; border:1px solid #30363d; border-radius:6px; "
            f"color:#58a6ff; font-size:0.75rem; font-weight:600; text-decoration:none; "
            f'letter-spacing:0.03em; transition:border-color 0.2s;"'
            f' onmouseover="this.style.borderColor=\'#58a6ff\'"'
            f' onmouseout="this.style.borderColor=\'#30363d\'"'
            f">{name} ↗</a>"
            for name, url in gi["links"].items()
        )
        st.markdown(
            f"""
            <div style="background:#161b22; border:1px solid #30363d; border-radius:12px;
                padding:20px 24px; height:100%;">
                <div style="font-size:0.72rem; font-weight:700; letter-spacing:0.1em;
                    text-transform:uppercase; color:#8b949e; margin-bottom:14px;">
                    External Databases
                </div>
                <div style="line-height:2.2;">
                    {link_buttons}
                </div>
                <div style="margin-top:20px; padding-top:14px; border-top:1px solid #30363d;">
                    <div style="font-size:0.72rem; font-weight:700; letter-spacing:0.1em;
                        text-transform:uppercase; color:#8b949e; margin-bottom:10px;">
                        Target Classification
                    </div>
                    <div style="font-size:0.8rem; color:#c9d1d9; line-height:1.9;">
                        <span style="color:#8b949e;">Pathway:</span>&nbsp;
                        <b>{cfg['pathway']}</b><br>
                        <span style="color:#8b949e;">Experiment:</span>&nbsp;
                        <b>{cfg['experiment']}</b><br>
                        <span style="color:#8b949e;">Model:</span>&nbsp;
                        <b>{cfg['model_type']}</b><br>
                        <span style="color:#8b949e;">PDB Structure:</span>&nbsp;
                        <a href="https://www.rcsb.org/structure/{cfg['pdb']}"
                           target="_blank"
                           style="color:#58a6ff; font-family:'JetBrains Mono',monospace;
                                  font-weight:600; text-decoration:none;">
                           {cfg['pdb']} ↗
                        </a>
                    </div>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

st.markdown("<br>", unsafe_allow_html=True)

# ── Load data ──────────────────────────────────────────────────────────────────
df_full = load_data(cfg["csv"])

if df_full is None:
    st.error(
        f"CSV dosyası bulunamadı: `{cfg['csv']}`\n\n"
        "Lütfen dosya yolunu kontrol edin."
    )
    st.stop()

df_filtered = df_full[df_full["prediction_score"] >= min_score].reset_index(drop=True)
df_top100 = df_full.head(100)

# ── Data & Features section ────────────────────────────────────────────────────
st.markdown('<div class="section-header">Data & Features</div>', unsafe_allow_html=True)

DATA_META = {
    "MYC (Transcription Factor)": {
        "active": 483, "inactive": 177, "total": 660,
        "train": 495, "test": 165,
        "fp_full": 2048, "fp_selected": 150,
    },
    "CTNNB1 (Wnt Signaling)": {
        "active": 516, "inactive": 199, "total": 715,
        "train": 536, "test": 179,
        "fp_full": 2048, "fp_selected": 30,
    },
}
dm = DATA_META[selected_gene]

data_col, feat_col = st.columns([1, 1], gap="large")

with data_col:
    active_pct   = dm["active"]   / dm["total"] * 100
    inactive_pct = dm["inactive"] / dm["total"] * 100
    imbalance    = round(dm["active"] / dm["inactive"], 2)

    st.markdown(
        "**Training Dataset Composition**",
    )

    tbl_data = {
        "Metric":  ["Active compounds", "Inactive compounds", "Total samples",
                    "Train / Test split", "Class imbalance ratio"],
        "Count":   [str(dm["active"]), str(dm["inactive"]), str(dm["total"]),
                    f"{dm['train']} / {dm['test']}", f"{imbalance}x"],
        "Share":   [f"{active_pct:.1f}%", f"{inactive_pct:.1f}%", "—",
                    "75 / 25%", "balanced weights"],
    }
    st.dataframe(
        pd.DataFrame(tbl_data),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown(f"Active **{active_pct:.1f}%** &nbsp;|&nbsp; Inactive **{inactive_pct:.1f}%**",
                unsafe_allow_html=True)
    st.progress(active_pct / 100)

with feat_col:
    reduction_pct = (1 - dm["fp_selected"] / dm["fp_full"]) * 100

    st.markdown("**Feature Engineering Pipeline**")

    pipeline_rows = {
        "Step": ["1 — Morgan Fingerprints",
                 "2 — Feature Importance",
                 "3 — Denoising (Top-k)",
                 "4 — Final Feature Set"],
        "Detail": [
            f"Radius=2, full bit-vector: {dm['fp_full']:,} bits",
            "Tree-based impurity importance (Gini) on full feature set",
            "Low-importance bits discarded; model retrained on signal-only features",
            f"{dm['fp_selected']} bits selected — {reduction_pct:.1f}% noise removed",
        ],
    }
    st.dataframe(
        pd.DataFrame(pipeline_rows),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown(
        f"Signal bits retained: **{dm['fp_selected']}** / {dm['fp_full']} "
        f"&nbsp;·&nbsp; Model: **{cfg['model_type']}**",
        unsafe_allow_html=True,
    )
    st.progress((dm["fp_selected"] / dm["fp_full"]))

st.markdown("<br>", unsafe_allow_html=True)

# ── Metric cards ───────────────────────────────────────────────────────────────
st.markdown('<div class="section-header">Overview Metrics</div>', unsafe_allow_html=True)

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.markdown(
        """
        <div class="metric-card">
            <div class="metric-label">Compounds Screened</div>
            <div class="metric-value">119M</div>
            <div class="metric-sub">PubChem database</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col2:
    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Model AUC Score (Test)</div>
            <div class="metric-value">{cfg['auc_test']:.4f}</div>
            <div class="metric-sub">{cfg['model_type']} · {cfg['experiment']}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col3:
    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Total Hits Retained</div>
            <div class="metric-value">{len(df_full):,}</div>
            <div class="metric-sub">After all filters</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col4:
    avg_score = df_full["prediction_score"].mean()
    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Avg Prediction Score</div>
            <div class="metric-value">{avg_score:.4f}</div>
            <div class="metric-sub">Top 1000 candidates</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

# ── AUC breakdown row ──────────────────────────────────────────────────────────
st.markdown(
    f"""
    <div style="display:flex; gap:12px; margin-top:12px; margin-bottom:4px;">
        <div style="flex:1; background:#161b22; border:1px solid #30363d; border-radius:8px;
            padding:10px 16px; display:flex; justify-content:space-between; align-items:center;">
            <span style="font-size:0.75rem; color:#8b949e; font-weight:600; text-transform:uppercase; letter-spacing:0.08em;">
                Model
            </span>
            <span style="font-size:0.88rem; color:#e6edf3; font-weight:600;">
                {cfg['model_type']}
            </span>
        </div>
        <div style="flex:1; background:#161b22; border:1px solid #30363d; border-radius:8px;
            padding:10px 16px; display:flex; justify-content:space-between; align-items:center;">
            <span style="font-size:0.75rem; color:#8b949e; font-weight:600; text-transform:uppercase; letter-spacing:0.08em;">
                CV AUC
            </span>
            <span style="font-size:0.88rem; color:#58a6ff; font-weight:700; font-family:'JetBrains Mono',monospace;">
                {cfg['auc_cv']:.4f}
            </span>
        </div>
        <div style="flex:1; background:#161b22; border:1px solid #30363d; border-radius:8px;
            padding:10px 16px; display:flex; justify-content:space-between; align-items:center;">
            <span style="font-size:0.75rem; color:#8b949e; font-weight:600; text-transform:uppercase; letter-spacing:0.08em;">
                Train AUC
            </span>
            <span style="font-size:0.88rem; color:#58a6ff; font-weight:700; font-family:'JetBrains Mono',monospace;">
                {cfg['auc_train']:.4f}
            </span>
        </div>
        <div style="flex:1; background:#161b22; border:1px solid #388bfd; border-radius:8px;
            padding:10px 16px; display:flex; justify-content:space-between; align-items:center;">
            <span style="font-size:0.75rem; color:#8b949e; font-weight:600; text-transform:uppercase; letter-spacing:0.08em;">
                Test AUC ✓
            </span>
            <span style="font-size:0.88rem; color:#3fb950; font-weight:700; font-family:'JetBrains Mono',monospace;">
                {cfg['auc_test']:.4f}
            </span>
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("<br>", unsafe_allow_html=True)

# ── 3D Protein structure & Histogram (side by side) ───────────────────────────
left_col, right_col = st.columns([1, 1], gap="large")

with left_col:
    st.markdown('<div class="section-header">3D Protein Structure</div>', unsafe_allow_html=True)
    st.markdown(
        f"""
        <div class="pdb-card">
            Target: <b style="color:#e6edf3;">{gene_short}</b> &nbsp;|&nbsp;
            PDB ID: <span class="pdb-code">{cfg['pdb']}</span> &nbsp;|&nbsp;
            Pathway: <b style="color:#e6edf3;">{cfg['pathway']}</b>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.markdown("<br>", unsafe_allow_html=True)

    try:
        import py3Dmol
        from stmol import showmol

        with st.spinner(f"Loading PDB structure {cfg['pdb']}..."):
            viewer = py3Dmol.view(query=f"pdb:{cfg['pdb']}", width=550, height=400)
            viewer.setStyle({"cartoon": {"color": "spectrum"}})
            viewer.addSurface(
                py3Dmol.VDW,
                {"opacity": 0.12, "color": "white"},
                {"hetflag": False},
            )
            viewer.setBackgroundColor("#0d1117")
            viewer.zoomTo()
            showmol(viewer, height=420, width=580)

    except ImportError:
        st.markdown(
            f"""
            <div style="background:#161b22; border:1px dashed #30363d; border-radius:10px;
                padding:40px; text-align:center; color:#8b949e;">
                <div style="font-size:2rem;">🧬</div>
                <div style="font-weight:600; color:#e6edf3; margin:8px 0;">PDB: {cfg['pdb']}</div>
                <div style="font-size:0.82rem; line-height:1.6;">
                    3D viewer requires <code>stmol</code> and <code>py3Dmol</code>.<br>
                    Install with: <code>pip install stmol py3Dmol</code><br><br>
                    <a href="https://www.rcsb.org/structure/{cfg['pdb']}"
                       style="color:#58a6ff;" target="_blank">
                       View on RCSB PDB →
                    </a>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    except Exception as e:
        st.warning(f"3D yapı yüklenemedi: {e}")

with right_col:
    st.markdown(
        '<div class="section-header">Prediction Score Distribution</div>',
        unsafe_allow_html=True,
    )

    fig_hist = px.histogram(
        df_full,
        x="prediction_score",
        nbins=40,
        title=f"{gene_short} — Score Distribution (n={len(df_full):,})",
        labels={"prediction_score": "Prediction Score", "count": "Count"},
        color_discrete_sequence=[cfg["color"]],
        template="plotly_dark",
    )
    fig_hist.update_layout(
        plot_bgcolor="#161b22",
        paper_bgcolor="#161b22",
        font={"family": "Inter", "color": "#e6edf3"},
        title_font_size=13,
        bargap=0.08,
        xaxis=dict(gridcolor="#30363d", linecolor="#30363d"),
        yaxis=dict(gridcolor="#30363d", linecolor="#30363d"),
        margin=dict(l=20, r=20, t=50, b=20),
    )
    if min_score > df_full["prediction_score"].min():
        fig_hist.add_vline(
            x=min_score,
            line_dash="dash",
            line_color="#f85149",
            annotation_text=f"Filter: {min_score:.3f}",
            annotation_font_color="#f85149",
        )
    st.plotly_chart(fig_hist, use_container_width=True)

st.divider()

# ── Top Hits table ─────────────────────────────────────────────────────────────
st.markdown(
    '<div class="section-header">Top 100 Drug Candidates</div>',
    unsafe_allow_html=True,
)

display_df = df_filtered.head(100).copy()
display_df.index = range(1, len(display_df) + 1)

info_col, _ = st.columns([3, 1])
with info_col:
    st.markdown(
        f'<div style="font-size:0.8rem; color:#8b949e; margin-bottom:8px;">'
        f"Showing top {len(display_df)} candidates with score ≥ {min_score:.3f} "
        f"(filtered from {len(df_full):,} total hits)"
        f"</div>",
        unsafe_allow_html=True,
    )

st.dataframe(
    display_df[["CID", "prediction_score", "Canonical_SMILES"]].rename(
        columns={
            "CID": "PubChem CID",
            "prediction_score": "Prediction Score",
            "Canonical_SMILES": "SMILES",
        }
    ),
    use_container_width=True,
    height=340,
    column_config={
        "Prediction Score": st.column_config.ProgressColumn(
            "Prediction Score",
            min_value=0.0,
            max_value=1.0,
            format="%.4f",
        )
    },
)

st.divider()

# ── Molecule detail panel ──────────────────────────────────────────────────────
st.markdown(
    '<div class="section-header">Molecule Inspector</div>', unsafe_allow_html=True
)

if display_df.empty:
    st.info("No molecules match the current score filter. Lower the minimum score in the sidebar.")
else:
    cid_options = display_df["CID"].tolist()
    selected_cid = st.selectbox(
        "Select a PubChem CID to inspect:",
        options=cid_options,
        format_func=lambda c: f"CID {c}",
    )

    mol_row = display_df[display_df["CID"] == selected_cid].iloc[0]
    smiles = mol_row["Canonical_SMILES"]
    score = mol_row["prediction_score"]

    mol_col, prop_col = st.columns([1, 1], gap="large")

    with mol_col:
        st.markdown(
            f'<div style="font-size:0.8rem; color:#8b949e; margin-bottom:8px;">'
            f"CID <b style='color:#58a6ff;'>{selected_cid}</b> &nbsp;|&nbsp; "
            f"Score: <b style='color:#3fb950;'>{score:.4f}</b>"
            f"</div>",
            unsafe_allow_html=True,
        )

        img = draw_molecule_2d(smiles)
        if img is not None:
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            buf.seek(0)
            st.image(buf, caption=f"2D Structure — CID {selected_cid}", use_container_width=True)
        else:
            st.markdown(
                """
                <div style="background:#161b22; border:1px dashed #30363d; border-radius:8px;
                    padding:30px; text-align:center; color:#8b949e;">
                    2D structure unavailable.<br>
                    Install rdkit: <code>pip install rdkit-pypi</code>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown(
            f'<div style="font-size:0.7rem; color:#484f58; margin-top:8px; '
            f'word-break:break-all; font-family:JetBrains Mono, monospace;">'
            f"SMILES: {smiles}</div>",
            unsafe_allow_html=True,
        )

    with prop_col:
        st.markdown(
            '<div style="font-size:0.85rem; font-weight:600; color:#e6edf3; '
            'margin-bottom:12px;">Molecular Properties</div>',
            unsafe_allow_html=True,
        )

        if "CTNNB1" in selected_gene:
            st.info(
                "**CTNNB1 is a Protein-Protein Interaction (PPI) target.**  \n"
                "PPI inhibitors are structurally larger and frequently exceed standard "
                "Lipinski thresholds (Rule of Five). The screening pipeline applied a "
                "relaxed filter (≤ 2 violations). Compounds passing **beyond Rule-of-Five "
                "(bRo5)** criteria are biologically expected and valid for this target."
            )

        props = compute_lipinski(smiles)
        if props is not None:
            render_lipinski_card(props)
        else:
            st.info("Molecular properties unavailable — install `rdkit-pypi`.")

        st.markdown("<br>", unsafe_allow_html=True)

        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{selected_cid}"
        st.markdown(
            f'<div style="text-align:center;">'
            f'<a href="{pubchem_url}" target="_blank" '
            f'style="background:#1f6feb; color:white; padding:8px 20px; '
            f'border-radius:6px; text-decoration:none; font-size:0.82rem; '
            f'font-weight:600; display:inline-block;">'
            f"🔗 View on PubChem</a></div>",
            unsafe_allow_html=True,
        )

st.divider()

# ── Score distribution by top 100 (zoomed) ────────────────────────────────────
st.markdown(
    '<div class="section-header">Score Distribution — Top 100 Candidates</div>',
    unsafe_allow_html=True,
)

fig2 = px.bar(
    df_top100.reset_index(),
    x="index",
    y="prediction_score",
    color="prediction_score",
    color_continuous_scale=["#1f6feb", "#58a6ff", "#3fb950"],
    labels={"index": "Rank", "prediction_score": "Prediction Score"},
    title=f"{gene_short} — Ranked Score Profile (Top 100)",
    template="plotly_dark",
)
fig2.update_layout(
    plot_bgcolor="#161b22",
    paper_bgcolor="#161b22",
    font={"family": "Inter", "color": "#e6edf3"},
    title_font_size=13,
    coloraxis_showscale=False,
    xaxis=dict(gridcolor="#30363d", linecolor="#30363d", title="Rank"),
    yaxis=dict(gridcolor="#30363d", linecolor="#30363d"),
    margin=dict(l=20, r=20, t=50, b=20),
)
st.plotly_chart(fig2, use_container_width=True)

# ── Molecular Docking Results ──────────────────────────────────────────────────
st.divider()
st.markdown('<div class="section-header">Molecular Docking Results</div>', unsafe_allow_html=True)

_BASE      = Path(__file__).parent
_DOCK_JSON = _BASE / "docking_summary_exp02.json"
_EXP02_DIR = _BASE
_PROT_DIR  = _BASE

# exp02: her iki protein de tamamlandı
_GENE_TO_PROT = {
    "MYC (Transcription Factor)": "MYC",
    "CTNNB1 (Wnt Signaling)"    : "CTNNB1",
}
_GENE_TO_PDB = {
    "MYC (Transcription Factor)": "6G6K",
    "CTNNB1 (Wnt Signaling)"    : "1JPW",
}
_GENE_COLOR = {
    "MYC (Transcription Factor)": "#1f6feb",   # mavi — MYC
    "CTNNB1 (Wnt Signaling)"    : "#388bfd",   # açık mavi — CTNNB1
}
_GENE_IFACE = {
    "MYC (Transcription Factor)": "MYC–MAX interface (6G6K, chains A+C × B+D)",
    "CTNNB1 (Wnt Signaling)"    : "β-catenin–TCF4 interface (1JPW, chain A × chain D)",
}

_prot_key  = _GENE_TO_PROT.get(selected_gene, "MYC")
_pdb_code  = _GENE_TO_PDB.get(selected_gene, "6G6K")
_rec_color = _GENE_COLOR.get(selected_gene, "#1f6feb")

if not _DOCK_JSON.exists():
    st.warning(f"exp02 JSON bulunamadı: `{_DOCK_JSON}`  — `docking_pipeline2.py` çalıştırın.")
else:
    with open(_DOCK_JSON) as _f:
        _dock_raw = json.load(_f)

    _gene_dock = sorted(
        [d for d in _dock_raw if d.get("protein") == _prot_key],
        key=lambda x: x["best_affinity"],
    )

    if not _gene_dock:
        st.info(f"JSON'da {_prot_key} sonucu bulunamadı.")
    else:
        _RT       = 0.001987 * 298   # RT @ 298 K (kcal/mol)
        _ML_SCORE = {155548927: 0.9492, 137355896: 0.9405, 124117970: 0.9164,
                     91515241:  0.8944, 46206285:  0.8941}
        _ML_RANK  = {155548927: 1, 137355896: 2, 124117970: 3, 91515241: 4, 46206285: 5}
        _COLORS   = ["#3fb950", "#58a6ff", "#f0883e", "#da3633", "#8b949e"]

        # ── exp02 badge + intro + metric cards ────────────────────────────────
        _intro_l, _intro_r = st.columns([3, 1], gap="large")

        with _intro_l:
            _iface_txt = _GENE_IFACE.get(selected_gene, "")
            st.markdown(
                f"""
                <div style="background:#161b22; border:1px solid #30363d; border-radius:10px;
                    padding:18px 22px; font-size:0.83rem; color:#c9d1d9; line-height:1.85;">
                    <span style="background:#238636; color:#e6edf3; border-radius:4px;
                        padding:2px 8px; font-size:0.72rem; font-weight:600;
                        margin-bottom:10px; display:inline-block;">exp02</span>
                    &nbsp;
                    <b style="color:#e6edf3;">AutoDock Vina v1.2.5</b> — top-5 ML candidates
                    docked against the
                    <b style="color:#58a6ff;">{_iface_txt}</b>.<br>
                    Grid: <code>30×30×30 Å³</code> &nbsp;·&nbsp;
                    Exhaustiveness: <code>64</code> &nbsp;·&nbsp;
                    Poses per ligand: <code>9</code> &nbsp;·&nbsp;
                    Run in parallel (MYC ∥ CTNNB1).
                </div>
                """,
                unsafe_allow_html=True,
            )

        with _intro_r:
            _best    = _gene_dock[0]
            _best_ki = math.exp(_best["best_affinity"] / _RT) * 1e6
            st.metric("Best Affinity",  f"{_best['best_affinity']:.3f} kcal/mol")
            st.metric("Best Candidate", f"CID {_best['cid']}")
            st.metric("Estimated Kᵢ",   f"~{_best_ki:.1f} µM")

        st.markdown("<br>", unsafe_allow_html=True)

        # ── Results table ──────────────────────────────────────────────────────
        st.markdown(
            '<div style="font-size:0.85rem; font-weight:600; color:#e6edf3; margin-bottom:10px;">'
            'Docking Affinity Ranking</div>',
            unsafe_allow_html=True,
        )

        _rows_html = ""
        for _i, _mol in enumerate(_gene_dock):
            _cid   = _mol["cid"]
            _aff   = _mol["best_affinity"]
            _ki    = math.exp(_aff / _RT) * 1e6
            _mlr   = _ML_RANK.get(_cid, "?")
            _mls   = _ML_SCORE.get(_cid, 0)
            _poses = _mol["all_poses"][1:]
            _close = sum(1 for p in _poses if p["rmsd_ub"] <= 4.0)
            _conf  = "High" if _close >= 4 else ("Medium" if _close >= 2 else "Low")
            _bc    = _COLORS[_i]
            _cc    = "#3fb950" if _conf == "High" else ("#f0883e" if _conf == "Medium" else "#8b949e")
            _star  = " ★" if _i == 0 else ""
            _rows_html += f"""
            <tr style="border-bottom:1px solid #21262d;">
                <td style="padding:10px 14px; font-weight:700; color:{_bc};">{_i+1}{_star}</td>
                <td style="padding:10px 14px; font-family:'JetBrains Mono',monospace;
                    color:#58a6ff;">{_cid}</td>
                <td style="padding:10px 14px; font-weight:600; color:#e6edf3;">{_aff:.3f}</td>
                <td style="padding:10px 14px; color:#8b949e;">~{_ki:.1f} µM</td>
                <td style="padding:10px 14px; color:#8b949e;">ML:{_mlr} ({_mls:.4f})</td>
                <td style="padding:10px 14px; font-weight:600; color:{_cc};">{_conf}</td>
            </tr>"""

        st.markdown(
            f"""
            <div style="overflow-x:auto; border-radius:10px; border:1px solid #30363d;">
            <table style="width:100%; border-collapse:collapse; background:#161b22;">
                <thead>
                <tr style="background:#21262d; border-bottom:2px solid #30363d;">
                    <th style="padding:10px 14px; text-align:left; font-size:0.75rem;
                        color:#8b949e; font-weight:600;">Rank</th>
                    <th style="padding:10px 14px; text-align:left; font-size:0.75rem;
                        color:#8b949e; font-weight:600;">PubChem CID</th>
                    <th style="padding:10px 14px; text-align:left; font-size:0.75rem;
                        color:#8b949e; font-weight:600;">Affinity (kcal/mol)</th>
                    <th style="padding:10px 14px; text-align:left; font-size:0.75rem;
                        color:#8b949e; font-weight:600;">Est. Kᵢ</th>
                    <th style="padding:10px 14px; text-align:left; font-size:0.75rem;
                        color:#8b949e; font-weight:600;">ML Score</th>
                    <th style="padding:10px 14px; text-align:left; font-size:0.75rem;
                        color:#8b949e; font-weight:600;">Pose Confidence</th>
                </tr>
                </thead>
                <tbody>{_rows_html}</tbody>
            </table>
            </div>
            """,
            unsafe_allow_html=True,
        )

        st.markdown("<br>", unsafe_allow_html=True)

        # ── Charts ─────────────────────────────────────────────────────────────
        _ch_l, _ch_r = st.columns(2, gap="large")

        with _ch_l:
            _cids_s     = [str(d["cid"]) for d in _gene_dock]
            _affs       = [d["best_affinity"] for d in _gene_dock]
            _bar_colors = [_COLORS[i] for i in range(len(_gene_dock))]
            _fig_aff = go.Figure(go.Bar(
                x=_affs, y=_cids_s, orientation="h",
                marker_color=_bar_colors,
                text=[f"{a:.3f}" for a in _affs],
                textposition="auto",
                textfont=dict(color="#e6edf3", size=11),
            ))
            _fig_aff.update_layout(
                title=f"Binding Affinity — {_pdb_code} (exp02)",
                title_font_size=12,
                template="plotly_dark",
                plot_bgcolor="#161b22", paper_bgcolor="#161b22",
                font=dict(family="Inter", color="#e6edf3", size=11),
                xaxis=dict(gridcolor="#30363d", title="kcal/mol", autorange="reversed"),
                yaxis=dict(gridcolor="#30363d", title="CID"),
                margin=dict(l=10, r=10, t=40, b=10),
                height=300,
            )
            st.plotly_chart(_fig_aff, use_container_width=True)

        with _ch_r:
            _sc_fig = go.Figure()
            for _i, _mol in enumerate(_gene_dock):
                _sc_fig.add_trace(go.Scatter(
                    x=[p["rmsd_ub"]  for p in _mol["all_poses"]],
                    y=[p["affinity"] for p in _mol["all_poses"]],
                    mode="markers+text",
                    name=f"CID {_mol['cid']}",
                    marker=dict(color=_COLORS[_i], size=9, opacity=0.85),
                    text=[str(p["mode"]) for p in _mol["all_poses"]],
                    textposition="top center",
                    textfont=dict(size=8),
                ))
            _sc_fig.add_vline(x=2.0, line_dash="dash", line_color="#484f58",
                              annotation_text="2 Å threshold", annotation_font_size=9)
            _sc_fig.update_layout(
                title="Pose Landscape (RMSD ub vs Affinity)",
                title_font_size=12,
                template="plotly_dark",
                plot_bgcolor="#161b22", paper_bgcolor="#161b22",
                font=dict(family="Inter", color="#e6edf3", size=11),
                xaxis=dict(gridcolor="#30363d", title="RMSD upper bound (Å)"),
                yaxis=dict(gridcolor="#30363d", title="Affinity (kcal/mol)"),
                legend=dict(font=dict(size=9), bgcolor="#161b22"),
                margin=dict(l=10, r=10, t=40, b=10),
                height=300,
            )
            st.plotly_chart(_sc_fig, use_container_width=True)

        # ── 3D Viewer ──────────────────────────────────────────────────────────
        st.markdown(
            '<div style="font-size:0.85rem; font-weight:600; color:#e6edf3; margin-bottom:10px;">'
            '3D Docking Pose Viewer</div>',
            unsafe_allow_html=True,
        )

        _cid_opts   = [d["cid"] for d in _gene_dock]
        _cid_labels = {d["cid"]: f"CID {d['cid']}  ({d['best_affinity']:.3f} kcal/mol)"
                       for d in _gene_dock}
        _sel_cid = st.selectbox(
            f"Select ligand to visualize — {_pdb_code} mode 1 (best pose):",
            options=_cid_opts,
            format_func=lambda x: _cid_labels[x],
            key=f"sel_cid_{_prot_key}",
        )

        _rec_path = _PROT_DIR / f"{_pdb_code}_receptor.pdbqt"
        _lig_path = _EXP02_DIR / f"{_prot_key}_CID_{_sel_cid}_out.pdbqt"

        if not _rec_path.exists() or not _lig_path.exists():
            _missing = []
            if not _rec_path.exists(): _missing.append(f"`{_rec_path.name}`")
            if not _lig_path.exists(): _missing.append(f"`{_lig_path.name}`")
            st.markdown(
                f"""
                <div style="background:#161b22; border:1px dashed #30363d; border-radius:8px;
                    padding:24px; text-align:center; color:#8b949e;">
                    PDBQT file(s) not found: {', '.join(_missing)}<br>
                    <small>Run <code>docking_pipeline2.py</code> to generate exp02 outputs.</small>
                </div>
                """,
                unsafe_allow_html=True,
            )
        else:
            # ── PDBQT → temiz PDB dönüşümü ────────────────────────────────────
            # 3Dmol.js 'pdbqt' formatını TANIMAZ → gri küre hatası.
            # Çözüm: Python'da PDBQT'ye özgü satırları (ROOT/BRANCH/TORSDOF/REMARK)
            # at, ATOM sütunlarını 66 karakterde kes (kısmi yük + AD4 tipi kaldır).

            def _pdbqt_to_pdb(path: Path, ligand_mode: bool = False) -> str:
                """PDBQT → 3Dmol.js-uyumlu PDB string."""
                _keep   = ("ATOM", "HETATM", "TER", "END")
                _skip   = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH",
                           "TORSDOF", "REMARK", "MODEL", "ENDMDL")
                out     = []
                in_m1   = not ligand_mode   # reseptörde MODEL etiketi yok

                for ln in path.read_text(errors="ignore").splitlines():
                    tag = ln[:6].strip()
                    if ligand_mode:
                        if tag == "MODEL" and not in_m1:
                            in_m1 = True
                            continue
                        if tag == "ENDMDL" and in_m1:
                            break
                    if not in_m1:
                        continue
                    if tag in _skip:
                        continue
                    if tag in ("ATOM", "HETATM"):
                        # sütun 66 sonrasını kes (AD4 atom tipi + kısmi yük)
                        clean = ln[:66].rstrip().ljust(66)
                        if ligand_mode:
                            clean = "HETATM" + clean[6:]  # liganı HETATM yap
                        out.append(clean)
                    elif tag in ("TER", "END"):
                        out.append(ln[:6])

                # JS template literal içinde backtick ve ${ sorun çıkarır → kaçış
                return "\n".join(out).replace("`", "'").replace("${", "$ {")

            _rec_pdb = _pdbqt_to_pdb(_rec_path, ligand_mode=False)
            _lig_pdb = _pdbqt_to_pdb(_lig_path, ligand_mode=True)

            _html_3d = f"""<!DOCTYPE html>
<html><head>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ margin:0; padding:0; background:#0d1117; overflow:hidden; }}
  #viewer {{ width:100%; height:470px; position:relative; }}
</style>
</head>
<body><div id="viewer"></div>
<script>
(function() {{
  var viewer = $3Dmol.createViewer(
    document.getElementById('viewer'),
    {{ backgroundColor: '0x0d1117', antialias: true }}
  );

  // ── Reseptör: PDB formatında yükle, cartoon ────────────────────────────
  viewer.addModel(`{_rec_pdb}`, 'pdb');
  viewer.setStyle(
    {{ model: 0 }},
    {{ cartoon: {{ color: '{_rec_color}', opacity: 0.85, thickness: 0.4 }} }}
  );

  // ── Ligand: HETATM olarak yükle, stick + sphere ────────────────────────
  viewer.addModel(`{_lig_pdb}`, 'pdb');
  viewer.setStyle(
    {{ model: 1 }},
    {{
      stick:   {{ colorscheme: 'Jmol', radius: 0.18 }},
      sphere:  {{ colorscheme: 'Jmol', radius: 0.32 }}
    }}
  );

  // ── Kamerayı ligandı merkez alacak şekilde konumlandır ─────────────────
  viewer.zoomTo({{ model: 1 }});
  viewer.zoom(0.85);
  viewer.render();
}})();
</script>
</body></html>"""

            components.html(_html_3d, height=478, scrolling=False)
            st.markdown(
                f'<div style="font-size:0.7rem; color:#484f58; text-align:center; margin-top:4px;">'
                f'Cartoon: {_iface_txt} &nbsp;·&nbsp; '
                f'Colored sticks + spheres: CID {_sel_cid} — mode 1 &nbsp;·&nbsp; '
                f'Powered by <a href="https://3dmol.org" target="_blank" '
                f'style="color:#484f58;">3Dmol.js</a></div>',
                unsafe_allow_html=True,
            )

        # ── Download buttons ───────────────────────────────────────────────────
        st.markdown("<br>", unsafe_allow_html=True)
        _dl1, _dl2, _dl3 = st.columns(3)

        _best_lig_dl = _EXP02_DIR / f"{_prot_key}_CID_{_gene_dock[0]['cid']}_out.pdbqt"
        with _dl1:
            if _best_lig_dl.exists():
                st.download_button(
                    f"⬇ Best Pose PDBQT ({_prot_key})",
                    data=_best_lig_dl.read_text(),
                    file_name=f"{_prot_key}_CID_{_gene_dock[0]['cid']}_exp02_best.pdbqt",
                    mime="text/plain",
                    use_container_width=True,
                )
        with _dl2:
            st.download_button(
                "⬇ Docking Summary (exp02)",
                data=_DOCK_JSON.read_text(),
                file_name="docking_summary_exp02.json",
                mime="application/json",
                use_container_width=True,
            )
        with _dl3:
            if _rec_path.exists():
                st.download_button(
                    f"⬇ Receptor PDBQT ({_pdb_code})",
                    data=_rec_path.read_text(),
                    file_name=f"{_pdb_code}_receptor.pdbqt",
                    mime="text/plain",
                    use_container_width=True,
                )

# ── Footer ─────────────────────────────────────────────────────────────────────
st.markdown(
    """
    <div class="footer-bar">
        <span class="footer-icon">⏳</span>
        <b>MD Simulations (100ns) results are currently being processed.</b>
        &nbsp; Results will be incorporated into the candidate ranking upon completion.
        <br><br>
        <span style="font-size:0.72rem; color:#484f58;">
            Virtual Screening Pipeline &nbsp;·&nbsp;
            119M PubChem Compounds &nbsp;·&nbsp;
            TÜBİTAK 2209-A Project
        </span>
    </div>
    """,
    unsafe_allow_html=True,
)
