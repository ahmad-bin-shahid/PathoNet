import streamlit as st
import pandas as pd
import requests
import numpy as np
from pyvis.network import Network
import streamlit.components.v1 as components
import plotly.express as px
import base64

# --- 1. RESEARCH-GRADE DARK UI ARCHITECTURE ---
st.set_page_config(
    page_title="PathoNet | Systems Biology Navigator",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Global Consistent Dark Theme CSS
st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&display=swap');
    
    /* Global Background & Text */
    html, body, [class*="css"], .stApp {
        background-color: #0d1117 !important;
        color: #c9d1d9 !important;
        font-family: 'Inter', sans-serif;
    }
    
    /* Typography */
    h1, h2, h3, h4, h5 { font-weight: 500; color: #f0f6fc !important; margin-bottom: 10px; }
    p, li { line-height: 1.8; color: #8b949e; font-size: 15px; }
    strong { color: #f0f6fc; }

    /* Unified Consistent Panels */
    .rationale-panel {
        background-color: #161b22;
        border: 1px solid #30363d;
        padding: 25px;
        margin-bottom: 20px;
    }
    
    .guidance-panel {
        background-color: #0d1117;
        border-left: 4px solid #58a6ff;
        padding: 18px 25px;
        margin-bottom: 30px;
        font-size: 14px;
        color: #c9d1d9;
    }

    .interpretation-panel {
        border-top: 1px solid #30363d;
        margin-top: 40px;
        padding-top: 25px;
    }

    /* Consistent Tab Navigation */
    .stTabs [data-baseweb="tab-list"] { gap: 40px; border-bottom: 1px solid #30363d; }
    .stTabs [data-baseweb="tab"] { 
        background-color: transparent !important; 
        border: none !important; 
        color: #8b949e !important; 
        padding: 12px 0px; 
    }
    .stTabs [aria-selected="true"] { 
        color: #58a6ff !important; 
        border-bottom: 2px solid #58a6ff !important; 
    }

    /* Tables & Frames */
    [data-testid="stDataFrame"] { border: 1px solid #30363d; background-color: #0d1117; }
    
    /* Footer */
    .data-source-footer {
        background-color: #010409;
        border-top: 1px solid #30363d;
        padding: 50px 20px;
        margin-top: 80px;
    }

    /* Inputs */
    .stTextInput input {
        background-color: #0d1117 !important;
        color: white !important;
        border: 1px solid #30363d !important;
    }

    #MainMenu, footer, header {visibility: hidden;}
    </style>
    """, unsafe_allow_html=True)

# --- 2. CONDENSED HEADER ---
try:
    with open("logo.png", "rb") as f:
        logo_b64 = base64.b64encode(f.read()).decode()
    st.markdown(f"""
        <div style="text-align:center; margin-bottom: 0px;">
            <img src="data:image/png;base64,{logo_b64}" width="450">
            <p style="font-weight: 500; color: #58a6ff; margin-top: -15px; letter-spacing: 0.5px; font-size: 14px;">
                A systems-level navigator for disease genomics and therapeutic repurposing.
            </p>
        </div>
    """, unsafe_allow_html=True)
except:
    st.markdown("<h2 style='text-align: center; color: #f0f6fc; margin-bottom:0;'>PATHONET</h2>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: center; color: #58a6ff; margin-top:-15px; font-size:14px;'>A systems-level navigator for disease genomics and therapeutic repurposing.</p>", unsafe_allow_html=True)

st.markdown("<div style='max-width: 900px; margin: 0 auto; text-align: justify; padding: 20px 0 40px 0;'>PathoNet is an integrated computational environment designed for the precision characterization of human disease states. By synthesizing multi-omic evidence from high-throughput repositories, the platform enables researchers to transition from phenotypic queries to actionable molecular landscapes.</div>", unsafe_allow_html=True)

# --- 3. ANALYTICAL LOGIC ---
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

def get_disease_id(query):
    query_body = """query search($q: String!) { search(queryString: $q, entityNames: ["disease"]) { hits { id name } } }"""
    try:
        r = requests.post(OT_URL, json={'query': query_body, 'variables': {"q": query}})
        hits = r.json()['data']['search']['hits']
        return hits[0] if hits else None
    except: return None

def get_molecular_data(efo_id):
    query_body = """
    query molData($efoId: String!) {
      disease(efoId: $efoId) {
        associatedTargets(page: {index: 0, size: 40}) {
          rows { target { approvedSymbol approvedName } score }
        }
        knownDrugs(size: 80) {
          rows {
            drug { name maximumClinicalTrialPhase }
            target { approvedSymbol }
          }
        }
      }
    }"""
    try:
        r = requests.post(OT_URL, json={'query': query_body, 'variables': {"efoId": efo_id}})
        data = r.json()['data']['disease']
        genes = [{"Symbol": row['target']['approvedSymbol'], "Name": row['target']['approvedName'], "Score": round(row['score'], 3)} for row in data['associatedTargets']['rows']]
        raw_drugs = [{"Drug": row['drug']['name'], "Phase": row['drug']['maximumClinicalTrialPhase'], "Target": row['target']['approvedSymbol']} for row in data['knownDrugs']['rows']]
        
        df_d = pd.DataFrame(raw_drugs)
        drugs = df_d.groupby(['Drug', 'Phase'])['Target'].apply(lambda x: ', '.join(sorted(set(x)))).reset_index().sort_values('Phase', ascending=False).to_dict('records') if not df_d.empty else []
        return genes, drugs
    except: return [], []

def get_pathways(symbols):
    try:
        r = requests.post("https://reactome.org/AnalysisService/identifiers/projection", headers={'Content-Type': 'text/plain'}, data="\n".join(symbols))
        token = r.json()['summary']['token']
        res = requests.get(f"https://reactome.org/AnalysisService/token/{token}?pageSize=10").json()
        return [{"Pathway": p['name'][0] if isinstance(p['name'], list) else p['name'], "Significance (-log10 P)": round(-np.log10(p['entities']['pValue']), 2)} for p in res.get('pathways', [])]
    except: return []

def generate_network(symbols):
    try:
        r = requests.get("https://string-db.org/api/json/network", params={"identifiers": "%0d".join(symbols), "species": 9606})
        net = Network(height="600px", width="100%", bgcolor="#0d1117", font_color="#c9d1d9")
        net.barnes_hut(gravity=-5000)
        for edge in r.json():
            net.add_node(edge['preferredName_A'], label=edge['preferredName_A'], color="#58a6ff", size=20)
            net.add_node(edge['preferredName_B'], label=edge['preferredName_B'], color="#58a6ff", size=20)
            net.add_edge(edge['preferredName_A'], edge['preferredName_B'], color="#30363d")
        net.save_graph("network.html")
        return "network.html"
    except: return None

# --- 4. MAIN INTERFACE ---
query = st.text_input("Clinical Phenotype Query Input", placeholder="Enter condition (e.g. Chronic Kidney Disease)")
st.sidebar.markdown("### Technical Parameters")
min_score = st.sidebar.slider("Evidence Stringency", 0.0, 1.0, 0.20)

if query:
    target = get_disease_id(query)
    if target:
        genes, drugs = get_molecular_data(target['id'])
        filtered_genes = [g for g in genes if g['Score'] >= min_score]
        symbols = [g['Symbol'] for g in filtered_genes]

        tab1, tab2, tab3, tab4 = st.tabs(["Genomic Identity", "Pathway Enrichment", "Interaction Network", "Drug Repurposing"])

        with tab1:
            st.markdown('<div class="rationale-panel"><strong>Scientific Rationale: Multi-Evidence Target Identification</strong><br>Prioritization of genomic loci via harmonic weighting of GWAS, somatic mutations, and clinical evidence.</div>', unsafe_allow_html=True)
            st.markdown('<div class="guidance-panel"><strong>Operating Instructions:</strong> Review targets sorted by Association Score. Adjust the sidebar slider to refine stringency.</div>', unsafe_allow_html=True)
            
            if filtered_genes:
                st.dataframe(pd.DataFrame(filtered_genes), use_container_width=True, hide_index=True)
            
            st.markdown('<div class="interpretation-panel"><strong>Clinical Interpretation:</strong> Scores > 0.40 represent high-confidence genetic drivers likely central to the disease etiology.</div>', unsafe_allow_html=True)

        with tab2:
            st.markdown('<div class="rationale-panel"><strong>Scientific Rationale: Over-Representation Analysis</strong><br>Mapping gene sets to the Reactome hierarchy to determine systemic biological failure.</div>', unsafe_allow_html=True)
            st.markdown('<div class="guidance-panel"><strong>Operating Instructions:</strong> Examine the bar chart for statistical significance. -log10 P > 1.30 is significant.</div>', unsafe_allow_html=True)
            
            p_data = get_pathways(symbols)
            if p_data:
                
                df_p = pd.DataFrame(p_data)
                fig = px.bar(df_p, x='Significance (-log10 P)', y='Pathway', orientation='h', color_discrete_sequence=['#58a6ff'])
                fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', font_color="#c9d1d9", yaxis={'autorange':'reversed'})
                st.plotly_chart(fig, use_container_width=True)
            
            st.markdown('<div class="interpretation-panel"><strong>Clinical Interpretation:</strong> Significantly enriched pathways highlight biological mechanisms (e.g., cytokine signaling) compromised in the disease state.</div>', unsafe_allow_html=True)

        with tab3:
            st.markdown('<div class="rationale-panel"><strong>Scientific Rationale: Interaction Topology</strong><br>Visualization of Protein-Protein Interaction (PPI) clusters and regulatory hubs via STRING-DB.</div>', unsafe_allow_html=True)
            st.markdown('<div class="guidance-panel"><strong>Operating Instructions:</strong> Drag nodes to explore connectivity. Densely clustered nodes represent master regulators.</div>', unsafe_allow_html=True)
            
            net_path = generate_network(symbols)
            if net_path:
                
                with open(net_path, 'r', encoding='utf-8') as f:
                    components.html(f.read(), height=600)
            
            st.markdown('<div class="interpretation-panel"><strong>Clinical Interpretation:</strong> Hub proteins (nodes with many connections) are critical therapeutic targets due to their broad regulatory influence.</div>', unsafe_allow_html=True)

        with tab4:
            st.markdown('<div class="rationale-panel"><strong>Scientific Rationale: Pharmacological Cross-Referencing</strong><br>Mapping validated targets to established therapeutic agents in the ChEMBL database.</div>', unsafe_allow_html=True)
            st.markdown('<div class="guidance-panel"><strong>Operating Instructions:</strong> Focus on Phase 4 agents for established safety profiles in drug repurposing contexts.</div>', unsafe_allow_html=True)
            
            if drugs:
                
                df_d = pd.DataFrame(drugs)
                df_d.columns = ["Agent Name", "Max Phase", "Target Loci"]
                st.dataframe(df_d, use_container_width=True, hide_index=True)
            
            st.markdown('<div class="interpretation-panel"><strong>Clinical Interpretation:</strong> Approved molecules targeting relevant genomic drivers represent the most efficient path for off-label clinical trials.</div>', unsafe_allow_html=True)

# --- 5. FOOTER ---
st.markdown("""
<div class="data-source-footer">
    <div style="max-width: 1000px; margin: 0 auto; color: #8b949e; font-size: 12px;">
        <h5 style="color: #f0f6fc;">Data Provenance</h5>
        Open Targets v24.09 | Reactome v89 | STRING-DB v12.0 | ChEMBL v34<br>
        Disclaimer: This tool is for research purposes and does not constitute clinical advice.
    </div>
</div>
""", unsafe_allow_html=True)