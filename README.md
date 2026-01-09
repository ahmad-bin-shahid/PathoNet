# Overview
PathoNet is a web-based bioinformatics application designed to translate disease phenotypes into integrated molecular and therapeutic insights. The platform implements a systems biology workflow by combining genomic association evidence, functional pathway enrichment, protein–protein interaction networks, and drug–target relationships. It is intended for exploratory research, hypothesis generation, and translational bioinformatics applications.
The application is implemented using Python and Streamlit and relies on curated biomedical databases to ensure interpretability and biological validity.

# Core Design Philosophy
PathoNet follows a disease-centric, multi-layer analytical model, where each analytical layer builds upon the previous one. Rather than analyzing genes in isolation, the platform contextualizes disease biology across molecular functions, interaction networks, and pharmacological space.

# Application Workflow
# 1. Disease Input
--> The workflow begins with a user-defined disease or clinical condition. This disease is mapped to standardized ontology identifiers to ensure compatibility with downstream databases and APIs.

# 2. Associated Genes Module
Objective: Identify genes strongly associated with the selected disease.
Data Source:
--> Open Targets Platform (GraphQL API)
Methodology:
--> Disease–gene associations are retrieved using curated genetic, experimental, clinical, and literature evidence.
--> Each gene is assigned an association score (0–1), reflecting the strength and consistency of evidence across multiple data types.
--> Genes are ranked in descending order of association score to prioritize biologically relevant targets.
Output:
--> Ranked gene association table with evidence-backed relevance scores.

# 3. Functional Enrichment Pathways
Objective: Translate gene-level signals into biological mechanisms.
Data Source:
--> Reactome Knowledgebase (REST API)
Methodology:
--> The top disease-associated genes are subjected to over-representation analysis (ORA).
--> Enrichment is evaluated against curated pathway gene sets using statistical significance testing.
--> Pathways passing significance thresholds are retained and ranked.
Output:
--> Statistically enriched pathways linked to disease biology.
--> Graph-based visualization representing pathway relationships and shared gene overlap.

# 4. Gene Interaction Networks
Objective: Identify molecular interactions and regulatory hubs.
Data Source:
--> STRING Database (REST API)
Methodology:
--> Protein–protein interaction data are retrieved for disease-associated genes.
--> Only interactions above a defined confidence threshold are included.
A network graph is constructed where:
--> Nodes represent genes
--> Edges represent functional or physical interactions
Interpretation:
--> Highly connected nodes (network hubs) may represent key regulators or control points in disease progression.
Output:
--> Interactive gene interaction network visualization.

# 5. Drug Repurposing Analysis
Objective: Identify existing drugs that can modulate disease-associated targets.
Data Source:
ChEMBL via Open Targets integration
Methodology:
--> Disease-associated genes are cross-referenced with known drug targets.
--> Drugs are aggregated across multiple targets to highlight multi-target therapeutic potential.
Prioritization considers:
--> Target relevance
--> Strength of drug–gene interaction
Clinical development status
--> Output:
Ranked list of candidate drugs with associated target information.

# Core Libraries
pandas
numpy
requests
plotly
pyvis
External APIsOpen Targets
Reactome
STRING
ChEMBL

# Intended Use
PathoNet is designed for:
--> Systems biology exploration
--> Translational bioinformatics research
--> Drug repurposing hypothesis generation
--> Educational and academic demonstration
--> The platform is not intended for direct clinical decision-making.
