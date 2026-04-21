import streamlit as st
import controller as controller_st
from controller import Controller
import builder
import time
from dotenv import load_dotenv
from st_link_analysis import st_link_analysis
from st_link_analysis.component.layouts import LAYOUTS

# Load environment variables
load_dotenv()

st.markdown("<b>Predict kinases for a phosphosite and view biological connections between the kinase and the substrate protein</b>",unsafe_allow_html=True)

cut_off_col1, cut_off_col2, cut_off_col3 = st.columns([0.5,0.3,1])
with cut_off_col1:
    st.markdown("<div style='padding-top:35px;text-align:right;'>Kinases with prediction probability <b>>=</b></div>",unsafe_allow_html=True)
with cut_off_col2:
    cut_off_prob_thresh = st.number_input('',value=0.8,min_value=0.0,max_value=1.0)
with cut_off_col3:
    st.markdown("<div style='padding-top: 35px;font-style: italic;'>To view kinases with lower prediction probability, update this value.",unsafe_allow_html=True)

# Get substrate genes for the dropdown
kg_substrate_genes = controller_st.get_substrate_genes(cut_off_prob_thresh)

# Session maintenance
if 'show_results_pred_ks' not in st.session_state: st.session_state.show_results_pred_ks = False
if 'selected_substrate_gene' not in st.session_state:
    st.session_state['selected_substrate_gene'] = kg_substrate_genes[0]
if 'selected_substrate_protein' not in st.session_state:
    st.session_state['selected_substrate_protein'] = ''
if 'selected_site' not in st.session_state:
    st.session_state['selected_site'] = ''


controller = Controller()
col1, col2, col3, col4 = st.columns(4)
motif = None

if st.button('Load an example phosphosite data for prediction'):
    st.session_state.selected_substrate_gene = 'MED14'
    st.session_state.selected_substrate_protein = 'O60244'
    st.session_state.selected_site = 'S1112'

with col1:
    # Substrate gene dropdown
    substrate_gene = st.selectbox('Substrate (Gene Name):',
                                            options=kg_substrate_genes,
                                            key="selected_substrate_gene"
                                            )
with col2:
    substrate_protein_options = ['']
    disable_sp_dd = True
    # Populate substrate protein dropdown based on selected substrate gene
    if substrate_gene:
        disable_sp_dd = False
        substrate_proteins = controller.get_proteins(substrate_gene,cut_off_prob_thresh)
        substrate_protein_options.extend(substrate_proteins)
    substrate_protein = st.selectbox('Substrate (UniProt ID):',
                                     options=substrate_protein_options,
                                     disabled=disable_sp_dd,
                                     key="selected_substrate_protein"
                                     )
with col3:
    site_options = ['']
    disable_sm_dd = True
    # Populate site dropdown based on selected substrate protein
    if substrate_protein:
        disable_sm_dd = False
        site_motif_pairs = controller.get_phosphosite(substrate_gene, substrate_protein,cut_off_prob_thresh)
        site_options.extend(site_motif_pairs.keys())
        sorted_sites = sorted(site_options[1:], key=lambda x: (int(x[1:])))
        site_options = [''] + sorted_sites
    site = st.selectbox('Sites:',
                        options=site_options,
                        disabled=disable_sm_dd,
                        key="selected_site",
                        )
with col4:
    # Display motif of the selected site
    if site:
        motif = site_motif_pairs[site]
        phosphoacceptor = motif[4]
        st.text('')
        st.text('')
        st.markdown(f'''Motif at site: {motif[:4]}:red[**{phosphoacceptor}**]{motif[5:]}''')
        
if not motif:
    st.session_state.show_results_pred_ks = False

# On click of Predict kinases button, get predicted kinases and set session state to show results
if st.button('Predict kinases',type='primary'):
    result_df = controller.get_predicted_kinases(substrate_protein,motif,site,cut_off_prob_thresh)
    st.session_state.show_results_pred_ks = True
    grid_key = str(time.time())
    st.session_state.results_ks = (result_df, grid_key)

# If session state is set to show results, display the predicted kinases grid
if st.session_state.show_results_pred_ks:
    data_df, grid_key = st.session_state.results_ks
    if '::auto_unique_id::' in data_df.columns:
        data_df.drop(columns=['::auto_unique_id::'],axis=1,inplace=True)
    result_grid = builder.build_grid(data_df,grid_key)
    selected_row = result_grid['selected_rows']
    
    # If any row is selected, get the biological connections between the kinase and the substrate protein
    if selected_row is not None:
        kinase_uniprot_id = selected_row['Kinase UniProt ID'].iloc[0]
        substrate_uniprot_id = substrate_protein
        kinase_gene = selected_row['Kinase Gene'].iloc[0]
        motif = motif
        with st.container():
            st.subheader(f"Biological connections between {kinase_gene} and {substrate_gene}", divider=True)
            try:
                # Get kinase-substrate links for the selected kinase and the input substrate, and display the network
                ks_triples = controller.get_kinase_substrate_links(kinase_gene,substrate_gene)
                elements = builder.get_graph_elements(ks_triples)

                node_styles, edge_styles = builder.get_graph_style(elements, 'label')

                # Display the network results
                r_col1, r_col2 = st.columns(2)
                with r_col1:
                    layout = st.selectbox('Layout:',options=sorted(LAYOUTS.keys()),index=4)
               
                st_link_analysis(elements, layout, node_styles, edge_styles)

            except Exception as e:
                st.error(f"Cannot retrive connections between these proteins:{e}")
                st.session_state.show_results_pred_ks = False
                st.session_state.selected_substrate_gene = ''
                st.session_state.selected_substrate_protein = ''
                st.session_state.selected_site = ''