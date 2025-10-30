import streamlit as st
import controller as controller_st
from controller import Controller
import validator
from validator import CustomError
import pandas as pd
import builder
import time
from dotenv import load_dotenv
from st_link_analysis import st_link_analysis
from st_link_analysis.component.layouts import LAYOUTS

load_dotenv()

st.set_page_config(
    page_title="KSMoFinder",
    layout="wide"
)

st.markdown(
    "<h2>KSMoFinder - <span style='color:orange;'>K</span>inase \
        <span style='color:orange;'>S</span>ubstrate \
            <span style='color:orange;'>Mo</span>tif Finder</h2>",
    unsafe_allow_html=True
)

about_page = st.Page('about.py',title='About KSMoFinder')
predict_kinases_page = st.Page('predict_kinases.py',title='Predict Kinases')
kinase_inference_page = st.Page('kinase_inference.py',title='Infer Kinases')

nav_menu = st.navigation([about_page,predict_kinases_page,kinase_inference_page])
nav_menu.run()


# st.markdown(
#     "<h1>KSMoFinder - <span style='color:blue;'>K</span>inase \
#         <span style='color:blue;'>S</span>ubstrate \
#             <span style='color:blue;'>Mo</span>tif Finder</h1>",
#     unsafe_allow_html=True
# )
# st.subheader('Predicting kinases of human phosphosites')
# st.markdown("KSMoFinder learns from biological contexts of proteins and sequence context of motifs surrounding the phosphosites \
#             by embedding a knowledge graph which integrates multiple relationships including\
#             gene ontology (biological processes, cellular location, molecular function), participating pathways and complexes, tissue of expression of proteins, \
#             kinase domain information, and motif specificity of kinases.",unsafe_allow_html=True)


# pred_kinase_tab, kinase_inference_tab = st.tabs(["**Predict kinase**", "**Kinase Inference**"])

# with pred_kinase_tab:

    
#     st.markdown("<b>Predict kinases for a phosphosite and view biological connections between the kinase and the substrate protein</b>",unsafe_allow_html=True)

#     # Session maintenance
#     if 'show_results_pred_ks' not in st.session_state: st.session_state.show_results_pred_ks = False

#     cut_off_col1, cut_off_col2, cut_off_col3 = st.columns([0.9,0.5,1])
#     with cut_off_col1:
#         st.markdown("<div style='padding-top: 35px;'>Kinases with prediction probability <b>>=</b></div>",unsafe_allow_html=True)
#     with cut_off_col2:
#         cut_off_prob_thresh = st.number_input('',value=0.7,min_value=0.0,max_value=1.0)
#     with cut_off_col3:
#         st.markdown("<div style='padding-top: 20px;font-style: italic;'>To view kinases with lower prediction probability, update this value.",unsafe_allow_html=True)
#     kg_substrate_genes = controller_st.get_substrate_genes(0.9)

#     controller = Controller()
#     col1, col2, col3, col4 = st.columns(4)
#     motif = None
#     with col1:
#         substrate_gene = st.selectbox('Substrate (Gene Name):',
#                                         options=kg_substrate_genes)
#     with col2:
#         substrate_protein_options = ['']
#         disable_sp_dd = True
#         if substrate_gene:
#             disable_sp_dd = False
#             substrate_proteins = controller.get_proteins(substrate_gene,cut_off_prob_thresh)
#             substrate_protein_options.extend(substrate_proteins)
#         substrate_protein = st.selectbox('Substrate (UniProt ID):',options=substrate_protein_options,disabled=disable_sp_dd)
#     with col3:
#         site_options = ['']
#         disable_sm_dd = True
#         if substrate_protein:
#             disable_sm_dd = False
#             site_motif_pairs = controller.get_phosphosite(substrate_gene, substrate_protein,cut_off_prob_thresh)
#             site_options.extend(site_motif_pairs.keys())
#         site = st.selectbox('Sites:',options=site_options,disabled=disable_sm_dd)
#     with col4:
#         if site:
#             motif = site_motif_pairs[site]
#             phosphoacceptor = motif[4]
#             st.text('')
#             st.text('')
#             st.markdown(f'''Motif at site: {motif[:4]}:red[**{phosphoacceptor}**]{motif[5:]}''')
        
#     if not motif:
#         st.session_state.show_results_pred_ks = False

#     if st.button('Predict kinases',type='primary'):
#         result_df = controller.get_predicted_kinases(substrate_protein,motif,cut_off_prob_thresh)
#         st.session_state.show_results_pred_ks = True
#         grid_key = str(time.time())
#         st.session_state.results_ks = (result_df, grid_key)


#     def my_call_back() -> None:
#         val = st.session_state["mygraph"]
#         if val["action"] == "expand": 
#             node_ids = val["data"]["node_ids"]
#             # .. handle expand - currently only one node allowed

#     if st.session_state.show_results_pred_ks:
#         data_df, grid_key = st.session_state.results_ks
#         if '::auto_unique_id::' in data_df.columns:
#             data_df.drop(columns=['::auto_unique_id::'],axis=1,inplace=True)
#         result_grid = builder.build_grid(data_df,grid_key)
#         selected_row = result_grid['selected_rows']

            
#         if selected_row is not None:
#             kinase_uniprot_id = selected_row['Kinase UniProt ID'].iloc[0]
#             substrate_uniprot_id = substrate_protein
#             kinase_gene = selected_row['Kinase Gene'].iloc[0]
#             motif = motif
#             with st.container():
#                 st.subheader(f"Biological connections between {kinase_gene} and {substrate_gene}", divider=True)
#                 intermediate_max = 3
#                 try:
#                     ks_triples = controller.get_kinase_substrate_links(kinase_gene,substrate_gene,intermediate_max)
#                     elements = builder.get_graph_elements(ks_triples)

#                     node_styles, edge_styles = builder.get_graph_style(elements)
#                     r_col1, r_col2 = st.columns(2)
#                     with r_col1:
#                         layout = st.selectbox('Layout:',options=sorted(LAYOUTS.keys()),index=4)
                                
#                     st_link_analysis(elements, layout, node_styles, edge_styles,node_actions=['expand'], on_change=my_call_back, key="mygraph")
#                 except Exception as e:
#                     st.error(f"Cannot retrive connections between these proteins:{e}")

# with kinase_inference_tab:

#     st.markdown('''<b>Infer differentially regulated kinases for a condition state using KSMoFinder's predictions</b><br/> 
#                 Enriched kinases are determined using Fisher's exact test.<br/>
#                 For each kinase, <br/>
#                 &nbsp;&nbsp;&nbsp;&nbsp;background includes <u>all sites of the kinase</u> <br/>
#                 &nbsp;&nbsp;&nbsp;&nbsp;foreground includes <u>all sites of the kinase</u> that are in the input file <br/>
#                 <u>'all sites of the kinase'</u> includes KSMoFinder's predicted sites of the kinase and sites reported as catalyzed by the kinase in iPTMnet and PhosphositePlus.   
#                 ''',
#                 unsafe_allow_html=True)

#     if 'show_results_ki' not in st.session_state: st.session_state.show_results_ki = False
#     if 'results_file_id' not in st.session_state: st.session_state.results_file_id = ''
#     if 'input_file_id' not in st.session_state: st.session_state.input_file_id = ''
    
#     st.markdown("""Upload a CSV file with at least two columns - protein and motif.\nThe two columns must be delimited by '|'.<br/>
#                 Optionally, include log2FC and p-value and choose thresholds for the two columns <br/><br/>""",unsafe_allow_html=True)

#     if st.button('Load a sample file'):
#         sample_uploaded_file = 'data/sample_file_kinase_enrich.csv'
#         st.session_state.sample_file = sample_uploaded_file
#         st.session_state.input_file_id = 'sample_file_id'
    
#     user_uploaded_file = st.file_uploader("Choose a CSV file",type=["csv"])

#     if (user_uploaded_file is None) and (st.session_state.input_file_id == 'sample_file_id'):
#         st.markdown("""<b>Data in sample file is displayed below:</b>""",unsafe_allow_html=True)
#         df = pd.read_csv(st.session_state.sample_file,sep='|')
#         st.dataframe(df,hide_index=True)
#         st.session_state.input_psite_df = df

#     if (user_uploaded_file is None) and (st.session_state.input_file_id != 'sample_file_id'):
#         st.session_state.show_results_ki = False

#     if user_uploaded_file is not None:
#         st.markdown("""<b>Data from your file is displayed below:</b>""",unsafe_allow_html=True)
#         df = pd.read_csv(user_uploaded_file,sep='|')
#         st.dataframe(df,hide_index=True)
#         st.session_state.input_psite_df = df
#         st.session_state.sample_file = False
#         st.session_state.input_file_id = user_uploaded_file.file_id

#     if 'input_psite_df' in st.session_state:
#         df = st.session_state.input_psite_df
#         df_columns = df.columns
#         try:
#             log2FC_thresh = None
#             pval_thresh = None
#             validator.validate_input_data(df)
#             if 'log2FC' in df_columns:
#                 log2FC_thresh = st.text_input('Enter absolute threshold of log2FC:')
#             if 'pval' in df_columns:
#                 pval_thresh = st.text_input('Enter pvalue threshold:')
#             if st.button('Infer kinases',type='primary'):
#                 input_fg_sites, result_df, fg_sites = controller.get_dysregulated_kinases(df,log2FC_thresh,pval_thresh)
#                 # Session maintenance
#                 st.session_state.show_results_ki = True
#                 st.session_state.inference_results = result_df
#                 st.session_state.input_fg_sites = input_fg_sites
#                 st.session_state.fg_sites = fg_sites
#                 st.session_state.results_file_id = st.session_state.input_file_id

#             if (st.session_state.show_results_ki) & \
#                     (st.session_state.results_file_id == st.session_state.input_file_id):
#                 result_df = st.session_state.inference_results
#                 input_fg_sites = st.session_state.input_fg_sites
#                 fg_sites = st.session_state.fg_sites

#                 st.text(f'Total number of phosphorylation sites in the input: {input_fg_sites}')
#                 st.text(f'Total number of phosphorylation sites used for inference: {fg_sites}')
#                 builder.build_inference_grid(result_df)
                
#         except CustomError as e:
#             st.error(e)
#             st.session_state.show_results_ki = False
#             st.session_state.inference_results = None
#             st.session_state.input_fg_sites = None
#             st.session_state.fg_sites = None


st.markdown(
    f"""
    <hr style="margin-top: 50px;">
    <div style="text-align: left; color: gray; font-size: 0.9em;">
        GitHub repository: https://github.com/manju-anandakrishnan/KSMoFinder<br/>
        Data repository: https://doi.org/10.5281/zenodo.15730847<br/><br/>
        <b>If you use KSMoFinder's predictions, cite our preprint:</b><br/>
        KSMoFinder - Knowledge graph embedding of proteins and motifs for predicting kinases of human phosphosites<br/>
        Manju Anandakrishnan, Karen E. Ross, Chuming Chen, K. Vijay-Shanker, Cathy H. Wu<br/>
        bioRxiv 2025.10.21.683733; doi: https://doi.org/10.1101/2025.10.21.683733
    </div>
    """,
    unsafe_allow_html=True
)
        


            
