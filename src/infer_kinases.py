import streamlit as st
from src.controller import Controller
import src.validator as validator
from src.validator import CustomError
import pandas as pd
import src.builder as builder
from datetime import datetime
from st_link_analysis import st_link_analysis
from st_link_analysis.component.layouts import LAYOUTS
import src.util as util

st.markdown(
    "<h4><b>Infer differentially regulated kinases for a condition state using one of the compiled libraries or upload a custom library.</b></h4>",
    unsafe_allow_html=True
)

st.markdown('''Enriched kinases are determined using Fisher's exact test.<br/>
                For each kinase, the below contingency table template is used to compute p-value<br/>
                In the below table, <b>'phosphosites targeted by a kinase'</b> includes sites associated to the kinase in the selected/uploaded library.   
                ''',
                unsafe_allow_html=True)

st.image("images/FET_2X2_template.png")

# Session maintenance
if 'show_results_ki' not in st.session_state: st.session_state.show_results_ki = False
if 'visualize_subnetwork' not in st.session_state: st.session_state.visualize_subnetwork = False
if 'results_file_id' not in st.session_state: st.session_state.results_file_id = ''
if 'input_file_id' not in st.session_state: st.session_state.input_file_id = ''
if 'inference_results' not in st.session_state: st.session_state.inference_results = None
if 'network_input_hash' not in st.session_state: st.session_state.network_input_hash = None
if 'ks_triples' not in st.session_state: st.session_state.ks_triples = None

# Experiment phosphoproteomic data input section
with st.container(border=True):
    st.markdown(f"""<b><p style="color:#0000FF;">Upload a phosphoproteome dataset:</p></b>""",unsafe_allow_html=True)

    # File upload section    
    st.markdown("""Upload a CSV file with at least two columns - protein and site.<br/>
                    You file must include a header line with labels - protein, site. <br/>
                    Protein must be UniProt Accession Number. It must be a reviewed entry for organsim, 'Human'. <br/>
                    The column, 'site' must be of the format \<residue\>\<position\>. For example: Y1135, T202, Y1000 <br/>
                    Optionally, include logFC and p-value and choose thresholds for the two columns. <br/>
                    All other columns in the phosphoproteome file will be ignored. <br/>""",unsafe_allow_html=True)

    # Sample file upload
    if st.button('Load a sample data file'):
        sample_uploaded_file = 'data/sample_file_kinase_enrich.csv'
        st.session_state.sample_file = sample_uploaded_file
        st.session_state.input_file_id = 'sample_file_id'

    st.markdown(f"""<b>Upload your CSV file:</b>""",unsafe_allow_html=True) 
    user_uploaded_file = st.file_uploader("experiment_phosphoproteome_file",type=["csv"],label_visibility='collapsed')

    # Display sample data from the uploaded file or sample file
    if (user_uploaded_file is None) and (st.session_state.input_file_id == 'sample_file_id'):
        st.markdown("""<b>Data from the phosphoproteome file:</b>""",unsafe_allow_html=True)
        df = pd.read_csv(st.session_state.sample_file,sep=',')
        st.dataframe(df,hide_index=True)
        st.session_state.input_psite_df = df

    # If there is no uploaded file, don't show previous results
    if (user_uploaded_file is None) and (st.session_state.input_file_id != 'sample_file_id'):
        st.session_state.show_results_ki = False

    # If user uploads a file, read and display the data
    if user_uploaded_file is not None:
        st.markdown("""<b>Data from the phosphoproteome file:</b>""",unsafe_allow_html=True)
        df = pd.read_csv(user_uploaded_file,sep=',')
        st.dataframe(df,hide_index=True)
        st.session_state.input_psite_df = df
        st.session_state.sample_file = False
        st.session_state.input_file_id = user_uploaded_file.file_id

# Background kinase-substrate library file upload
with st.container(border=True):
    st.markdown(f"""<b><p style="color:#0000FF;">Choose (or) upload a background kinase-substrate library:</p></b>""",unsafe_allow_html=True)
    st.markdown(f"""Baseline<b>(BL)</b> refers to curated kinase-substrate relationships from iPTMnet and PhosphositePlus""",unsafe_allow_html=True)
    bg_ks_library_key = st.radio("Kinase-substrate library", 
                            ["**Baseline (BL)**", 
                            "**BL-KSMo**", 
                            "**BL-KA**", 
                            "**BL-KSMoKA**",
                            "**BL-nKIN**",
                            "**Custom library**"], 
                            captions=[
                                "Curated \n\n kinase-substrate \n\n library",
                                "BL + \n\n KSMoFinder's \n\n predictions",
                                "BL + \n\n Kinome Atlas \n\n &nbsp;",
                                "BL + \n\n KSMoFinder's predictions (ST) + \n\n Kinome Atlas (Y)",
                                "BL + \n\n NetworKIN predictions \n\n &nbsp;",
                                "Upload a custom \n\n kinase-substrate library \n\n &nbsp;"], index=1, horizontal=True, label_visibility='collapsed')

    custom_bg_df = None
    if bg_ks_library_key == "**Custom library**":
        st.markdown("""<b>Custom library format specifications:</b><br/>
                    Upload a background file (CSV format with the mandatory columns: kinase, substrate, site). <br/>
                    You file must include a header line with labels - kinase, substrate, site. <br/>
                    The values for kinase and substrate must be UniProt Accession Numbers of the human proteins. The UniProt Accession Numbers must be reviewed entries in UniProt.<br/>
                    The column, 'site' must be of the format \<residue\>\<position\>. For example: Y1135, T202, Y1000 <br/>
                    All other columns in the background file will be ignored.<br/>""",unsafe_allow_html=True)
                                        
        background_uploaded_file = st.file_uploader("background_file",type=["csv"],label_visibility='collapsed',key='background_file')

        if background_uploaded_file is not None:
            st.write("*Background file uploaded successfully. It will be used for kinase enrichment analysis.*")
            st.write("*Background file preview:*")
            custom_bg_df = pd.read_csv(background_uploaded_file)
            st.dataframe(custom_bg_df.head(5),hide_index=True)


# Instantiate the controller
controller = Controller()

# If there is input data, validate and infer kinases
if 'input_psite_df' in st.session_state:
    df = st.session_state.input_psite_df
    df_columns = df.columns
    try:
        logFC_thresh = None
        pval_thresh = None
        validator.validate_input_data(df)
        if 'logFC' in df_columns:
            logFC_thresh = st.text_input('Enter absolute threshold of logFC:',value = '1')
        if 'p-value' in df_columns:
            pval_thresh = st.text_input('Enter p-value threshold:', value='0.05')
        
        if (bg_ks_library_key != '**Custom library**') or ((bg_ks_library_key == '**Custom library**') & (custom_bg_df is not None)):

            # On click of Infer kinases button, perform kinase enrichment analysis and store results in session state
            if st.button('Infer kinases',type='primary'):
                if 'site' in df_columns:
                    df['site'] = df['site'].apply(lambda x:x.upper())
                try:
                    input_exp_sites, result_df, exp_sites, in_exp_df = controller.get_enriched_kinases(df,bg_ks_library_key,bg_df=custom_bg_df,logFC=logFC_thresh,pval=pval_thresh)
                    # Session maintenance
                    st.session_state.show_results_ki = True
                    st.session_state.inference_results = result_df
                    st.session_state.input_exp_sites = input_exp_sites
                    st.session_state.exp_sites = exp_sites
                    st.session_state.results_file_id = st.session_state.input_file_id
                    st.session_state.in_exp_df = in_exp_df

                    # On click of Infer kinases button
                    st.session_state.result_grid_selected_indices = []
                    st.session_state.visualize_subnetwork = False

                    # Unique key for inference results grid. Regenerated on each inference
                    st.session_state.grid_key_ki = f'grid_ki_{datetime.now().strftime("%Y%m%d%H%M%S")}'
                except Exception as e:
                    # In case of exception
                    st.error('An exception occurred when conducting enrichment analysis for your data. \n\n Make sure all records in your phosphoproteomics data aligns with the given instructions. \n\n If you are using a custom library, verify if all records in the library follows the given instructions. \n\n If the issue persits after you have validated, please contact manjua@udel.edu')
        
        # If results are to be shown and the input file has not changed, display the results
        if (st.session_state.show_results_ki) & \
                (st.session_state.results_file_id == st.session_state.input_file_id):
            result_df = st.session_state.inference_results
            input_exp_sites = st.session_state.input_exp_sites
            exp_sites = st.session_state.exp_sites
            in_exp_df = st.session_state.in_exp_df
            
            if '::auto_unique_id::' in result_df.columns:
                result_df.drop(columns=['::auto_unique_id::'],axis=1,inplace=True)
            
            # Display the enrichment table and subsequent results only if the result_df has at least one enriched kinase
            if result_df.shape[0] >= 1:

                st.markdown(f"""<b>Enrichment Results:</b><br/>
                        Total number of phosphorylation sites in the input: {input_exp_sites}<br/>
                        Total number of phosphorylation sites used for inference: {exp_sites}
                         """,unsafe_allow_html=True)

                # Build and display the results grid
                result_grid_df = builder.build_inference_grid(result_df, key=st.session_state.grid_key_ki)
                
                # Download enrichment results
                st.download_button(
                    label="Download enrichment results",
                    data=result_df.iloc[:,:-1].to_csv(index=False).encode('utf-8'),
                    file_name=f'enrichment_results_{datetime.now().strftime("%Y%m%d%H%M%S")}.csv',
                    mime='text/csv',
                    icon=":material/download:"
                )

                st.markdown(f'<b>To visualize the subnetwork of enriched kinases and substrates, select one or more kinases from the table using the checkboxes on the left.</b>',unsafe_allow_html=True)

                # If any rows are selected, store the selected indices and show options to visualize the subnetwork
                if result_grid_df['selected_rows'] is not None:
                    st.session_state.result_grid_selected_indices = result_grid_df['selected_rows'].index.tolist()
                    selected_kinase_uniprot_ids = result_grid_df['selected_rows']['Kinase UniProt ID'].values.tolist()

                    if st.button('Visualize the subnetwork of enriched kinases and substrates',type='primary'):
                        st.session_state.visualize_subnetwork = True

                    if st.session_state.visualize_subnetwork:
                        network_relationship_types = st.multiselect(
                            "Select the relationship type(s) to display in the network:",
                            ("Pathway", "Complex", "Biological Process", "Cellular Component", "Molecular Function", "Tissue of Expression"),
                            #default=["Pathway", "Complex", "Biological Process", "Cellular Component", "Molecular Function", "Tissue of Expression"]
                            default=["Pathway"]
                        )

                        if len(network_relationship_types) >= 1:
                            
                            # Get the kinase-substrate links for the selected kinases
                            ks_network_df = in_exp_df[in_exp_df['kinase'].isin(selected_kinase_uniprot_ids)].copy().reset_index(drop=True)
                            network_kinases = ks_network_df['kinase'].unique().tolist()
                            network_kinases.sort()

                            # Get the substrates connected to the selected kinases
                            network_substrates = ks_network_df['substrate'].unique().tolist()
                            network_substrates.sort()

                            network_relationship_types.sort()

                            with st.container():
                                st.subheader(f"The network of the enriched kinase(s) and the substrate proteins", divider=True)
                                try:
                                    # hash input for session maintenance
                                    network_input_hash = hash(tuple(network_kinases + network_substrates + network_relationship_types))
                                    
                                    # If the input has changed, retrieve the network data from the controller
                                    if ('network_input_hash' not in st.session_state) | (st.session_state.network_input_hash != network_input_hash):                                
                                        ks_triples = []

                                        # Initialize progress bar
                                        progress_text = "Retrieving connections. Please wait..."
                                        percent_complete = 0
                                        ks_network_progress_bar = st.progress(percent_complete, text=progress_text)
                                        
                                        # For the selected relationship types, get connections of the selected kinases to all substrates
                                        total_iterations = len(network_relationship_types) * len(network_substrates)
                                        iteration_count = 0
                                        for network_relationship_type in network_relationship_types:                            
                                            for substrate in network_substrates:
                                                ks_triples.extend(controller.get_kinase_substrate_network(network_kinases,substrate,util.RELATIONSHIP_TYPE_DICT[network_relationship_type]))
                                                # Update progress bar
                                                iteration_count +=1
                                                percent_complete = float(iteration_count / total_iterations)
                                                ks_network_progress_bar.progress(percent_complete, text=progress_text)
                                        
                                        # Complete the progress bar
                                        ks_network_progress_bar.empty()

                                        # Store the input hash in session state
                                        st.session_state.network_input_hash = network_input_hash
                                        st.session_state.ks_triples = ks_triples
                                    else:
                                        ks_triples = st.session_state.ks_triples

                                    elements = builder.get_graph_elements(ks_triples)

                                    # Pass edge caption as label to display labels on edges
                                    if len(network_relationship_types) == 1:
                                        node_styles, edge_styles = builder.get_graph_style(elements, show_edge_label=False, highlight_kinase_ids=network_kinases, highlight_substrate_ids=network_substrates)
                                    else:
                                        node_styles, edge_styles = builder.get_graph_style(elements, show_edge_label=True, highlight_kinase_ids=network_kinases, highlight_substrate_ids=network_substrates)
                                    
                                    # Display legend for the network
                                    st.markdown(f"""
                                        {builder.get_ks_network_legend(network_relationship_types)}
                                        """, 
                                        unsafe_allow_html=True)

                                    r_col1, r_col2 = st.columns(2)
                                    with r_col1:
                                        layout = st.selectbox('Layout:',options=sorted(LAYOUTS.keys()),index=4)
                                    st_link_analysis(elements, layout, node_styles, edge_styles)
                                    
                                    # Display the network
                                    #st_link_analysis(elements, 'cose', node_styles, edge_styles)

                                    # Display table view of the network data
                                    st.markdown(f"""<b>Table view of the network data</b>
                                        """, 
                                        unsafe_allow_html=True)
                                    
                                    # From the triples, get a tabular format of the kinase-substrate network
                                    ks_network_df = util.get_table_format(ks_triples, network_kinases, network_substrates)
                                    
                                    # Display the table grid
                                    builder.build_ksnetwork_grid(ks_network_df)

                                    # Download ks-network results
                                    st.download_button(
                                        label="Download network results",
                                        data=ks_network_df.iloc[:,:-1].to_csv(index=False).encode('utf-8'),
                                        file_name=f'network_results_{datetime.now().strftime("%Y%m%d%H%M%S")}.csv',
                                        mime='text/csv',
                                        icon=":material/download:"
                                    )

                                except Exception as e:
                                    st.error(f"Cannot retrive connections between these proteins:{e}")
            else:
                st.markdown(f"""<p style="color:#FF4D00;">There are no enriched kinase based on the selected kinase-substrate library. <br/>
                            Choose an alternate background and/or validate your phosphoproteome file, and thresholds if applicable.</p>""",unsafe_allow_html=True)
    except CustomError as e:
        st.error(e)
        st.session_state.show_results_ki = False
        st.session_state.inference_results = None
        st.session_state.input_exp_sites = None
        st.session_state.exp_sites = None
        st.session_state.visualize_subnetwork = False
        st.session_state.results_file_id = ''
        st.session_state.input_file_id = ''


        