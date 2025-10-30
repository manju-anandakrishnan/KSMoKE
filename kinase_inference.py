import streamlit as st
from controller import Controller
import validator
from validator import CustomError
import pandas as pd
import builder


controller = Controller()

st.markdown('''<b>Infer differentially regulated kinases for a condition state using KSMoFinder's predictions</b><br/> 
                Enriched kinases are determined using Fisher's exact test.<br/>
                For each kinase, the below contigency table template is used to compute p-value<br/>
                In the below table, <b>'phosphosites targeted by a kinase'</b> includes KSMoFinder's predicted sites of the kinase and sites reported as catalyzed by the kinase in iPTMnet and PhosphositePlus.   
                ''',
                unsafe_allow_html=True)

st.image("images/FET_contingency_table.png")

if 'show_results_ki' not in st.session_state: st.session_state.show_results_ki = False
if 'results_file_id' not in st.session_state: st.session_state.results_file_id = ''
if 'input_file_id' not in st.session_state: st.session_state.input_file_id = ''
    
st.markdown("""Upload a CSV file with at least two columns - protein and motif.\nThe two columns must be delimited by '|'.<br/>
                Optionally, include log2FC and p-value and choose thresholds for the two columns <br/><br/>""",unsafe_allow_html=True)

if st.button('Load an example data file'):
    sample_uploaded_file = 'data/sample_file_kinase_enrich.csv'
    st.session_state.sample_file = sample_uploaded_file
    st.session_state.input_file_id = 'sample_file_id'
    
user_uploaded_file = st.file_uploader("Choose a CSV file",type=["csv"])

if (user_uploaded_file is None) and (st.session_state.input_file_id == 'sample_file_id'):
    st.markdown("""<b>Data in sample file is displayed below:</b>""",unsafe_allow_html=True)
    df = pd.read_csv(st.session_state.sample_file,sep='|')
    st.dataframe(df,hide_index=True)
    st.session_state.input_psite_df = df

if (user_uploaded_file is None) and (st.session_state.input_file_id != 'sample_file_id'):
    st.session_state.show_results_ki = False

if user_uploaded_file is not None:
    st.markdown("""<b>Data from your file is displayed below:</b>""",unsafe_allow_html=True)
    df = pd.read_csv(user_uploaded_file,sep='|')
    st.dataframe(df,hide_index=True)
    st.session_state.input_psite_df = df
    st.session_state.sample_file = False
    st.session_state.input_file_id = user_uploaded_file.file_id

if 'input_psite_df' in st.session_state:
    df = st.session_state.input_psite_df
    df_columns = df.columns
    try:
        log2FC_thresh = None
        pval_thresh = None
        validator.validate_input_data(df)
        if 'log2FC' in df_columns:
            log2FC_thresh = st.text_input('Enter absolute threshold of log2FC:')
        if 'pval' in df_columns:
            pval_thresh = st.text_input('Enter pvalue threshold:')
        if st.button('Infer kinases',type='primary'):
            input_fg_sites, result_df, fg_sites = controller.get_dysregulated_kinases(df,log2FC_thresh,pval_thresh)
            # Session maintenance
            st.session_state.show_results_ki = True
            st.session_state.inference_results = result_df
            st.session_state.input_fg_sites = input_fg_sites
            st.session_state.fg_sites = fg_sites
            st.session_state.results_file_id = st.session_state.input_file_id
        
        if (st.session_state.show_results_ki) & \
                (st.session_state.results_file_id == st.session_state.input_file_id):
            result_df = st.session_state.inference_results
            input_fg_sites = st.session_state.input_fg_sites
            fg_sites = st.session_state.fg_sites

            st.text(f'Total number of phosphorylation sites in the input: {input_fg_sites}')
            st.text(f'Total number of phosphorylation sites used for inference: {fg_sites}')

            if '::auto_unique_id::' in result_df.columns:
                result_df.drop(columns=['::auto_unique_id::'],axis=1,inplace=True)

            builder.build_inference_grid(result_df)
      
    except CustomError as e:
        st.error(e)
        st.session_state.show_results_ki = False
        st.session_state.inference_results = None
        st.session_state.input_fg_sites = None
        st.session_state.fg_sites = None


        