import streamlit as st
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

st.set_page_config(
    page_title="KSMoFinder",
    layout="wide"
)

about_page = st.Page('src/about.py',title='About KSMoFinder')
predict_kinases_page = st.Page('src/predict_kinases.py',title='Predict Kinases')
kinase_inference_page = st.Page('src/infer_kinases.py',title='Infer Kinases')

# Navigation menu
nav_menu = st.navigation([about_page,predict_kinases_page,kinase_inference_page])
nav_menu.run()
        


            
