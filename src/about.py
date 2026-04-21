import streamlit as st

st.markdown(
    "<h3><u>About KSMoFinder</h3></u>",
    unsafe_allow_html=True
)

st.markdown("""
            <p>Protein phosphorylation is a common post-translational modification (PTM) that regulates many cellualar process 
            and signaling pathways. This crucial PTM is modulated by intricate regulatory factors mediating interaction of kinase and substrate
            protein. Most research in this space study kinase-substrate relationship with a reductionist approach focusing on sequence similarities 
            between kinases and motifs (short amino acid residues surrounding the phosphorylation site). However, the intricate dependencies between 
            proteins in cellular environment demands looking beyond motif similarities.</p>
            <p>Here, we report KSMoFinder, a supervised neural network classifier model that learns from the biological contexts of proteins and sequence context of motifs 
            to predict kinases of human phosphosites. </p>
            <p>KSMoFinder is trained on the embeddings of a knowledge graph which integrates multiple relationships including gene ontology (
            biological processes, cellular location, molecular function), participating pathways and complexes, tissue of expression of proteins, 
            kinase domain information, and motif specificity of kinases. The semantics of the knowledge graph is learned using different knowledge 
            graph embedding algorithms and the entities are embedded in latent vector space. By bilinear transformation, the vectors of the entities 
            (kinases, substrate(protein)s and motifs) are combined and a supervised neural network classifier is trained to predict phosphorylation 
            probability for a given kinase-phosphosite pair. KSMoFinder offers a broad kinase coverage including 430 human kinases across 9 kinase groups, Atypical, AGC, CAMK, CK1, CMGC, STE, TKL, TK and Other.</p>
            <p>The figure below presents the overview of KSMoFinder's methodology.</p>
            """,unsafe_allow_html=True)
            
st.image("images/KSMoFinder_Overview.png")

st.markdown("""
            <br/>
            <p>The supervised classifier model is trained on 16,137 positive kinase-phosphosite relationships from iPTMnet and PhosphositePlus. Negative training samples of classifier include non-interacting protein pairs and non specific motifs of kinases.</p>
            <p>KSMoFinder's knowledge graph includes 4,881,408 triples including 32 unique relationship types and 360,098 nodes.
            The types of relationships in our knowledge graph and count of triples <head node, relation, tail node> for each relationship type is provided in the table below:</b></p>
            """,unsafe_allow_html=True)

st.image("images/KSMo-KG-data.png")

st.markdown("""
            <p>For additional details about KSMoFinder, refer to our <a href="https://doi.org/10.1093/bioadv/vbaf289">article</a>.</p>
            """,unsafe_allow_html=True)


st.markdown(
    f"""
    <hr style="margin-top: 50px;">
    <div style="text-align: left; color: gray; font-size: 0.9em;">
        GitHub repository: https://github.com/manju-anandakrishnan/KSMoFinder<br/>
        Data repository: https://doi.org/10.5281/zenodo.15730847<br/><br/>
        <b>If you use KSMoFinder's predictions, please cite:</b><br/>
        Manju Anandakrishnan, Karen E Ross, Chuming Chen, K Vijay-Shanker, Cathy H Wu, <br/>
        KSMoFinder—Knowledge graph embedding of proteins and motifs for predicting kinases of human phosphosites, <br/>
        Bioinformatics Advances, 2025;, vbaf289, https://doi.org/10.1093/bioadv/vbaf289
    </div>
    """,
    unsafe_allow_html=True
)