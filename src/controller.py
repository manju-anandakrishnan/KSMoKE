import src.service as service
from src.service import KinaseSubstrateService, KinaseInferenceService
import streamlit as st
import pandas as pd
import os

# Constants of directories
PRED_DIR = 'data/ks_libraries/BL-KSMo'
EVIDENCE_DIR = 'data/ks_evidence'

def get_evidence(substrate,site,kinase):
    '''
    Docstring for get_evidence: This method fetches evidence URL for a given kinase-substrate-site triplet based on data in iPTMnet data
    
    :param substrate: This is the substrate UniProt ID
    :param site: This is the phosphorylation site on the substrate
    :param kinase: This is the kinase UniProt ID
    :return: evidence_URL: This is the evidence URL string if found else empty string
    '''
    evidence_URL = ''
    ks_evidence_df = pd.read_csv(os.path.join(EVIDENCE_DIR,'iPTMnet_ks.tsv'),sep='\t')
    evidence = ks_evidence_df[(ks_evidence_df['substrate_acc']==substrate)
                              & (ks_evidence_df['site']==site)
                              & (ks_evidence_df['kinase_acc']==kinase)]
    if not evidence.empty:
        evidence_URL = substrate
    return evidence_URL


def get_predictions(substrate,motif,site,prob):
    '''
    Docstring for get_predictions: This method fetches predicted kinases for a given substrate-motif pair above the input probability threshold.
    
    :param substrate: This is the substrate UniProt ID of interest
    :param motif: This is the motif sequence surrounding the phosphorylation site
    :param site: This is the phosphorylation site on the substrate
    :param prob: This is the prediction probability threshold
    :return: pred_df: This is the dataframe of predicted kinases with prediction probabilities and evidence URLs (if there is a evidence of kinase-substrate relation in iPTMnet)
    '''
    file_path = os.path.join(PRED_DIR,f'{substrate}.feather')
    pred_df = pd.read_feather(file_path)
    pred_df = pred_df[(pred_df['substrate_motif']==substrate+'_'+motif)
                      & (pred_df['ksf2_pred']>=prob)].copy()
    pred_df = pred_df[['kinase','ksf2_pred']].copy()
    pred_df['ksf2_pred'] = pred_df['ksf2_pred'].round(3)
    pred_df['Kinase Gene'] =  pred_df['kinase'].apply(lambda x:get_protein_gene_dict()[x])
    pred_df['Evidence'] = pred_df.apply(lambda row: get_evidence(substrate,site,row['kinase']), axis=1)
    pred_df.rename({'kinase':'Kinase UniProt ID',
                    'ksf2_pred':'Prediction Probability'}, inplace=True, axis=1)
    pred_df = pred_df.reset_index(drop=True)
    pred_df.drop_duplicates(inplace=True)
    pred_df.sort_values(by='Prediction Probability',ascending=False,inplace=True)
    pred_df = pred_df[['Kinase Gene','Kinase UniProt ID','Prediction Probability','Evidence']]
    return pred_df


def get_substrate_genes(cut_off_prob):
    '''
    Docstring for get_substrate_genes: This method fetches all substrate gene names from the kinase-substrate prediction dataset above the input probability threshold.

    :param cut_off_prob: This is the prediction probability threshold
    :return: substrate_genes: This is the list of substrate gene names
    '''
    substrate_genes = ['']
    substrate_genes.extend(service.get_substrate_genes(cut_off_prob))
    substrate_genes.sort()
    return substrate_genes

@st.cache_data
def get_protein_gene_dict():
    '''
    Docstring for get_protein_gene_dict: This method fetches a dictionary mapping protein UniProt IDs to gene names from the gene_protein.csv file in data folder.
    :return: protein_gene_dict: This is the dictionary mapping protein UniProt IDs to gene names
    '''
    df = pd.read_csv('data/gene_protein.csv',sep='|')
    protein_gene_dict = dict(zip(df['protein'].to_list(),df['gene'].to_list()))
    return protein_gene_dict

class Controller:

    def __init__(self):
        '''
        Docstring for __init__: This is the initialization method for Controller class.
        
        :param self: This is the instance of the class
        '''
        self.ks_service = KinaseSubstrateService()
        self.ki_service = KinaseInferenceService()

    def get_proteins(self,gene,cut_off_prob):
        '''
        Docstring for get_proteins: This method invokes the service to fetch all protein UniProt IDs for a given gene name from the kinase-substrate prediction dataset above the input probability threshold.
        
        :param self: This is the instance of the class
        :param gene: This is the gene name
        :param cut_off_prob: This is the prediction probability threshold
        :return: proteins: This is the list of protein UniProt IDs for the given gene name
        '''
        protein_phosphosites = service.get_protein_phosphosites(cut_off_prob)
        protein_phosphosites = protein_phosphosites[gene]
        proteins = protein_phosphosites.keys()
        return proteins

    def get_phosphosite(self,gene, protein, cut_off_prob):
        '''
        Docstring for get_phosphosite: This method invokes the service to fetch all phosphosites and their motifs for a given gene name and protein UniProt ID from the kinase-substrate prediction dataset above the input probability threshold.
        
        :param self: This is the instance of the class
        :param gene: This is the gene name
        :param protein: This is the protein UniProt ID
        :param cut_off_prob: This is the prediction probability threshold
        :return: site_motif_pairs: This is the dictionary of phosphosites and their motifs for the given gene name and protein UniProt ID
        '''
        phosphosites = service.get_protein_phosphosites(cut_off_prob)
        protein_phosphosites = phosphosites[gene][protein]
        site_motif_pairs = {site:m for site,m in protein_phosphosites}
        return site_motif_pairs

    def get_predicted_kinases(self,substrate_protein,motif,site,prediction_prob_cutoff):
        return get_predictions(substrate_protein,motif,site,prediction_prob_cutoff)
    
    def get_kinase_substrate_links(self, kinase, substrate):
        '''
        Docstring for get_kinase_substrate_links: This method calls the service method to fetch kinase-substrate links for a given kinase and substrate.

        :param self: This is the instance of the class
        :param kinase: This is the kinase UniProt ID
        :param substrate: This is the substrate UniProt ID
        '''
        return self.ks_service.get_kinase_substrate_links(kinase,substrate)
    
    def get_kinase_substrate_network(self, kinase_list, substrate_list, relationship_type=None):
        '''
        Docstring for get_kinase_substrate_network: This method calls the service method to fetch kinase-substrate network for given lists of kinases and list of substrates.

        :param self: This is the instance of the class
        :param kinase_list: This is the list of kinase UniProt IDs
        :param substrate_list: This is the list of substrate UniProt IDs
        :param relationship_type: This is the type of relationship to filter by (e.g., 'participates_in', 'part_of_complex')
        :return: triples: This is the list of kinase-substrate triples representing the network
        '''
        return self.ks_service.get_kinase_substrate_network(kinase_list, substrate_list, relationship_type)
    
    def get_enriched_kinases(self, df, bg_ks_library_key, bg_df=None, logFC=None, pval=None):
        '''
        Docstring for get_enriched_kinases: This method calls the service method to perform kinase enrichment analysis on the input experimental dataframe with optional logFC and p-value thresholds.

        :param self: This is the instance of the class
        :param df: This is the experimental dataframe
        :param logFC: This is the log fold change threshold
        :param pval: This is the p-value threshold
        return: input_exp_df_shape: This is the number of rows in the input experimental dataframe
                result_df: This is the dataframe of enriched kinases with enrichment statistics
                exp_df_shape: This is the number of rows in the filtered experimental dataframe used for enrichment analysis
                in_exp_df: This is the filtered experimental dataframe used for enrichment analysis
        '''
        input_exp_df_shape = df.shape[0]
        enrichment_results, exp_df_shape, in_exp_df = self.ki_service.get_enriched_kinases(df,bg_ks_library_key,bg_df,logFC,pval)
        result_df = pd.DataFrame.from_dict(enrichment_results, orient='index')
        result_df = result_df.reset_index()
        result_df = result_df.rename(columns={'index': 'Kinase UniProt ID'})
        result_df['Kinase Gene'] =  result_df['Kinase UniProt ID'].apply(lambda x:get_protein_gene_dict().get(x))
        result_df.sort_values(by=['adj_p-value'],inplace=True)
        result_df = result_df[['Kinase Gene','Kinase UniProt ID','predicted_sites','p-value','adj_p-value']]
        result_df = result_df[result_df['adj_p-value'] < 0.05].copy().reset_index(drop=True)
        result_df['p-value'] = result_df['p-value'].apply(lambda x: '{:.2e}'.format(x))
        result_df['adj_p-value'] = result_df['adj_p-value'].apply(lambda x: '{:.2e}'.format(x))
        result_df['predicted_sites'] = result_df['predicted_sites'].astype('str')
        result_df.rename(columns={'predicted_sites':'Count of sites targeted by the kinase'}, inplace=True)
        return input_exp_df_shape, result_df, exp_df_shape, in_exp_df
    
