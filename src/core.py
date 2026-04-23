import streamlit as st
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

class KinaseEnrichmentCore:

    '''
    This method loads the background BL-KSMo kinase-substrate library to the cache
    '''
    @st.cache_data
    def get_ks_library_ksmo():
        ks_library_df = pd.read_feather('data/ks_libraries/BL-KSMo.feather')
        ks_library_df = ks_library_df[['kinase','substrate','site']]
        ks_library_df.drop_duplicates(inplace=True)
        return ks_library_df
    
    '''
    This method loads the background BL-stKSMoyKA kinase-substrate library to the cache
    '''
    @st.cache_data
    def get_ks_library_stksmoyka():
        ks_library_df = pd.read_feather('data/ks_libraries/BL-stKSMoyKA.feather')
        ks_library_df = ks_library_df[['kinase','substrate','site']]
        ks_library_df.drop_duplicates(inplace=True)
        return ks_library_df
    
    '''
    This method loads the background BL-nKIN kinase-substrate library to the cache
    '''
    @st.cache_data
    def get_ks_library_nkin():
        ks_library_df = pd.read_feather('data/ks_libraries/BL-nKIN.feather')
        ks_library_df = ks_library_df[['kinase','substrate','site']]
        ks_library_df.drop_duplicates(inplace=True)
        return ks_library_df
    
    '''
    This method loads the background BL-KA kinase-substrate library to the cache
    '''
    @st.cache_data
    def get_ks_library_ka():
        ks_library_df = pd.read_feather('data/ks_libraries/BL-KA.feather')
        ks_library_df = ks_library_df[['kinase','substrate','site']]
        ks_library_df.drop_duplicates(inplace=True)
        return ks_library_df
    
    '''
    This method loads the background BL kinase-substrate library to the cache
    '''
    @st.cache_data
    def get_ks_library_bl():
        ks_library_df = pd.read_feather('data/ks_libraries/BL.feather')
        ks_library_df = ks_library_df[['kinase','substrate','site']]
        ks_library_df.drop_duplicates(inplace=True)
        return ks_library_df

    def get_ks_library(ks_library_key, custom_bg_df = None):
        ks_library_df = None
        if ks_library_key == "**BL-KSMo**":
            ks_library_df = KinaseEnrichmentCore.get_ks_library_ksmo()
        elif ks_library_key == "**BL-KA**":
            ks_library_df = KinaseEnrichmentCore.get_ks_library_ka()
        elif ks_library_key == "**BL-nKIN**":
            ks_library_df = KinaseEnrichmentCore.get_ks_library_nkin()
        elif ks_library_key == "**BL-KSMoKA**":
            ks_library_df = KinaseEnrichmentCore.get_ks_library_stksmoyka()
        elif ks_library_key == "**Baseline (BL)**":
            ks_library_df = KinaseEnrichmentCore.get_ks_library_bl()
        elif ks_library_key == '**Custom library**':
            ks_library_df = custom_bg_df
        else:
            raise ValueError(f"Invalid background kinase-substrate library key: {ks_library_key}")
        substrate_library = ks_library_df.copy()
        substrate_library.drop(columns=['kinase'],inplace=True,axis=1)
        substrate_library.drop_duplicates(inplace=True)
        library_sites_count = substrate_library.shape[0]
        return ks_library_df, library_sites_count
    

    def __init__(self,exp_df,logFC,pval):
        '''
        This is the initialization method for KinaseEnrichmentCore.
        
        :param self: instance of the class
        :param exp_df: dataframe of experimental phosphorylation sites
        :param logFC: threshold for logFC column
        :param pval: threshold for p-value column
        '''

        # Filter experimental dataframe based on thresholds
        if logFC:
            logFC = float(logFC)
            exp_df = exp_df[exp_df['logFC'].abs() > logFC].copy().reset_index(drop=True)
        if pval:
            pval = float(pval)
            exp_df = exp_df[exp_df['p-value'] < pval].copy().reset_index(drop=True)
        exp_df.drop_duplicates(inplace=True)
        experiment_sites_count = exp_df.shape[0]
        self.experiment_sites_count = experiment_sites_count
        self.exp_df = exp_df


    def _merge_exp_library(self,ks_library_df):
        '''
        This method merges experimental dataframe with kinase-substrate library dataframe 
        and creates two dictionaries with counts of sites targeted by kinases in experimental data and not in the experimental data.
        
        :param self: instance of the class
        :param ks_library_df: kinase-substrate library dataframe
        ''' 
        exp_merged = self.exp_df.merge(ks_library_df,how='outer',left_on=['protein','site'],right_on=['substrate','site'],indicator=True)
        
        in_exp_df = exp_merged[exp_merged._merge.isin(['both'])].copy().reset_index(drop=True)
        not_in_exp_df = exp_merged[exp_merged._merge=='right_only'].copy().reset_index(drop=True)

        self.in_exp_df = in_exp_df
        self.in_exp_kinase_cnt_dict = in_exp_df['kinase'].value_counts().to_dict()
        self.not_in_exp_kinase_cnt_dict = not_in_exp_df['kinase'].value_counts().to_dict()


    def get_contingency_tbl_cnts(self,kinase):
        '''
        This method returns the counts for a contingency table for a given kinase.
        
        :param self: instance of the class
        :param kinase: kinase name
        '''
        a_cnt = self.in_exp_kinase_cnt_dict.get(kinase,0)
        c_cnt = self.experiment_sites_count - a_cnt
        b_cnt = self.not_in_exp_kinase_cnt_dict.get(kinase,0)
        d_cnt = self.library_sites_count - a_cnt - c_cnt - b_cnt
        return a_cnt, b_cnt, c_cnt, d_cnt


    def get_enrichment_results(self, bg_ks_library_key, custom_bg_df = None):

        '''
        This method performs enrichment analysis using Fisher's exact test and returns a dictionary of enrichment results.
        
        :param self: instance of the class
        '''
        kinase_library_df, self.library_sites_count = KinaseEnrichmentCore.get_ks_library(bg_ks_library_key, custom_bg_df)
        self._merge_exp_library(kinase_library_df)

        enrich_results = {}

        kinases = []
        kinases_pval = []

        # For each kinase in the library, compute Fisher's exact test p-value
        for lib_kinase in kinase_library_df['kinase'].unique():

            a_cnt, b_cnt, c_cnt, d_cnt = self.get_contingency_tbl_cnts(lib_kinase)

            if (a_cnt <= 0):
                continue

            if any(x < 0 for x in (a_cnt, b_cnt, c_cnt, d_cnt)) or any(x ==0 for x in (b_cnt+d_cnt, a_cnt+c_cnt)):
                continue

            res = fisher_exact([[a_cnt,b_cnt],[c_cnt,d_cnt]], alternative='greater')
            p_value = res.pvalue
            enrich_results[lib_kinase] = {'p-value':p_value,'adj_p-value':None,
                                        'predicted_sites':a_cnt}
            kinases.append(lib_kinase)
            kinases_pval.append(p_value) 

        # Multiple testing correction using Benjamini-Hochberg
        adjusted_p_values = false_discovery_control(kinases_pval, method='bh')
        for idx,adj_pval in enumerate(adjusted_p_values):
            kinase = kinases[idx]
            enrich_results[kinase]['adj_p-value'] = adj_pval

        return enrich_results, self.exp_df.shape[0], self.in_exp_df


