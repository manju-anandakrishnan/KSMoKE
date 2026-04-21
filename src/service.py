import src.dao as dao
from src.core import KinaseEnrichmentCore
from src.dao import DatabaseAccess as db_access, KinaseSubstrateDAO as ks_dao

class ServiceError(Exception):
    pass

def get_substrate_genes(cut_off_prob):
    '''
    This method invokes dao and fetches all substrate gene names from the kinase-substrate prediction dataset above the input probability threshold.
    
    :param cut_off_prob: The input probability threshold to filter substrate genes
    :return: substrate_genes: The set of substrate gene names
    '''
    kg_substrates = dao.get_kg_substrates(cut_off_prob)
    substrate_genes = {x[0] for x in kg_substrates}
    return substrate_genes


def get_protein_phosphosites(cut_off_prob):
    '''
    This method invokes dao and fetches all protein UniProt IDs and their phosphosites for a given gene name from the kinase-substrate prediction dataset above the input probability threshold.
    
    :param cut_off_prob: The input probability threshold to filter protein phosphosites
    :return: protein_sites: The dictionary of protein UniProt IDs and their phosphosites
    '''
    kg_substrates = dao.get_kg_substrates(cut_off_prob)
    protein_sites = {}
    for gene,p,site,m in kg_substrates:
        protein_sites.setdefault(gene,{})
        protein_sites[gene].setdefault(p,[])
        protein_sites[gene][p].append((site,m))
    return protein_sites

class KinaseSubstrateService:
    '''
    This class provides services related to kinase-substrate relationships.
    '''

    def __init__(self):
        self.dao = dao.KinaseSubstrateDAO(db_access.get_driver())

    def get_kinase_substrate_links(self,kinase,substrate):
        '''
        This method invokes dao to fetch kinase-substrate links for a given kinase and substrate.
        
        :param self: This is the instance of the class
        :param kinase: The kinase UniProt ID
        :param substrate: The substrate UniProt ID
        '''
        return self.dao.get_ks_links(kinase,substrate)

    def get_kinase_substrate_network(self, kinase_list, substrate_list, relationship_type=None):
        '''
        This method invokes dao to fetch kinase-substrate network for given lists of kinases and list of substrates.

        :param self: This is the instance of the class
        :param kinase_list: The list of kinase UniProt IDs
        :param substrate_list: The list of substrate UniProt IDs
        :param relationship_type: The relationship type to filter by (optional)
        '''
        return self.dao.get_ks_network(kinase_list, substrate_list, relationship_type)

class KinaseInferenceService:
    '''
    This class provides services related to kinase enrichment inference.
    '''

    def __init__(self):
        pass

    def get_enriched_kinases(self, df, bg_ks_library_key, custom_bg_df=None, logFC=None, pval=None):
        '''
        This method performs kinase enrichment analysis using the input dataframe.

        :param self: This is the instance of the class
        :param df: The input dataframe
        :param logFC: The log fold change threshold (optional)
        :param pval: The p-value threshold (optional)
        '''
        ki_core = KinaseEnrichmentCore(df,logFC,pval)
        return ki_core.get_enrichment_results(bg_ks_library_key,custom_bg_df)



