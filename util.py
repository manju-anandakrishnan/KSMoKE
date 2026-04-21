
import pandas as pd

RELATIONSHIP_TYPE_DICT = {
    "Pathway": "PARTICIPATES_IN_PATHWAY",
    "Complex": "PART_OF_COMPLEX",
    "Biological Process": "BIOLOGICAL_PROCESS",
    "Cellular Component": "CELLULAR_COMPONENT",
    "Molecular Function": "MOLECULAR_FUNCTION",
    "Tissue of Expression": "EXPRESSED_IN"
}

node_category_colors = {'Protein':'#a3b5c7',\
                        'Biological Process':'#4d2f9e',\
                            'Molecular Function':"#3adde6",\
                            'Cellular Component':"#f2d9e3",\
                            'Pathway':'#14ad9b',\
                            'Complex':'#ab3fd5',\
                            'Domain':'#1e6723',\
                            'Homologous Superfamily':'#5638e7',\
                            'Kinase Family':'#94cd85',\
                            'Kinase Family Group':'#af853e',
                            'Kinase Protein':"#ff6600",\
                            'Substrate Protein':"#15ff00",\
                            'Tissue of Expression':"#ddff00"}
                            #'Tissue of Expression':"#ef4c2c53"}

normalized_category_names = {'BioProcess':'Biological Process',\
                            'MolFunc':'Molecular Function',\
                            'CellComp':'Cellular Component',
                            'Tissue':'Tissue of Expression',
                            'HomologousSuperFam':'Homologous Superfamily',
                            'KinaseFam':'Kinase Family',
                            'KinaseFamGroup':'Kinase Family Group',
                            'KinaseProtein':'Kinase Protein',
                            'SubstrateProtein':'Substrate Protein'
                            }


def get_circle(color='black'):
    return f"""<div style='
        width: 20px;
        height: 20px;
        background-color: {color};
        border-radius: 50%;'>
        </div>"""


# def get_table_format(triples, kinases, substrates):
#     '''
#     This method generates a dataframe table from the given triples containing kinase-substrate relationships.
#     :param triples: This is the list of triples containing kinase-substrate relationships
#     :param kinases: This is the set of kinase UniProt IDs
#     :param substrates: This is the set of substrate UniProt IDs
#     :return: table_df: This is the dataframe containing the kinases and substrates grouped by relationship type and related entity'''

#     table_rows = []

#     # Treat each triple as a row in the table and build a table of kinase, substrate, relationship type, and related to columns
#     for triple in triples:
#         head, rel, tail = triple.get()
#         head_category = head.get_category()

#         kinase = substrate = None
#         if head_category == 'Protein' and head.get_id() in kinases:
#             kinase = head.get_name()
#         elif head_category == 'Protein' and head.get_id() in substrates:
#             substrate = head.get_name()

#         data_row = [kinase, substrate, rel.get_label(), tail.get_name()]
#         table_rows.append(data_row)
    
#     table_df = pd.DataFrame(table_rows, columns=['Kinase UniProt ID','Substrate UniProt ID','Relationship Type','Related To'])
#     table_df.dropna(how='all',inplace=True)
#     table_df.reset_index(drop=True,inplace=True)
    
#     # Group by Relationship Type and Related To, aggregating Kinases and Substrates into sets
#     table_df = table_df.groupby(['Relationship Type','Related To']).agg(Kinases=('Kinase UniProt ID',set), Substrates=('Substrate UniProt ID',set)).reset_index()
    
#     # Sort table by number of substrates in descending order
#     table_df.sort_values(by="Substrates",key=lambda x: x.apply(len), ascending=False,inplace=True)
#     table_df = table_df[['Kinases','Substrates','Relationship Type','Related To']]
    
#     # Convert sets of kinases and substrates to comma-separated strings
#     table_df['Kinases'] = table_df['Kinases'].apply(lambda x: ', '.join(kinase for kinase in x if kinase is not None))
#     table_df['Substrates'] = table_df['Substrates'].apply(lambda x: ', '.join(substrate for substrate in x if substrate is not None))
    
#     table_df['Kinase count'] = table_df['Kinases'].apply(lambda x: len(x.split(', ')))
#     table_df['Substrate count'] = table_df['Substrates'].apply(lambda x: len(x.split(', ')))

#     return table_df

def get_table_format(triples, kinases, substrates):
    '''
    This method generates a dataframe table from the given triples containing kinase-substrate relationships.
    :param triples: This is the list of triples containing kinase-substrate relationships
    :param kinases: This is the set of kinase UniProt IDs
    :param substrates: This is the set of substrate UniProt IDs
    :return: table_df: This is the dataframe containing the kinases and substrates grouped by relationship type and related entity'''

    table_rows = []

    # Treat each triple as a row in the table and build a table of kinase, substrate, relationship type, and related to columns
    for triple in triples:
        head, rel, tail = triple.get()
        head_category = head.get_category()

        kinase = substrate = None
        if head_category == 'Protein' and head.get_id() in kinases:
            kinase = head.get_name()
        elif head_category == 'Protein' and head.get_id() in substrates:
            substrate = head.get_name()

        data_row = [kinase, substrate, rel.get_label(), tail.get_name(), tail.get_id()]
        table_rows.append(data_row)
    
    table_df = pd.DataFrame(table_rows, columns=['Kinase UniProt ID','Substrate UniProt ID','Relationship Type','Related To','Related To ID'])
    table_df.dropna(how='all',inplace=True)
    table_df.reset_index(drop=True,inplace=True)
    
    # Group by Relationship Type and Related To, aggregating Kinases and Substrates into sets
    table_df = table_df.groupby(['Relationship Type','Related To', 'Related To ID']).agg(Kinases=('Kinase UniProt ID',set), Substrates=('Substrate UniProt ID',set)).reset_index()
    
    # Sort table by number of substrates in descending order
    table_df.sort_values(by="Substrates",key=lambda x: x.apply(len), ascending=False,inplace=True)
    table_df = table_df[['Kinases','Substrates','Relationship Type','Related To','Related To ID']]
    
    # Convert sets of kinases and substrates to comma-separated strings
    table_df['Kinases'] = table_df['Kinases'].apply(lambda x: ', '.join(kinase for kinase in x if kinase is not None))
    table_df['Substrates'] = table_df['Substrates'].apply(lambda x: ', '.join(substrate for substrate in x if substrate is not None))
    
    table_df['Kinase count'] = table_df['Kinases'].apply(lambda x: len(x.split(', ')))
    table_df['Substrate count'] = table_df['Substrates'].apply(lambda x: len(x.split(', ')))

    return table_df
    