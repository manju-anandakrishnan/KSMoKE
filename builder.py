from st_aggrid import GridOptionsBuilder, AgGrid, JsCode
from st_link_analysis import NodeStyle, EdgeStyle
import util


def build_grid(df,key):
    '''
    Docstring for build_grid: This method builds the grid for predicted kinases of a substrate with a clickable link in the Evidence column.

    :param df: The dataframe containing predicted kinases and probabilities
    :param key: The unique key for the grid
    :return: AgGrid object with the built grid
    '''
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_column(
        "Evidence", "Link to Evidence",
        cellRenderer=JsCode("""
            class UrlCellRenderer {
            init(params) {            
                this.eGui = document.createElement('div');
                if (params.value) {
                    const evidence_URL = `https://research.bioinformatics.udel.edu/iptmnet/entry/${params.value}/#asSub`;
                    this.eGui.innerHTML = `<a href='${evidence_URL}' target='_blank'>Review</a>`;
                } else {
                    this.eGui.innerHTML = '';
                }
            }
            getGui() {
                return this.eGui;
            }
            }
        """)
        )
    gb.configure_selection(selection_mode="single", use_checkbox=False)
    gb.configure_pagination(enabled=True,paginationAutoPageSize=True)

    return AgGrid(df,                  
                #height=min(27+len(df)*35, 400),                
                gridOptions=gb.build(),
                enable_enterprise_modules=False,              
                fit_columns_on_grid_load=True,
                key=str(key),
                custom_css={
                    "#gridToolBar": {
                        "padding-bottom": "0px !important",
                    }
                },
                allow_unsafe_jscode=True
                )


def get_graph_elements(triples):
    '''
    Docstring for get_graph_elements: This method builds the graph elements (nodes and edges) from the given triples.
    
    :param triples: This is the list of triples representing kinase-substrate relationships
    :return: elements: This is the dictionary containing nodes and edges for the graph
    '''
    elements = {'nodes':[], 'edges':[]}
    for triple in triples:
        head, rel, tail = triple.get()
        head_node = {'data':{'id':head.get_id(),'name':head.get_name(),'label':head.get_category()}}
        tail_node = {'data':{'id':tail.get_id(),'name':tail.get_name(),'label':tail.get_category()}}
        elements['nodes'].append(head_node)
        elements['nodes'].append(tail_node)

        ht_edge = {'data':{'id':rel.id,'label':rel.get_label(),'source':head.get_id(),'target':tail.get_id()}}
        elements['edges'].append(ht_edge)
    return elements

def get_edge_labels(edges):
    '''
    Docstring for get_edge_labels: This method extracts unique edge labels from the list of edges.

    :param edges: This is the list of edges in the graph
    :return: edge_labels: This is the set of unique edge labels
    '''
    edge_labels = set()
    for edge in edges:
        edge_labels.add(edge['data']['label'])
    return edge_labels

def get_node_labels(nodes):
    '''
    This method extracts unique node labels from the list of nodes.
    
    :param nodes: This is the list of nodes in the graph
    :return: node_categories: This is the set of unique node labels
    '''
    node_categories = set()
    for node in nodes:
        node_categories.add(node['data']['label'])
    return node_categories

def get_node_id(nodes):
    '''
    This method extracts unique node ids from the list of nodes.
    
    :param nodes: This is the list of nodes in the graph
    :return: node_ids: This is the set of unique node ids   
    '''
    node_ids = set()
    for node in nodes:
        node_ids.add(node['data']['id'])
    return node_ids


def get_graph_style(elements, show_edge_label=True, highlight_kinase_ids=None, highlight_substrate_ids=None):
    '''
    This method configures the node and edge styles for the graph based on the elements provided.
    '''
    
    # If highlight lists are provided, update node labels accordingly for legend display
    if highlight_substrate_ids:
        for node in elements['nodes']:
            if node['data']['id'] in highlight_substrate_ids:
                node['data']['label'] = 'SubstrateProtein'
    if highlight_kinase_ids:
        for node in elements['nodes']:
            if node['data']['id'] in highlight_kinase_ids:
                node['data']['label'] = 'KinaseProtein'
    
    # Update node labels to normalized category names for readability
    for node in elements['nodes']:
        label = node['data']['label']
        node['data']['label'] = util.normalized_category_names.get(label,label)

    # Get unique edge and node labels for styling
    edge_labels = get_edge_labels(elements['edges'])
    node_labels = get_node_labels(elements['nodes'])
    
    # Configure node and edge styles based on labels
    node_styles = [NodeStyle(node_label, util.node_category_colors[node_label],'name') for node_label in node_labels]
    if show_edge_label:
        edge_styles = [EdgeStyle(edge_label, caption='label', color="#e7edf0", directed=True) for edge_label in edge_labels]
    else:
        edge_styles = [EdgeStyle('', color="#f0e7ed", directed=True) for edge_label in edge_labels]
   
    return node_styles, edge_styles

def build_inference_grid(data, key=None):
    '''
    This method builds the grid of kinase enrichment inference results.
    
    :param data: This is the dataframe containing kinase enrichment results
    :param key: This is the unique key for the grid
    :return: AgGrid object with the built grid
    '''
    data.columns = ['Kinase gene',\
                    'Kinase UniProt ID',\
                        'Count of sites targeted by the kinase',\
                              'p-value', \
                                'Adjusted p-value']
    
    gb = GridOptionsBuilder.from_dataframe(data)

    gb.configure_selection(
        selection_mode='multiple',
        use_checkbox=True
    )

    gb.configure_pagination(
        enabled=True,
        paginationAutoPageSize=True
        )

    return AgGrid(
            data, 
            gridOptions=gb.build(), 
            autoHeight=True,
            width='100%',
            allow_unsafe_jscode=True, #Set it to True to allow jsfunction to be injected
            enable_enterprise_modules=False, 
            fit_columns_on_grid_load=True,
            reload_data=False,
            checkboxSelection=True,
            key=key
            )

def build_ksnetwork_grid(data, key=None):
    '''
    This method builds the grid for kinase-substrate network relationships.

    :param data: This is the dataframe containing kinase-substrate network relationships
    :param key: This is the unique key for the grid
    :return: AgGrid object with the built grid
    '''
    #data.columns = ['Kinase(s)', 'Substrate(s)', 'Relationship Type', 'Related to', 'Kinase count', 'Substrate count']
    data.columns = ['Kinases', 'Substrates', 'Relationship Type', 'Related To', 'Related To ID', 'Kinase count', 'Substrate count']
    gb = GridOptionsBuilder.from_dataframe(data)
    
    gb.configure_pagination(
        enabled=True,
        paginationAutoPageSize=True
    )
    return AgGrid(
        data,
        gridOptions=gb.build(),
        autoHeight=True,
        width='100%',
        allow_unsafe_jscode=True,
        enable_enterprise_modules=False,
        fit_columns_on_grid_load=True,
        reload_data=False,
        checkboxSelection=False,
        key=key
    )

def get_ks_network_legend(network_relationship_types):
    '''
    This method generates the HTML legend for the kinase-substrate network graph.
    
    :param network_relationship_types: This is the set of relationship types in the network
    :return: legend_html: This is the HTML string for the legend
    '''
    legend_html = """<b>Node Legend:</b><br><div style="display: flex; gap: 6px;">"""
    legend_html += f"{util.get_circle(color=util.node_category_colors['Kinase Protein'])} : Enriched Kinases;&nbsp;&nbsp;&nbsp;"
    legend_html += f"{util.get_circle(color=util.node_category_colors['Substrate Protein'])} : Substrate Proteins;&nbsp;&nbsp;&nbsp;"
    for relationship_type in network_relationship_types:
        circle_html = util.get_circle(color=util.node_category_colors[relationship_type])
        legend_html += f"{circle_html} : {relationship_type};&nbsp;&nbsp;&nbsp;"
    
    legend_html += '</div>'
    return legend_html