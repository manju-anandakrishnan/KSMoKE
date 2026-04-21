import pandas as pd
import streamlit as st

from neo4j import GraphDatabase
from neo4j.exceptions import ClientError
import os

class DataAccessError(Exception):
    pass

class DatabaseAccess:

    '''
    This class establishes connection to the database and provides a driver instance.
    '''
    neo4j_driver = None

    @classmethod
    def get_driver(cls):
        '''
        This method creates a new database driver and verifies connectivity, if one doesn't exist. Else, it returns the existing driver object.
        '''
        if cls.neo4j_driver is None:
            cls.neo4j_driver = GraphDatabase.driver(os.getenv("DB_URI"), 
                                                    auth=(os.getenv("DB_USER"),os.getenv("DB_PWD")), 
                                                    database=os.getenv("DB_NAME"))
            cls.neo4j_driver.verify_connectivity()
        return cls.neo4j_driver
    
    @classmethod
    def close_driver(cls):
        '''
        This method closes the existing driver object
        '''
        if cls.neo4j_driver:
            cls.neo4j_driver.close()
            cls.neo4j_driver = None


class KinaseSubstrateDAO:

    '''
    This class provides data access methods for kinase-substrate relationships from the neo4j database.
    '''

    # This is the query template to fetch kinase-substrate paths
    KS_PATH_QUERY = """MATCH path=(k:Protein)-[*1..3]-(s:Protein)
            WHERE k.name='{k_name}' and s.name='{s_name}' 
                RETURN nodes(path) as n, relationships(path) as r limit 30
            """
    
    # This is the query template to fetch kinase-substrate network for multiple kinases and substrates
    # Hop is restricted to 2 for performance reasons
    KS_NETWORK_QUERY = """MATCH path=(k:Protein)-[*1..2]-(s:Protein)
            WHERE k.id in {k_ids} and s.id = '{s_ids}'
                RETURN nodes(path) as n, relationships(path) as r limit 1000
            """
    
    KS_RELTYPE_NETWORK_QUERY = """MATCH path=(k:Protein)-[:{rel_type}*1..2]-(s:Protein)
            WHERE k.id in {k_ids} and s.id = '{s_ids}'
                RETURN nodes(path) as n, relationships(path) as r limit 1000
            """
    

    def __init__(self, db_driver) -> None:
        '''
        This is the initialization method for KinaseSubstrateDAO. It requires a database driver object to be instantiated.
        
        Attributes:
        db_driver: neo4j driver object
        '''
        self.driver = db_driver

    def query_paths(self,tx,query):
        '''
        This method executes the provided query in a read transaction and returns the parsed triples.
        
        :param self: instance of the class
        :param tx: transaction object
        :param query: query string to be executed
        :return: triples: set of Triple objects representing kinase-substrate relationships
        '''
        nodes = []
        relationships = []
        for record in tx.run(query,timeout=1000):
            nodes.extend(record.get('n'))
            relationships.extend(record.get('r'))   
        result_parser = ResultParser(nodes, relationships)
        triples = result_parser.get_triples()
        return triples     


    def get_ks_links(self, kinase, substrate):
        '''
        This method fetches kinase-substrate links from the database for a given kinase and substrate.

        :param self: instance of the class
        :param kinase: kinase protein uniprotID
        :param substrate: substrate protein uniprotID
        :return: ks_triples: set of Triple objects representing kinase-substrate relationships
        '''
        query = KinaseSubstrateDAO.KS_PATH_QUERY
        query = query.format(k_name=kinase,s_name=substrate)
        ks_triples = set()
        with self.driver.session() as session:
            try:
                ks_triples = session.execute_read(self.query_paths, query)
            except ClientError as e:
                raise DataAccessError("Neo4j query failed: {e}") from e
        return ks_triples
    
    def get_ks_network(self, kinase_list, substrate_list, relationship_type=None):
        '''
        This method fetches kinase-substrate network from the database for given lists of kinases and substrates.

        :param self: instance of the class
        :param kinase_list: list of kinase protein uniprotIDs
        :param substrate_list: list of substrate protein uniprotIDs
        :return: ks_triples: set of Triple objects representing kinase-substrate relationships
        '''
        if relationship_type:
            query = KinaseSubstrateDAO.KS_RELTYPE_NETWORK_QUERY
            query = query.format(k_ids=kinase_list,s_ids=substrate_list,rel_type=relationship_type)
        else:
            query = KinaseSubstrateDAO.KS_NETWORK_QUERY
            query = query.format(k_ids=kinase_list,s_ids=substrate_list)

        ks_triples = set()
        with self.driver.session() as session:
            try:    
                ks_triples = session.execute_read(self.query_paths, query)
            except ClientError as e:
                raise DataAccessError("Neo4j query failed: {e}") from e
        return ks_triples


class ResultParser:
    '''
    This class parses the neo4j query results into Triple, Node and Edge objects.
    '''

    def __init__(self,nodes, relationships):
        '''
        This is the initialization method for ResultParser.

        :param self: instance of the class
        :param nodes: list of node objects returned by the query
        :param relationships: list of relationship objects returned by the query
        '''
        self.nodes = nodes
        self.relationships = relationships
        self.triples = set()

    def _parse_to_triples(self):
        '''
        This method parses the nodes and relationships into Triple objects.

        :param self: instance of the class
        '''
        node_dict = {}
        triples = set()
        
        for node in self.nodes:
            node_obj = Node(node['id'],node['name'],next(iter(node.labels)))
            node_dict[node.element_id] = node_obj #key is startNodeElementId

        for relationship in self.relationships:
            start_node = relationship.nodes[0]
            end_node = relationship.nodes[1]
            rel_obj = Edge(relationship.element_id,relationship.type)
            triple = Triple(node_dict[start_node.element_id],rel_obj,node_dict[end_node.element_id])
            triples.add(triple)
        
        return triples

    def get_triples(self):
        if not self.triples:
            self.triples = self._parse_to_triples()
        return self.triples


class Triple:
    '''
    This class represents a triple consisting of head node, edge and tail node.
    '''

    def __init__(self, head, edge, tail):
        '''
        This is the initialization method for Triple.

        :param self: instance of the class
        :param head: head node object
        :param edge: edge object
        :param tail: tail node object
        '''
        self.head = head
        self.rel = edge
        self.tail = tail

    def get(self):
        '''
        This method returns the head, edge and tail of the triple.

        :param self: instance of the class
        '''
        return self.head, self.rel, self.tail
    
    def __eq__(self,other):
        '''
        This method checks equality between two Triple objects.

        :param self: instance of the class
        :param other: other Triple object to compare with
        '''
        if isinstance(other,Triple):
            return (self.head == other.head)  \
                    & (self.tail == other.tail) \
                    & (self.rel == other.rel)
        return False

    def __hash__(self):
        '''
        This method returns the hash value of the Triple object.

        :param self: instance of the class
        '''
        return hash((self.head,self.rel,self.tail))
    
class Node:

    '''
    This class represents a node in the graph.
    '''

    def __init__(self,id,name,category):
        '''
        This is the initialization method for Node.

        :param self: instance of the class
        :param id: id of the node
        :param name: name of the node
        :param category: category of the node
        '''
        self.id = id
        self.name = name
        self.category = category

    def get_id(self):
        return self.id
    
    def get_name(self):
        return self.name
    
    def get_category(self):
        return self.category
    
    def __eq__(self,other):
        if isinstance(other,Node):
            return self.id == other.id
        return False
    
    def __hash__(self):
        return hash(self.id)

    
class Edge:
    '''
    This class represents an edge in the graph.
    '''

    def __init__(self,id,label):
        '''
        This is the initialization method for Edge.

        :param self: instance of the class
        :param id: id of the edge
        :param label: label of the edge
        '''
        self.id = id
        self.label = label

    def get_id(self):
        return self.id
    
    def get_label(self):
        return self.label
    
    def __eq__(self,other):
        if isinstance(other,Edge):
            return self.id == other.id
        return False
    
    def __hash__(self):
        return hash(self.id)


# This method is cached to improve performance
@st.cache_data
def get_kg_substrates(cut_off_prob):
    '''
    This method returns all kg substrates as a list of tuples. 
    Each tuple represents a substrate protein, the phosphorylation site and the -/+4 motif
    :param cut_off_prob: This is the cut-off probability threshold for including a substrate
    :return: List of tuples of (substrate_gene, substrate_protein, site, motif)
    '''
    kg_substrates_df = pd.read_csv('data/substrates_motif.csv', sep='|')
    kg_substrates_df = kg_substrates_df[kg_substrates_df['max_score'] >= cut_off_prob].copy().reset_index(drop=True)
    substrate_genes, substrate_proteins,sites,motifs = kg_substrates_df['gene'].to_list(), kg_substrates_df['Seq_Substrate'].to_list(), kg_substrates_df['Site'].to_list(), kg_substrates_df['Motif'].to_list()
    return [(substrate_genes[idx],substrate_proteins[idx],sites[idx],motifs[idx]) for idx in range(0,kg_substrates_df.shape[0])]
