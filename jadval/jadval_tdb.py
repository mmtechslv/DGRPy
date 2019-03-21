#!/usr/bin/env python3

"""
DGRPy is a python package for DGRP microbiome data analysis
The project is under early development and this is my initial release.
NOTE: This is an unstable early release so as soon as I publish stable version I will also add documentation.
Many functions are not at final form at all. For example, JadvalOTU can only load OTU data only in format as OTUs_sample.csv

@author: Farid MUSA(mmtechslv)
"""

__author__ = "Farid MUSA"
__copyright__ = "Copyright (C) 2019, DGRPy Project"
__credits__ = ["Farid MUSA"]
__license__ = "GPLv3"
__version__ = "1.1"
__maintainer__ = "Farid MUSA"
__email__ = "farid.musa.h@gmail.com"
__status__ = "Development"


import copy
import ete3
import pandas as pd
import numpy as np

from .jadval_loader import *

class TaxonomyMain(JadvalPrime):
    """TaxonomyMain DGRPy class contain methods and variables shared amoung TaxonomyGG, TaxonomySILVA, TaxonomyNCBI, TaxonomyOTL, TaxonomyRDP """
    def __init__(self):
        super().__init__()
        self.taxonomy_df = pd.DataFrame(columns=(['lineage']+jConst.MAIN_RANKS)) #Dataframe with OTU taxonomy with rows as OTU ids and columns with various formatted taxonomic information
        self.taxonomy_map = None
        self.tree = None
        self.loaded_lineages = None #Pandas Series with original lineages
        self.init_state = {'taxa':False,'tree':False}

    #Following method loads greengenes taxonomic database
    @jErrorHandler(7)
    def _load_tree(self,iTreeFile,iTreeFormat='newick'): #iTreeFile - Tree file in format given in iTreeFormat, by default is newick;
        with open(iTreeFile, 'r') as tree_file:
            content = tree_file.read()
        self.tree = ete3.Tree(content,quoted_node_names=True,format=1)
        return
    
    #This method prunes the database tree for given taxa and returns a copy
    @jErrorHandler(8)
    def prune_tree_for_taxa_ids(self,iTaxa_ids):
        if self.tree and (type(iTaxa_ids)==set):
            tmp_tree = copy.deepcopy(self.tree)
            tmp_tree.prune(iTaxa_ids, preserve_branch_length=True)
            return tmp_tree
        else:
            return False
    
    #This function generate desired lineages
    @jErrorHandler(1,'Overwritten by TaxonomyMain')
    def generate_lineages(self,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False 
        return self.generate_lineages_from_taxa(iTaxa=self.taxonomy_df,iMissingRank=iMissingRank,iDesiredRanks=iDesiredRanks,iDropRanks=iDropRanks)
    
    #This method returns list of available rank levels. Method checks highest and lowest rank available in self.taxonomy_df that is also present among jConst.MAIN_RANKS
    @jErrorHandler(9)
    def get_avail_ranks(self):
        return [rank for rank in jConst.MAIN_RANKS if self.taxonomy_df.loc[:,rank].notna().any()]
    
    #This method sets taxons with "" = None. For example, lineages sometimes contain ranks such as (...; g__; s__) and (...; g__;). If consensus lineages are only available taxonomic information then such taxa are basically same and must be fixed.
    @jErrorHandler(10)
    def _fix_missing_taxa(self):
        self.taxonomy_df.loc[:,jConst.MAIN_RANKS] = self.taxonomy_df.loc[:,jConst.MAIN_RANKS].applymap(lambda x: None if (x=='') else x)
        return
    
    #This method reconstructs self.taxonomy_df.loc[:,'lineage'] by initiated taxa
    @jErrorHandler(11)
    def reconstruct_internal_lineages(self):
        self.taxonomy_df.loc[:,'lineage'] = self.generate_lineages(iMissingRank=True,iDesiredRanks=self.get_avail_ranks(),iDropRanks=False)
        return 

class TaxonomyGG(TaxonomyMain):
    """TaxonomyGG DGRPy class that is responsible for loading, maintaining, and manipulating GreenGenes Taxonomy"""
    def __init__(self,iTaxaMapFile=None,iTreeFile=None,iTreeFormat='newick',iDBInfo = None,iStateFile=None):
        super().__init__()
        self.db_info = iDBInfo #Must be given as dictionary that contains information about study/research/data source. Such as author name, article title, publication year, date when data was collected, NGS device information, etc. This information will be used when combining jadval objects and printed into reports during report generation.
        if iStateFile:
            self.load_state(iStateFile)
        if iTreeFile and iTaxaMapFile:
            self.load_taxonomy_database(iTaxaMapFile,iTreeFile,iTreeFormat)
    
    #Following method loads taxonomic database
    @jErrorHandler(12)
    def load_taxonomy_database(self,iTaxaMapFile,iTreeFile,iTreeFormat='newick'):       
        self.__load_taxonomy_map(iTaxaMapFile)
        self._load_tree(iTreeFile,iTreeFormat)
        self.__init_internal_taxonomy()
        self.reconstruct_internal_lineages()
        return
    
    #Following method loads Greengenes taxonomic database
    @jErrorHandler(13)
    def __load_taxonomy_map(self,iTaxaMapFile): #iDatabaseName - Name of database; iTreeFile - Tree file in format given in iTreeFormat, by default is newick; iTaxaMapFile - file with ID - Taxonomy data, by default values are separated by \t
        tmp_taxa_map = self.read_csv(iTaxaMapFile,iSeparator='\t')
        self.taxonomy_map = pd.Series(data=[e[1] for e in tmp_taxa_map],index=[e[0] for e in tmp_taxa_map])
        return
    
    #This function initiates self.taxonomy_df by breaking down lineages and storing ranks separately (Vectorization is not an option. Reoptimization is necessary!)
    @jErrorHandler(14,'GreenGenes Initializer Failed')
    def __init_internal_taxonomy(self):
        def allocater(lineage,ref_ranks):
            taxa_dict = {e[0]:e[1] for e in [e.strip().split('__') for e in lineage.split(';') if ('__' in e)]} #Efficiently converts lineage into dictionary
            taxa_dict_allowed = {rank:taxa_dict[rank] for rank in taxa_dict.keys() if rank in ref_ranks} #Drops unallowed ranks and orders ranks based on jConst.MAIN_RANKS rank order
            #Following loop sets unavailable ranks to None in order to avoid later problems
            for key in ref_ranks:
                if not (key in taxa_dict_allowed.keys()):
                    taxa_dict_allowed[key] = None
            taxa_list_ordered = [taxa_dict_allowed[rank] for rank in ref_ranks] #Orders ranks based on jConst.MAIN_RANKS rank order
            return [lineage]+taxa_list_ordered
        allocater_vectorized = np.vectorize(allocater,excluded=['ref_ranks'],otypes=[list])
        self.taxonomy_df = pd.DataFrame(index=list(self.taxonomy_map.index),data=list(allocater_vectorized(lineage=list(self.taxonomy_map.values),ref_ranks=jConst.MAIN_RANKS)),columns=self.taxonomy_df.columns)
        self._fix_missing_taxa()
        return True

class TaxonomySILVA(TaxonomyMain):
    """TaxonomySILVA DGRPy class that is responsible for loading, maintaining, and manipulating SILVA Taxonomy"""
    def __init__(self,iTaxaMapFile=None,iTreeFile=None,iTreeFormat='newick',iDBInfo = None,iStateFile=None):
        super().__init__()
        self.db_info = iDBInfo #Must be given as dictionary that contains information about study/research/data source. Such as author name, article title, publication year, date when data was collected, NGS device information, etc. This information will be used when combining jadval objects and printed into reports during report generation.
        if iStateFile:
            self.load_state(iStateFile)
        if iTreeFile and iTaxaMapFile:
            self.load_taxonomy_database(iTaxaMapFile, iTreeFile,iTreeFormat)
    
    #Following method loads taxonomic database
    @jErrorHandler(12)
    def load_taxonomy_database(self, iTaxaMapFile, iTreeFile,iTreeFormat='newick'):       
        self.__load_taxonomy_map(iTaxaMapFile)
        self._load_tree(iTreeFile,iTreeFormat)
        self.__init_internal_taxonomy()   
        self.reconstruct_internal_lineages()
        return
    
    #Following method loads greengenes taxonomic database
    @jErrorHandler(13)
    def __load_taxonomy_map(self,iTaxaMapFile): #iTaxaMapFile - file with ID - Taxonomy data, by default values are separated by \t
        tmp_taxa_map = self.read_csv(iTaxaMapFile,iSeparator='\t')
        self.taxonomy_map = pd.DataFrame(data=[[e[0], e[2]] for e in tmp_taxa_map],index=[e[1] for e in tmp_taxa_map],columns=['lineage','level'])
        return
    
    #This function initiates self.taxonomy_df by breaking down lineages and storing ranks separately (This function is NOT OPTIMIZED. Cython can by considered)
    @jErrorHandler(14, 'SILVA Initializer Failed')
    def __init_internal_taxonomy(self):
        SILVA_RANKS = jConst.MAIN_RANKS[:-1]
        tmp_taxonomy_df = pd.DataFrame(columns=(['lineage']+SILVA_RANKS),index=self.taxonomy_map.index,data=None)
        rank = SILVA_RANKS[0]
        tmp_taxonomy_df.loc[:,rank] = self.taxonomy_map.loc[self.taxonomy_map.loc[:,'level']==jConst.ITS['r2rank'][rank],:].apply(lambda row: row['lineage'][:-1].split(';')[0],axis=1)
        for rank in SILVA_RANKS[1:]:    
            tmp_taxonomy_df.loc[:,rank] = self.taxonomy_map.loc[self.taxonomy_map.loc[:,'level']==jConst.ITS['r2rank'][rank],:].apply(lambda row: row['lineage'][:-1].split(';'),axis=1)
        def reassign(taxons,r):
            new_taxa = {r:taxons.pop()}
            for r_i in SILVA_RANKS[:SILVA_RANKS.index(r)]:
                for t_i in range(len(taxons)):
                    if taxons[t_i] in tmp_taxonomy_df.loc[tmp_taxonomy_df.loc[:,r_i].notna(),r_i].values:
                        new_taxa[r_i] = taxons[t_i]
                        taxons.pop(t_i)
                        break
            return pd.Series(new_taxa)
        for rank in SILVA_RANKS[1:]:    
            tmp_taxonomy_df.update(tmp_taxonomy_df.loc[tmp_taxonomy_df.loc[:,rank].notna(),rank].apply(reassign,r=rank))
        tmp_taxonomy_df.insert(loc=len(tmp_taxonomy_df.columns),column='s',value=None)
        tmp_taxonomy_df = tmp_taxonomy_df.applymap(lambda x: None if pd.isna(x) else x)
        self.taxonomy_df = tmp_taxonomy_df
        self._fix_missing_taxa()
        return
