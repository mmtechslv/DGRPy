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


import pandas as pd
import numpy as np
import itertools
import copy

from .jadval_loader import *

class OtuTable(JadvalPrime):
    """OtuTable is responsible for loading, forming and maintaining, and manipulating OTU/Taxa data and sample metadata
    This class can load classic OTU tables, separate sample metadata or BIOM format files.
    """
    def __init__(self,iStudyName,iStudyInfo = None):
        self.study_name = iStudyName #String name of this study
        self.study_info = iStudyInfo #Must be given as dictionary that contains information about study/research/data source. Such as author name, article title, publication year, date when data was collected, NGS device information, etc. This information will be used when combining jadval objects and printed into reports during report generation.
        self.otu_table = None #Pandas dataframe only with OTU reads. Dataframe has rows as OTU ids and columns as samples with sample ids
        self.taxon_groups = [] #List with dicts of Group ID-to-OTU ID associations and vice versa
        self.taxonomy = pd.DataFrame(columns=(['lineage']+jConst.MAIN_RANKS)) #Dataframe with OTU taxonomy with rows as OTU ids and columns with various formatted taxonomic information
        self.sample_metadata = None #Pandas dataframe for sample metadata with rows as samples ids and columns for datasets
        self.taxonomic_db = [] # List of loaded taxonomic databases
        self.default_tdb_id = None # List ID of taxonomic database that is present in self.taxonomic_db and will be used as default taxonomic database
        self.label_assoc = {} #Dict with otu_id - label associations
        self.__loader = {'otu_f_path':None,'taxa_f_path':None,'meta_f_path':None,'lineages':None} #Dict with original loaded data
        self.lost_otus = [] #List of Dict with OTUs that was lost
        self._init_state = [False,False,False] #List with init_state values. Order: [OTU Table,Taxonomy,Metadata]

    #Following method loads classic OTU table. Classic OTU table must have following structure. First row must contain headers. Columns: 1st for OTU id labels; 2nd to last column for OTU reads, and if specified last column for taxonomy consensus lineages
    @jErrorHandler(15)
    def load_OTU_table(self, iFilename, iSeparator=',', iConsensusLineage=True, iReversed= False): #iFilename - must be entered; iSeparators - characters that separates data, by default is comma; iQuote - character used for quoting in csv file, by default character is <">; iConsensusLineage - if last column contains consensus lineages for every OTU, default is True; iReversed - if table is reversed, by default OTUs rows of first column and samples are amoung columns of first row
        pre_data  = self.read_csv(iFilename,iSeparator)
        last_col = (len(pre_data[0])-1) if iConsensusLineage else (len(pre_data[0]))
        pre_sample_ids = pre_data[0][1:(last_col)]
        pre_otu_labels = [str(elem[0]) for elem in pre_data[1:]]
        pre_otu_ids = list(range(len(pre_otu_labels)))
        pre_otu_table = [elem[1:(last_col)] for elem in pre_data[1:]]
        self.label_assoc = {'id':pd.Series(data=pre_otu_labels,index=pre_otu_ids),'label':pd.Series(data=pre_otu_ids,index=pre_otu_labels)}
        self.otu_table = pd.DataFrame(data=pre_otu_table,index=pre_otu_ids,columns=pre_sample_ids)
        self.__loader['otu_f_path'] = iFilename
        self._init_state[0] = True
        if iConsensusLineage:
            pre_otu_lineages = [elem[last_col] for elem in pre_data[1:]]
            self.__loader['lineages'] = pd.Series(data=pre_otu_lineages, index=pre_otu_labels)
            self.taxonomy = pd.DataFrame(index=pre_otu_ids,columns=self.taxonomy.columns)
            self.__construct_taxonomy_from_lineages(pd.Series(data=pre_otu_lineages,index=pre_otu_ids))
            self.reconstruct_internal_lineages()
            self.__loader['taxa_f_path'] = iFilename
            self._init_state[1] = True
        return

    #This method checks state of initiation of OtuTable. state_i can be 1,2,3 or combinations like [1,2] or [2,3] or [1,2,3]. The last one is equivalent to simply self.check_init()
    @jErrorHandler(17)
    def check_init(self,state_i=False):
        if state_i:
            state_i = jShared.ensure_list(state_i)
            state = True
            for i in state_i:
                state = state and self._init_state[i-1]
            return state
        else:
            return all(self._init_state)

    #This method loads lineages from separate CSV file with first column as OTU id labels and second as lineages
    @jErrorHandler(18)
    def load_OTU_lineages(self,iFilename,iSeparator=',',iHeader=False):
        if self.check_init(1):
            pre_data = self.read_csv(iFilename, iSeparator)
            pre_otu_labels = [str(elem[0]) for elem in pre_data[(1 if iHeader else 0):]]
            pre_otu_lineages = [elem[1] for elem in pre_data[(1 if iHeader else 0):]]
            pre_otu_ids = list(self.label_assoc['label'].loc[pre_otu_labels])
            self.__loader['lineages'] = pd.Series(data=pre_otu_lineages, index=pre_otu_labels)
            self.taxonomy = pd.DataFrame(index=pre_otu_ids, columns=self.taxonomy.columns)
            self.__construct_taxonomy_from_lineages(pd.Series(data=pre_otu_lineages, index=pre_otu_ids))
            self.reconstruct_internal_lineages()
            self.__loader['taxa_f_path'] = iFilename

            self._init_state[1] = True
        else:
            return jRaiseError(16)
        return

    #This function accepts single or many OTU labels and returns corresponding taxonomy as dictionaries.
    @jErrorHandler(19)
    def get_taxonomy_by_otu_label(self,iOTU_labels):
        otu_ids = list(self.label_assoc['label'].loc[iOTU_labels])
        return self.get_taxonomy_by_otu_id(otu_ids)
    
    #This function accepts single or many OTU ids and returns corresponding taxonomy as dictionaries.
    @jErrorHandler(20)
    def get_taxonomy_by_otu_id(self,iOTU_ids):
        iOTU_ids = jShared.ensure_list(iOTU_ids) #IMPROVE!!! Check if iOTU_ids is a list and converts it into list if it is not
        all_taxa_dict = {} #Dictionary with multiple taxonomy
        for otu_id in iOTU_ids:
            taxa_dict = dict(self.taxonomy.loc[otu_id,jConst.MAIN_RANKS])
            #Following sets single letter rank dict keys to full ones that are defined in jConst.ITS
            taxa_dict = {jConst.ITS['r2Rank'][k]:v for k,v in taxa_dict.items()}
            all_taxa_dict[otu_id] = taxa_dict # Adds taxonomy of current OTU into dict with multiple taxonomies
        return all_taxa_dict if (len(iOTU_ids)>1) else taxa_dict #Returns multiple taxonomy dict if more than one OTU id was requested.
    
    #This method returns list of available rank levels. Method checks highest and lowest rank available in self.taxonomy that is also present among jConst.MAIN_RANKS
    @jErrorHandler(9)
    def get_avail_ranks(self):
        return [rank for rank in jConst.MAIN_RANKS if self.taxonomy.loc[:,rank].notna().any()]
    
    #This method reconstructs self.taxonomy.loc[:,'lineage'] by initiated taxa
    @jErrorHandler(21)
    def reconstruct_internal_lineages(self):
        self.taxonomy.loc[:,'lineage'] = self.generate_lineages(iMissingRank=True,iDesiredRanks=self.get_avail_ranks(),iDropRanks=False)
        return

    #This method removes OTUs based on ids
    @jErrorHandler(22)
    def drop_otus_by_ids(self,iOTU_ids):
        iOTU_ids  = jShared.ensure_list(iOTU_ids)
        all_ids = self.get_ids()
        if all(map(lambda id: id in all_ids,iOTU_ids)):
            return self.__drop_taxa(iOTU_ids,'manual-by-id')
        else:
            return False
        
    #This method removes OTUs based on labels
    @jErrorHandler(23)
    def drop_otus_by_labels(self,iOTU_labels):
        iOTU_labels  = jShared.ensure_list(iOTU_labels)
        all_labels = self.get_labels()
        if all(map(lambda id: id in all_labels,iOTU_labels)):
            iOTU_ids = list(self.label_assoc['label'][iOTU_labels].values)
            return self.__drop_taxa(iOTU_ids,'manual-by-label')
        else:
            return False

    #This method removes OTUs that do not have any assigned taxonomy from self.taxonomy dataframe and returns OTUs that was removed
    @jErrorHandler(24)
    def drop_otus_without_taxa(self):
        otu_indices_drop = self.taxonomy.loc[self.taxonomy.loc[:,jConst.MAIN_RANKS].agg(lambda rank: len(''.join(map(lambda x: (str(x or '')),rank))),axis=1) < 1].index
        removed_otus = self.__drop_taxa(otu_indices_drop,'drop-otus-without-taxa')
        return removed_otus
    
    #This method drops otus that do not have taxon at given rank or ranks in iRanks. If iRanks contains more than one rank and iMerge = True then OTUs with any missing rank given in iRanks will be removed. Otherwise, if iMerge = False then for OTU to be removed, all ranks given in iRanks must be missing
    @jErrorHandler(25)
    def drop_otus_without_ranks(self,iRanks,iMerge=False):
        selected_otus = self.taxonomy.loc[:,(iRanks if type(iRanks)==list else [iRanks])].isna()
        selected_otus_merged = selected_otus.any(axis=1) if iMerge else selected_otus.all(axis=1)
        if selected_otus_merged.any():
            otus_without_ranks = self.taxonomy.loc[selected_otus_merged]
            otu_indices_drop = otus_without_ranks.index
            removed_otus = self.__drop_taxa(otu_indices_drop,('drop-otus-without-ranks (%s)%s' %(','.join(iRanks),(' + merge' if iMerge else ''))))
        else:
            removed_otus = False
        return removed_otus
    
    #Merge OTU reads by duplicated lineages
    @jErrorHandler(26)
    def merge_duplicated_taxa(self):
        groups = self.taxonomy.groupby('lineage')
        if len(groups.groups)>1:
            tmp_otu_table  = []
            tmp_otu_lineage  = []
            tmp_groups = []
            group_indices = list(range(len(groups.groups)))
            for lineage,otu_ids in groups.groups.items():
                tmp_otu_table.append(self.otu_table.loc[otu_ids,:].agg(lambda x: sum(map(int,x))))
                tmp_otu_lineage.append(lineage)
                tmp_groups.append(list(otu_ids))
            self.otu_table = pd.DataFrame(data=tmp_otu_table,index=group_indices,columns=self.otu_table.columns)
            self.__construct_taxonomy_from_lineages(pd.Series(data=tmp_otu_lineage,index=group_indices))
            self.__record_group_relations(tmp_groups,group_indices,'merge-duplicates')
        return
    
    #This method groups and merges OTUs based on ranks up-to(including) iRank.
    @jErrorHandler(27)
    def merge_otus_by_rank(self,iRank):
        ranks = jConst.MAIN_RANKS[:jConst.MAIN_RANKS.index(iRank)+1]
        missing_ranks = len(jConst.MAIN_RANKS[jConst.MAIN_RANKS.index(iRank)+1:])
        groups = self.taxonomy.groupby(ranks)
        if len(groups.groups)>1:
            tmp_otu_table  = []
            tmp_taxonomy  = []
            tmp_groups = []
            group_indices = list(range(len(groups.groups)))
            for taxa,otu_ids in groups.groups.items():
                tmp_otu_table.append(self.otu_table.loc[otu_ids,:].agg(lambda x: sum(map(int,x))))
                tmp_taxonomy.append([None]+([rank if pd.notna(rank) else None for rank in taxa] if type(taxa)==tuple else [taxa])+[None]*missing_ranks)
                tmp_groups.append(list(otu_ids))
            self.otu_table = pd.DataFrame(data=tmp_otu_table,index=group_indices,columns=self.otu_table.columns)
            self.taxonomy = pd.DataFrame(data=tmp_taxonomy,index=group_indices,columns=self.taxonomy.columns)
            self.__record_group_relations(tmp_groups,group_indices,("merge-by-rank(%s)" %iRank))
            self.reconstruct_internal_lineages()
        else:
            return False
        return tmp_groups
    
    #This function generate desired lineages
    @jErrorHandler(1)
    def generate_lineages(self,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False 
        return self.generate_lineages_from_taxa(iTaxa=self.taxonomy,iMissingRank=iMissingRank,iDesiredRanks=iDesiredRanks,iDropRanks=iDropRanks)
    
    #This method loads JadvalDB and stores truncated database taxonomy and tree for further usage
    @jErrorHandler(28)
    def load_database(self,iDerivationName,iJadvalDB_obj, iRankTolerance=False):
        correlation_result = self.__derive_taxonomy_db(iJadvalDB_obj, iRankTolerance)
        if not correlation_result:
            return False
        pruned_tree = iJadvalDB_obj.prune_tree_for_taxa_ids(correlation_result['taxon_ids'])
        DBType = iJadvalDB_obj.__class__.__name__
        derivation_record = {'name':iDerivationName,'type':DBType,'correlations':correlation_result['corr_as_series'],'tree':pruned_tree,'details':correlation_result}
        self.taxonomic_db.append(derivation_record)
        return derivation_record

    #This function correlate taxon ids between available lineages and iLineages. If iTolRanks is set then correlation of uncorrelated taxa will continue until all ranks are assigned or iTolRanks are exausted.
    @jErrorHandler(29)
    def generate_taxa_correlations(self,iLineages,iForRanks,iTolRanks=False):
        self_lineages = self.generate_lineages(True,iForRanks)    
        correlation_results = self.correlate_lineages(self_lineages,iLineages)
        correlations = correlation_results['correlations']
        uncorrelated_lineages = correlation_results['unassigned']
        uncorrelated_otu_ids = uncorrelated_lineages.index
        correlated_lineages = correlation_results['assigned']
        if iTolRanks:
            rank_trials = {'init':{'corr':correlations.to_dict(),'uncorr_ids':uncorrelated_otu_ids,'added':sum(correlations.notna()),'lineages':correlation_results['assigned']}}
            for limiting_rank in [rank for rank in jConst.MAIN_RANKS[::-1] if rank in iTolRanks.keys()]:
                self_lineages = self.generate_lineages(True,iForRanks,iTolRanks[limiting_rank])
                unassigned_self_lineages = self_lineages.loc[uncorrelated_otu_ids]
                correlation_results = self.correlate_lineages(unassigned_self_lineages,iLineages)
                new_correlations = correlation_results['correlations']
                correlations.loc[uncorrelated_otu_ids] = new_correlations
                uncorrelated_lineages = correlation_results['unassigned']
                uncorrelated_otu_ids = uncorrelated_lineages.index
                correlated_lineages.append(correlation_results['assigned'])
                rank_trials[limiting_rank] = {'corr':new_correlations.to_dict(),'uncorr_ids':uncorrelated_otu_ids,'added':sum(new_correlations.notna()),'lineages':correlation_results['assigned'],'tolerated_ranks':iTolRanks[limiting_rank]}
        return {'lineages':correlated_lineages, 'corr_as_series':correlations,'corr':correlations.to_dict(),'taxon_ids':set(correlations[correlations.notna()].values),'otu_ids':list(correlations[correlations.notna()].keys()),'total':sum(correlations.notna()),'rank_trials':rank_trials if iTolRanks else None}

    #This method resets self and reconstructs it so that only OTUs with given ids are present
    @jErrorHandler(30)
    def reconstruct_for_otu_ids(self,iOTU_ids):
        iOTU_ids = jShared.ensure_list(iOTU_ids)
        indices = self.taxonomy.index
        to_remove = indices[~indices.isin(iOTU_ids)]
        self.__drop_taxa(to_remove, 'internal-reconstruct')
        self.taxon_groups = []
        self.taxonomic_db = []
        self.default_tdb_id = None
        return

    #This method retrieves new OtuTable with only OTUs/Taxa correlatated with iTaxonomy taxonomic database
    @jErrorHandler(31)
    def retrieve_custom(self,iTaxonomy = None):
        if iTaxonomy is not None:
            if jShared.assert_taxonomic_db(iTaxonomy):
                c_self = copy.deepcopy(self)
                iOTU_ids = list(iTaxonomy['correlations'][iTaxonomy['correlations'].notna()].index)
                c_self.reconstruct_for_otu_ids(iOTU_ids)
                c_self.taxonomic_db = [iTaxonomy]
                c_self.default_tdb_id = 0
                c_self.__set_internal_ids(iTaxonomy['correlations'][iTaxonomy['correlations'].notna()])
                return c_self
        return False


    #This method gets active local taxonomic database that is defined in self.default_tdb_id
    @jErrorHandler(32)
    def get_active_taxonomic_db(self):
        if (self.default_tdb_id is not None):
            if self.default_tdb_id < len(self.taxonomic_db):
                tdb =  self.taxonomic_db[self.default_tdb_id]
                if jShared.assert_taxonomic_db(tdb):
                    return tdb
        return False

    #This method saves OTU table into CSV file. If iPreffix is set then OTU ids will start with preffix_ID
    @jErrorHandler(33)
    def save_as_csv(self,iFilename,iPreffix=False,):
        headers = ['ID'] + list(self.otu_table.columns) + ['Consensus Lineages']
        s_data = pd.concat([self.otu_table,self.taxonomy],axis=1)
        data = list(map(lambda x:list(x),list(s_data.loc[:,list(self.otu_table.columns)+['lineage']].values)))
        otu_ids = list(map(lambda x: iPreffix + '_' + str(x) ,list(s_data.index.values))) if iPreffix else list(s_data.index.values)
        content = [headers]+list(map(lambda x:[x[0]]+x[1],zip(otu_ids,data)))
        self.write_csv(content,iFilename)
        return

    #This method exports both OTU table and Tree file. This method must be fixed for sample data.
    #CHANGE TAXAN IDS WITH LINEAGES IN PDF RENDER
    @jErrorHandler(34)
    def export_data(self,iFilename=False):
        iFilename = iFilename if iFilename else self.study_name
        dir_path = self.ensure_new_dir(iFilename)
        if dir_path:
            table_filename = dir_path+'/'+iFilename + '_OTU_TABLE.csv'
            tree_filename = dir_path+'/'+iFilename + '_TREE.tree'
            tree_pdf_filename = dir_path+'/'+iFilename + '_TREE_VIEW.pdf'
            self.save_as_csv(table_filename)
            tdb = self.get_active_taxonomic_db()
            tdb['tree'].write(outfile=tree_filename, format=1)
            tdb['tree'].render(tree_pdf_filename)
        return

    #This method makes a copy of self
    @jErrorHandler(53)
    def copy(self):
        return copy.deepcopy(self)

    #This method finds OTU ids that contain given pattern. iCase = True for case sensitive, iRegex = True if Pattern is regular expression (otherwise it is literal string)
    @jErrorHandler(54)
    def search(self,iPattern,iCase=False,iRegex=False):
        return list(self.taxonomy[self.taxonomy.loc[:, 'lineage'].str.contains('Wolbachia',case=iCase,regex=iRegex)].index)

    #This method get list of internal OTU ids
    @jErrorHandler(55)
    def get_ids(self):
        return list(self.taxonomy.index)

    # This method get list of OTU labels
    @jErrorHandler(56)
    def get_labels(self):
        return list(itertools.chain.from_iterable(self.label_assoc['id']))


#>>>>>>>>>>>>>>>>>>>Private Methods<<<<<<<<<<<<<<<<<<<

    #This method changes internal OTU indices. THIS METHOD NEED FIX FOR GROUP IDS
    @jErrorHandler(35)
    def __set_internal_ids(self,new_id_assoc):
        self.otu_table.index = new_id_assoc
        self.taxonomy.index = new_id_assoc
        self.label_assoc['id'].index = new_id_assoc
        self.label_assoc['label'] = self.label_assoc['label'].replace(new_id_assoc)
        #GROUPS ARE MISSING

    #This function initiates self.taxonomy by breaking down lineages and storing ranks separately for further use
    @jErrorHandler(36)
    def __construct_taxonomy_from_lineages(self,lineages):
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
        self.taxonomy = pd.DataFrame(index=list(lineages.index),data=list(allocater_vectorized(lineage=list(lineages.values),ref_ranks=jConst.MAIN_RANKS)),columns=self.taxonomy.columns)
        self.__fix_missing_taxa()
        return

    #This method loads taxonomy database(JadvalDB) and correlates internal OTU ids with database taxon ids with rank tolerance given in iRank
    @jErrorHandler(37)
    def __derive_taxonomy_db(self,iJadvalDB, iRankTolerance=False):
        shared_avail_ranks = [rank for rank in self.get_avail_ranks() if rank in iJadvalDB.get_avail_ranks()]
        db_lineages = iJadvalDB.generate_lineages(True,shared_avail_ranks)
        #Here rank lists are generated until it reaches iRank. For example, if iRankTolerance='f' and shared_avail_ranks=['c', 'o', 'f', 'g', 's'] then iTolRanks will contain {'f': ['g', 's'], 'g': ['s']}
        if not iRankTolerance:
           iTolRanks = False
        else:
            if (iRankTolerance in shared_avail_ranks) and ((shared_avail_ranks.index(iRankTolerance)+1)!=len(shared_avail_ranks)):
                tol_index = len(shared_avail_ranks) - shared_avail_ranks.index(iRankTolerance)
                iTolRanks = {shared_avail_ranks[-(e+1)]:shared_avail_ranks[-e:] for e in range(1,tol_index)}
            else:
                return False
        correlation_result = self.generate_taxa_correlations(db_lineages,shared_avail_ranks,iTolRanks)
        return correlation_result

    #This function generates taxon_groups dict. taxon_groups['group'] - contains group ids and corresponding otu ids; taxon_groups['otus'] - contains otu ids and corresponding group ids
    @jErrorHandler(38)
    def __record_group_relations(self,groups,group_indices,action_name):
        otu_ids = []
        group_ids = []
        for group_id in group_indices:
            for otu_id in groups[group_id]:
                group_ids.append(group_id)
                otu_ids.append(otu_id)
        group_relations = {'new-prior':pd.Series(data=groups,index=group_indices),'prior-new':pd.Series(data=group_ids,index=otu_ids)}
        self.taxon_groups.append({action_name:group_relations})
        self.__reassign_label_assoc(group_relations)
        return

    #This method reconstructs self.label_assoc after grouping
    @jErrorHandler(39)
    def __reassign_label_assoc(self,assign_group):
        def gen_labels(taxons):
            labels = []
            for tmp in taxons:
                taxon = self.label_assoc['id'][tmp]
                if type(taxon)==list:
                    for in_taxon in taxon:
                        labels.append(in_taxon)
                else:
                    labels.append(taxon)
            return labels
        new_id_to_labels = assign_group['new-prior'].apply(gen_labels)
        new_label_to_id = self.label_assoc['label'].apply(lambda x: assign_group['prior-new'][x])
        self.label_assoc = {'id':new_id_to_labels,'label':new_label_to_id}
        return

    # This method records taxa that is to be removed into self.lost_otus
    @jErrorHandler(40)
    def __record_removed_taxa(self, otu_ids, action_name):
        removed_otus = pd.concat([self.otu_table.loc[otu_ids, :], self.taxonomy.loc[otu_ids, :]], axis=1)
        removed_otus.index = self.label_assoc['id'][otu_ids]
        self.lost_otus.append({action_name: removed_otus})
        return removed_otus

    # This method removes given otu ids from self.label_assoc
    @jErrorHandler(41)
    def __remove_labels_by_ids(self, ids):
        labels_to_remove = self.label_assoc['id'][ids] if type(self.label_assoc['id'][ids].iloc[0]) != list else list(
            itertools.chain.from_iterable(self.label_assoc['id'][ids]))
        self.label_assoc['id'] = self.label_assoc['id'].drop(ids)
        self.label_assoc['label'] = self.label_assoc['label'].drop(labels_to_remove)
        return labels_to_remove

    # This is main method that removes OTU from this instance by recording action into self.lost_otus; removing from self.label_assoc and droping otu from self.taxonomy and self.otu_table
    @jErrorHandler(42)
    def __drop_taxa(self, ids, action_name):
        removed_otus = self.__record_removed_taxa(ids, action_name)
        self.__remove_labels_by_ids(ids)
        self.taxonomy = self.taxonomy.drop(ids)
        self.otu_table = self.otu_table.drop(ids)
        return removed_otus

    #This method sets taxons with "" = None. For example, lineages sometimes contain ranks such as (...; g__; s__) and (...; g__;). If consensus lineages are only available taxonomic information then such taxa are basically same and must be fixed.
    @jErrorHandler(10)
    def __fix_missing_taxa(self):
        self.taxonomy.loc[:,jConst.MAIN_RANKS] = self.taxonomy.loc[:,jConst.MAIN_RANKS].applymap(lambda x: None if (x=='') else x)
        return

#>>>>>>>>>>>>>>>>>>Properties<<<<<<<<<<<<<<<<<<

    @property
    @jErrorHandler(44)
    def jadval(self):
        return pd.concat([self.otu_table,self.taxonomy],axis=1)

    @property
    @jErrorHandler(45)
    def lineages(self):
        return self.taxonomy.loc[:,'lineage']
    
    @property
    @jErrorHandler(46)
    def LabelById(self):
        return self.label_assoc['id']
    
    @property
    @jErrorHandler(47)
    def IdByLabel(self):
        return self.label_assoc['label']

    @property
    @jErrorHandler(48)
    def loader_info(self):
        return self.__loader

    @property
    @jErrorHandler(49)
    def phylo_tree(self):
        tdb = self.get_active_taxonomic_db()
        return tdb['tree'] if tdb else jRaiseWarning(1,'Phylogenetic Tree is Not Found.')

    @property
    @jErrorHandler(50)
    def taxonomy_db_type(self):
        tdb = self.get_active_taxonomic_db()
        return tdb['type'] if tdb else jRaiseWarning(1,'Database Type N/A.')

    @property
    @jErrorHandler(51)
    def taxonomy_db_name(self):
        tdb = self.get_active_taxonomic_db()
        return tdb['name'] if tdb else jRaiseWarning(1,'Database Name is N/A.')

    ###################### NOT FINISHED


    #Following method loads sample metadata with rows as sample names and columns as various data sets.
    def load_sample_metadata(self,iFilename,iSeparator=',', iQuote='"', iHeader=True,iReversed = False): #iFilename - must be entered; iSeparators - characters that separates data, by default is comma; iQuote - character used for quoting in csv file, by default character is <">; iHeader - True if first row is headers, otherwise datasets will be named automatically as SMD_1, SMD_2, etc; iReversed - if table is reversed, by default samples names are given in the first row and datasets are given in next columns 
        self.file_meta_filename = iFilename
