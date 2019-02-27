#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DGRPy is a python package for DGRP microbiome data analysis
The project is under early development and this is my initial release.
NOTE: This is an unstable early release so as soon as I publish stable version I will also add documentation.
Many functions are not at final form at all. For example, JadvalOTU can only load OTU data only in format as OTUs_sample.csv

@author: Farid MUSA(mmtechslv)
"""

import pandas as pd
import numpy as np
import csv
import random
import datetime
import pickle
import itertools 
import ete3
import copy

__author__ = "Farid MUSA"
__copyright__ = "Copyright (C) 2019, DGRPy Project"
__credits__ = ["Farid MUSA"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Farid MUSA"
__email__ = "farid.musa.h@gmail.com"
__status__ = "Development"


###Following constants should not be modified
#JADVAL_ITS  links internal and external representation of values
JADVAL_ITS = {'lineage':'Consensus Lineage',
'r2Rank': {'d':'Domain','k':'Kingdom','p':'Phylum','c':'Class','o':'Order','f':'Family','g':'Genus','s':'Specie'},
'r2rank': {'d':'domain','k':'kingdom','p':'phylum','c':'class','o':'order','f':'family','g':'genus','s':'specie'},
'rank2r': {'domain':'d','kingdom':'k','phylum':'p','class':'c','order':'o','family':'f','genus':'g','specie':'s'},
'rank2Rank': {'domain':'Domain','kingdom':'Kingdom','phylum':'Phylum','class':'Class','order':'Order','family':'Family','genus':'Genus','specie':'Specie'},
'Rank2r': {'Domain':'d','Kingdom':'k','Phylum':'p','Class':'c','Order':'o','Family':'f','Genus':'g','Specie':'s'}}
JADVAL_MAIN_RANKS = ['d','k','p','c','o','f','g','s']  #JADVAL_MAIN_RANKS list of main taxonomic ranks used in Jadval. Do not modify. Order is important!


class JadvalMain:
    """Jadval is main class of DGRPy and parental class for all others. This class contains definitions, variables, constants and methods shared with other classes."""
    
    @staticmethod
    def generate_lineages_from_taxa(iTaxa,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False    
        if iDesiredRanks and not all(e in JADVAL_MAIN_RANKS for e in iDesiredRanks):
            print('Impossible characters found in iDesiredRanks. Please use: '+','.join(JADVAL_MAIN_RANKS))
            return False
        iDropRanks = iDropRanks if iDropRanks else []
        add_next = (lambda r,t: ((r+'__'+t+'; ') if t != None else r+'__'+'; ')) if iMissingRank else (lambda r,t: ((r+'__'+t+'; ') if t != None else '')) 
        def make_lineage(taxon):
            tmp_lineage = ''
            for rank in make_ranks:
                tmp_lineage = tmp_lineage + add_next(rank,taxon.loc[rank] if (rank not in iDropRanks) else None)   #Add succeeding rank to lineage
            return tmp_lineage[:-2]
        if not iDesiredRanks:
            make_ranks = JADVAL_MAIN_RANKS
            new_lineages = iTaxa.apply(make_lineage,axis=1)
        else:
            make_ranks = [rank for rank in JADVAL_MAIN_RANKS if rank in iDesiredRanks]
            new_lineages = iTaxa.apply(make_lineage,axis=1)
        return new_lineages
    
    #This function find common lineages between iPrimaryLineages and iSecondaryLineages.
    #Result contains dict with {'assigned'} - reserved for lineages that are common to both lineage series
    #{'unassigned'} - reserved for lineages that are not shared between two lineage series
    #{'correlations'} - contains correlation ids between two lineages
    #Both iPrimaryLineages and iSecondaryLineages must be pandas Series with index as taxon ids and values as lineages in format that is common to both of them
    @staticmethod
    def correlate_lineages(iPrimaryLineages,iSecondaryLineages):
        s_linages = iSecondaryLineages[iSecondaryLineages.isin(iPrimaryLineages)].drop_duplicates().reset_index()
        s_linages.columns = ['secondary_id','lineage']
        s_linages = s_linages.set_index('lineage')
        p_lineages = iPrimaryLineages.reset_index()
        p_lineages.columns = ['primary_id','lineage']
        p_lineages = p_lineages.set_index('lineage')
        correlated_lineages = p_lineages.join(s_linages)
        correlations = correlated_lineages.applymap(lambda x: x if pd.notna(x) else None).set_index('primary_id')
        return {'assigned':correlated_lineages.index,'unassigned':iPrimaryLineages[~iPrimaryLineages.isin(iSecondaryLineages)],'correlations':correlations.iloc[:,0]}
    
    #Cuts lineages by given levels starting from the begining or end
    @staticmethod
    def cut_lineages(iLineages,iLevels):#iLineages - pandas Series of lineages and ids as indices; iLevels - levels to cut; iReverse - by default methods starts cuting from begining if iReverse is set True count will start from end
        lineages = list(iLineages.values)
        lineage_indices = list(iLineages.index)
        iReverse = True if iLevels<0 else False
        iLevels = abs(iLevels)
        new_lineages = []
        for taxon_lineage in lineages:
            tmp_lineage = taxon_lineage
            for level in range(iLevels):
                tmp_lineage = tmp_lineage[(tmp_lineage.find(';')+2):] if (not iReverse) else tmp_lineage[:-(tmp_lineage[::-1].find(';')+1)]
            new_lineages.append(tmp_lineage)
        new_lineages_series = pd.Series(data=new_lineages,index=lineage_indices)  
        return new_lineages_series 
    
    #Reads CSV file and returns its content
    @staticmethod
    def _read_csv(iFilename,iSeparator=',', iQuote='"'): 
        file_content= []
        with open(iFilename, 'r') as otu_file:
            otu_file_reader = csv.reader(otu_file, delimiter=iSeparator, quotechar=iQuote)
            for row in otu_file_reader:
                file_content.append(row)
        return file_content
    
    #This method pickles and saves instance attributes
    def save_state(self,iFilename=False):
        c_dt = datetime.datetime.now()
        if iFilename:
            file_name = iFilename
        else:
            file_name = self.__class__.__name__+'_'+str(random.randint(1,1000))
        file_name = file_name + '('+c_dt.strftime('%Y-%m-%d(%H-%M-%S)')+')'
        with open(file_name,'wb') as state_file:
            pickle.dump(self.__dict__,state_file,3)
        return
    
    #This method unpickles and loads instance
    def load_state(self,iFilename):
        with open(iFilename,'rb') as state_file:
            tmp_dict = pickle.load(state_file)
        self.__dict__.update(tmp_dict)
        return

class JadvalOTU(JadvalMain):
    """JadvalOTU is responsible for loading, forming and maintaining, and manipulating OTU/Taxa data and sample metadata
    This class can load classic OTU tables, separate sample metadata or BIOM format files.
    """
    def __init__(self,iStudyInfo = None):
        self.study_info = iStudyInfo #Must be given as dictionary that contains information about study/research/data source. Such as author name, article title, publication year, date when data was collected, NGS device information, etc. This information will be used when combining jadval objects and printed into reports during report generation.
        self.otu_table = None #Pandas dataframe only with OTU reads. Dataframe has rows as OTU ids and columns as samples with sample ids
        self.taxon_groups = [] #List with dicts of Group ID-to-OTU ID associations and vice versa
        self.taxonomy = pd.DataFrame(columns=(['lineage']+JADVAL_MAIN_RANKS)) #Dataframe with OTU taxonomy with rows as OTU ids and columns with various formatted taxonomic information
        self.sample_metadata = None #Pandas dataframe for sample metadata with rows as samples ids and columns for datasets
        self.taxonomic_db = {} # Dictionary of loaded taxonomic databases
        self.label_assoc = {} #Dict with otu_id - label associations
        self.loaded_lineages = None #Pandas Series with original lineages
        self.lost_otus = [] #List of Dict with OTUs that was lost
        
    #Following method loads classic OTU table. Classic OTU table must have following structure. First row must contain headers. Columns: 1st for OTU ids; 2nd to last column for OTU reads, and if specified last column for taxonomy consensus lineages
    def load_OTU_table(self, iFilename, iSeparator=',', iQuote='"', iConsensusLineage=True, iReversed= False): #iFilename - must be entered; iSeparators - characters that separates data, by default is comma; iQuote - character used for quoting in csv file, by default character is <">; iConsensusLineage - if last column contains consensus lineages for every OTU, default is True; iReversed - if table is reversed, by default OTUs rows of first column and samples are amoung columns of first row
        pre_data  = self._read_csv(iFilename,iSeparator,iQuote)
        last_col = (len(pre_data[0])-1) if iConsensusLineage else len(pre_data[0])
        pre_sample_ids = pre_data[0][1:last_col]
        pre_otu_labels = [str(elem[0]) for elem in pre_data[1:]]
        pre_otu_ids = list(range(len(pre_otu_labels)))
        pre_otu_table = [elem[1:last_col] for elem in pre_data[1:]]
        pre_otu_lineages = [elem[last_col] for elem in pre_data[1:]] if iConsensusLineage else False
        self.label_assoc = {'id':pd.Series(data=pre_otu_labels,index=pre_otu_ids),'label':pd.Series(data=pre_otu_ids,index=pre_otu_labels)}
        self.otu_table = pd.DataFrame(data=pre_otu_table,index=pre_otu_ids,columns=pre_sample_ids)
        self.loaded_lineages  = pd.Series(data=pre_otu_lineages,index=pre_otu_labels)
        if iConsensusLineage:
            self.taxonomy = pd.DataFrame(index=pre_otu_ids,columns=self.taxonomy.columns)
            self.__construct_taxonomy_from_lineages(pd.Series(data=pre_otu_lineages,index=pre_otu_ids))
            self.reconstruct_internal_lineages()
        return
    
    #This function accepts single or many OTU labels and returns corresponding taxonomy as dictionaries.
    def get_taxonomy_by_otu_label(self,iOTU_labels):
        otu_ids = list(self.otu_label_assoc['label'].loc[iOTU_labels])
        return self.get_taxonomy_by_otu_id(otu_ids)
    
    #This function accepts single or many OTU ids and returns corresponding taxonomy as dictionaries.
    def get_taxonomy_by_otu_id(self,iOTU_ids):
        iOTU_ids = iOTU_ids if type(iOTU_ids)==list else [iOTU_ids] #IMPROVE!!! Check if iOTU_ids is a list and converts it into list if it is not
        all_taxa_dict = {} #Dictionary with multiple taxonomy
        for otu_id in iOTU_ids:
            taxa_dict = dict(self.taxonomy.loc[otu_id,JADVAL_MAIN_RANKS])
            #Following sets single letter rank dict keys to full ones that are defined in JADVAL_ITS
            taxa_dict = {JADVAL_ITS['r2Rank'][k]:v for k,v in taxa_dict.items()}
            all_taxa_dict[otu_id] = taxa_dict # Adds taxonomy of current OTU into dict with multiple taxonomies
        return all_taxa_dict if (len(iOTU_ids)>1) else taxa_dict #Returns multiple taxonomy dict if more than one OTU id was requested.
   
    #This function initiates self.taxonomy by breaking down lineages and storing ranks separately for further use
    def __construct_taxonomy_from_lineages(self,lineages):
        def allocater(lineage,ref_ranks):
            taxa_dict = {e[0]:e[1] for e in [e.strip().split('__') for e in lineage.split(';') if ('__' in e)]} #Efficiently converts lineage into dictionary
            taxa_dict_allowed = {rank:taxa_dict[rank] for rank in taxa_dict.keys() if rank in ref_ranks} #Drops unallowed ranks and orders ranks based on JADVAL_MAIN_RANKS rank order
            #Following loop sets unavailable ranks to None in order to avoid later problems
            for key in ref_ranks:
                if not (key in taxa_dict_allowed.keys()):
                    taxa_dict_allowed[key] = None
            taxa_list_ordered = [taxa_dict_allowed[rank] for rank in ref_ranks] #Orders ranks based on JADVAL_MAIN_RANKS rank order
            return [lineage]+taxa_list_ordered
        allocater_vectorized = np.vectorize(allocater,excluded=['ref_ranks'],otypes=[list])
        self.taxonomy = pd.DataFrame(index=list(lineages.index),data=list(allocater_vectorized(lineage=list(lineages.values),ref_ranks=JADVAL_MAIN_RANKS)),columns=self.taxonomy.columns)
        self.__fix_missing_taxa()
        return
    
    #This method sets taxons with "" = None. For example, lineages sometimes contain ranks such as (...; g__; s__) and (...; g__;). If consensus lineages are only available taxonomic information then such taxa are basically same and must be fixed.
    def __fix_missing_taxa(self):
        self.taxonomy.loc[:,JADVAL_MAIN_RANKS] = self.taxonomy.loc[:,JADVAL_MAIN_RANKS].applymap(lambda x: None if (x=='') else x)
        return
    
    #This method returns list of available rank levels. Method checks highest and lowest rank available in self.taxonomy that is also present among JADVAL_MAIN_RANKS
    def get_avail_ranks(self):
        return [rank for rank in JADVAL_MAIN_RANKS if self.taxonomy.loc[:,rank].notna().any()]
    
    #This method reconstructs self.taxonomy.loc[:,'lineage'] by initiated taxa
    def reconstruct_internal_lineages(self):
        self.taxonomy.loc[:,'lineage'] = self.generate_lineages(iMissingRank=True,iDesiredRanks=self.get_avail_ranks(),iDropRanks=False)
        return
    
    #This method records taxa that is to be removed into self.lost_otus
    def __record_removed_taxa(self,otu_ids,action_name):
        removed_otus = pd.concat([self.otu_table.loc[otu_ids,:],self.taxonomy.loc[otu_ids,:]],axis=1)
        removed_otus.index = self.label_assoc['id'][otu_ids]
        self.lost_otus.append({action_name:removed_otus})
        return removed_otus
    
    #This method removes given otu ids from self.label_assoc
    def __remove_labels_by_ids(self,ids):
        labels_to_remove = self.label_assoc['id'][ids] if type(self.label_assoc['id'][ids].iloc[0])!=list else list(itertools.chain.from_iterable(self.label_assoc['id'][ids]))
        self.label_assoc['id'] = self.label_assoc['id'].drop(ids)
        self.label_assoc['label'] = self.label_assoc['label'].drop(labels_to_remove)
        return labels_to_remove
    
    #This is main method that removes OTU from this instance by recording action into self.lost_otus; removing from self.label_assoc and droping otu from self.taxonomy and self.otu_table
    def __drop_taxa(self,ids,action_name):
        removed_otus = self.__record_removed_taxa(ids,action_name)
        self.__remove_labels_by_ids(ids)
        self.taxonomy = self.taxonomy.drop(ids)
        self.otu_table = self.otu_table.drop(ids)
        return removed_otus
        
        
    #This method removes OTUs based on ids
    def drop_otus_by_ids(self,iOTU_ids):
        iOTU_ids  = iOTU_ids if type(iOTU_ids)==list else [iOTU_ids]
        if self.label_assoc['label'].isin(iOTU_ids).sum() == len(iOTU_ids):
            return self.__drop_taxa(iOTU_ids,'manual-by-id')
        else:
            return False
        
    #This method removes OTUs based on labels
    def drop_otus_by_labels(self,iOTU_labels):
        iOTU_labels  = iOTU_labels if type(iOTU_labels)==list else [iOTU_labels]
        if self.label_assoc['id'].isin(iOTU_labels).sum() == len(iOTU_labels):
            iOTU_ids = list(self.label_assoc['label'][iOTU_labels].values)
            return self.__drop_taxa(iOTU_ids,'manual-by-label')
        else:
            return False

    #This method removes OTUs that do not have any assigned taxonomy from self.taxonomy dataframe and returns OTUs that was removed
    def drop_otus_without_taxa(self):
        otu_indices_drop = self.taxonomy.loc[self.taxonomy.loc[:,JADVAL_MAIN_RANKS].agg(lambda rank: len(''.join(map(lambda x: (str(x or '')),rank))),axis=1) < 1].index
        removed_otus = self.__drop_taxa(otu_indices_drop,'drop-otus-without-taxa')
        return removed_otus
    
    #This method drops otus that do not have taxon at given rank or ranks in iRanks. If iRanks contains more than one rank and iMerge = True then OTUs with any missing rank given in iRanks will be removed. Otherwise, if iMerge = False then for OTU to be removed, all ranks given in iRanks must be missing
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
    
    #This method reconstructs self.label_assoc after grouping
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
    
    #This function generates taxon_groups dict. taxon_groups['group'] - contains group ids and corresponding otu ids; taxon_groups['otus'] - contains otu ids and corresponding group ids
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
    
    #Merge OTU reads by duplicated lineages
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
    def merge_otus_by_rank(self,iRank):
        ranks = JADVAL_MAIN_RANKS[:JADVAL_MAIN_RANKS.index(iRank)+1]
        missing_ranks = len(JADVAL_MAIN_RANKS[JADVAL_MAIN_RANKS.index(iRank)+1:])
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
    def generate_lineages(self,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False 
        return self.generate_lineages_from_taxa(iTaxa=self.taxonomy,iMissingRank=iMissingRank,iDesiredRanks=iDesiredRanks,iDropRanks=iDropRanks)
    
    #This method loads JadvalDB and stores truncated database taxonomy and tree for further usage
    def load_database(self,iDerivationName,iJadvalDB_obj, iRankTolerance=False):
        correlation_result = self.__derive_taxonomy_db(iJadvalDB_obj, iRankTolerance)
        if not correlation_result:
            return False
        pruned_tree = iJadvalDB_obj.prune_tree_for_taxa(correlation_result['taxon_ids'])
        DBType = iJadvalDB_obj.__class__.__name__
        derivation_record = {'correlations':correlation_result['corr_as_series'],'tree':pruned_tree,'correlation_details':correlation_result}
        if DBType not in self.taxonomic_db.keys():
            self.taxonomic_db[DBType] = {iDerivationName:derivation_record}
        else:
            self.taxonomic_db[DBType][iDerivationName] = derivation_record
        return derivation_record
    
    
    
    #This method loads taxonomy database(JadvalDB) and correlates internal OTU ids with database taxon ids with rank tolerance given in iRank
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
        
    #This function correlate taxon ids between available lineages and iLineages. If iTolRanks is set then correlation of uncorrelated taxa will continue until all ranks are assigned or iTolRanks are exausted.
    def generate_taxa_correlations(self,iLineages,iForRanks,iTolRanks=False):
        self_lineages = self.generate_lineages(True,iForRanks)    
        correlation_results = self.correlate_lineages(self_lineages,iLineages)
        correlations = correlation_results['correlations']
        uncorrelated_lineages = correlation_results['unassigned']
        uncorrelated_otu_ids = uncorrelated_lineages.index
        correlated_lineages = correlation_results['assigned']
        if iTolRanks:
            rank_trials = {'init':{'corr':correlations.to_dict(),'uncorr_ids':uncorrelated_otu_ids,'added':sum(correlations.notna()),'lineages':correlation_results['assigned']}}
            for limiting_rank in [rank for rank in JADVAL_MAIN_RANKS[::-1] if rank in iTolRanks.keys()]:
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


    #This method checks if iJadvalOTU instance contains same taxonomy as self
    def is_jadval_taxa_alike(self,iJadvalOTU):
        self_lineage_sorted = self.taxonomy.loc[:,'lineage'].sort_values(axis=0).reset_index(drop=True)
        iJadvalOTU_lineage_sorted = iJadvalOTU.taxonomy.loc[:,'lineage'].sort_values(axis=0).reset_index(drop=True)    
        return self_lineage_sorted.equals(iJadvalOTU_lineage_sorted)
    
    
    ##### Properties 
    @property
    def jadval_(self):
        return pd.concat([self.otu_table,self.taxonomy],axis=1)
    
    @property
    def id2labels(self):
        return self.label_assoc['id']
    
    @property
    def label2id(self):
        return self.label_assoc['label']

    @property
    def lineages(self):
        return self.taxonomy.loc[:,'lineage']
    
    @property
    def original_lineages(self):
        return self.loaded_lineages
    

class JadvalTaxonomy(JadvalMain):
    """JadvalTaxonomy DGRPy class contain methods and variables shared amoung JadvalTaxonomyGreengenes, JadvalTaxonomySILVA, JadvalTaxonomyNCBI, JadvalTaxonomyOTL, JadvalTaxonomyRDP """
    def __init__(self):
        super().__init__()
        self.taxonomy_df = pd.DataFrame(columns=(['lineage']+JADVAL_MAIN_RANKS)) #Dataframe with OTU taxonomy with rows as OTU ids and columns with various formatted taxonomic information
        self.tree = None
        self.loaded_lineages = None #Pandas Series with original lineages

     #Following method loads greengeens taxonomic database
    def _load_tree(self,iTreeFile,iTreeFormat='newick'): #iTreeFile - Tree file in format given in iTreeFormat, by default is newick;
        with open(iTreeFile, 'r') as tree_file:
            content = tree_file.read()
        self.tree = ete3.Tree(content,quoted_node_names=True,format=1)
        return
    
    #This method prunes the database tree for given taxa and returns a copy
    def prune_tree_for_taxa(self,iTaxa_ids):
        if self.tree and (type(iTaxa_ids)==set):
            tmp_tree = copy.deepcopy(self.tree)
            tmp_tree.prune(iTaxa_ids, preserve_branch_length=True)
            return tmp_tree
        else:
            return False
    
    #This function generate desired lineages
    def generate_lineages(self,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False 
        return self.generate_lineages_from_taxa(iTaxa=self.taxonomy_df,iMissingRank=iMissingRank,iDesiredRanks=iDesiredRanks,iDropRanks=iDropRanks)
    
    #This method returns list of available rank levels. Method checks highest and lowest rank available in self.taxonomy_df that is also present among JADVAL_MAIN_RANKS
    def get_avail_ranks(self):
        return [rank for rank in JADVAL_MAIN_RANKS if self.taxonomy_df.loc[:,rank].notna().any()]
    
    #This method sets taxons with "" = None. For example, lineages sometimes contain ranks such as (...; g__; s__) and (...; g__;). If consensus lineages are only available taxonomic information then such taxa are basically same and must be fixed.
    def _fix_missing_taxa(self):
        self.taxonomy_df.loc[:,JADVAL_MAIN_RANKS] = self.taxonomy_df.loc[:,JADVAL_MAIN_RANKS].applymap(lambda x: None if (x=='') else x)
        return
    
    #This method reconstructs self.taxonomy_df.loc[:,'lineage'] by initiated taxa
    def reconstruct_internal_lineages(self):
        self.taxonomy_df.loc[:,'lineage'] = self.generate_lineages(iMissingRank=True,iDesiredRanks=self.get_avail_ranks(),iDropRanks=False)
        return 

class JadvalTaxonomyGreengenes(JadvalTaxonomy):
    """JadvalTaxonomyGreengenes DGRPy class that is responsible for loading, maintaining, and manipulating GreenGenes Taxonomy"""
    def __init__(self,iTaxaMapFile=None,iTreeFile=None,iTreeFormat='newick',iDBInfo = None,iStateFile=None):
        super().__init__()
        self.db_info = iDBInfo #Must be given as dictionary that contains information about study/research/data source. Such as author name, article title, publication year, date when data was collected, NGS device information, etc. This information will be used when combining jadval objects and printed into reports during report generation.
        if iStateFile:
            self.load_state(iStateFile)
        if iTreeFile and iTaxaMapFile:
            self.load_taxonomy_database(iTaxaMapFile,iTreeFile,iTreeFormat)
    
    #Following method loads taxonomic database
    def load_taxonomy_database(self,iTaxaMapFile,iTreeFile,iTreeFormat='newick'):       
        self.__load_taxonomy_map(iTaxaMapFile)
        self._load_tree(iTreeFile,iTreeFormat)
        self.__init_internal_taxonomy()
        self.reconstruct_internal_lineages()
        return
    
    #Following method loads greengeens taxonomic database
    def __load_taxonomy_map(self,iTaxaMapFile): #iDatabaseName - Name of database; iTreeFile - Tree file in format given in iTreeFormat, by default is newick; iTaxaMapFile - file with ID - Taxonomy data, by default values are separated by \t
        tmp_taxa_map = []
        with open(iTaxaMapFile, 'r') as map_file:
            taxa_map_file = csv.reader(map_file, delimiter='\t')
            for taxa in taxa_map_file:
                tmp_taxa_map.append(taxa)
        self.taxonomy_map = pd.Series(data=[e[1] for e in tmp_taxa_map],index=[e[0] for e in tmp_taxa_map])
        return
    
    #This function initiates self.taxonomy_df by breaking down lineages and storing ranks separately (Vectorization is not an option. Reoptimization is necessary!)
    def __init_internal_taxonomy(self):
        def allocater(lineage,ref_ranks):
            taxa_dict = {e[0]:e[1] for e in [e.strip().split('__') for e in lineage.split(';') if ('__' in e)]} #Efficiently converts lineage into dictionary
            taxa_dict_allowed = {rank:taxa_dict[rank] for rank in taxa_dict.keys() if rank in ref_ranks} #Drops unallowed ranks and orders ranks based on JADVAL_MAIN_RANKS rank order
            #Following loop sets unavailable ranks to None in order to avoid later problems
            for key in ref_ranks:
                if not (key in taxa_dict_allowed.keys()):
                    taxa_dict_allowed[key] = None
            taxa_list_ordered = [taxa_dict_allowed[rank] for rank in ref_ranks] #Orders ranks based on JADVAL_MAIN_RANKS rank order
            return [lineage]+taxa_list_ordered
        allocater_vectorized = np.vectorize(allocater,excluded=['ref_ranks'],otypes=[list])
        self.taxonomy_df = pd.DataFrame(index=list(self.taxonomy_map.index),data=list(allocater_vectorized(lineage=list(self.taxonomy_map.values),ref_ranks=JADVAL_MAIN_RANKS)),columns=self.taxonomy_df.columns)
        self._fix_missing_taxa()
        return True

class JadvalTaxonomySILVA(JadvalTaxonomy):
    """JadvalTaxonomySILVA DGRPy class that is responsible for loading, maintaining, and manipulating SILVA Taxonomy"""
    def __init__(self,iTaxaMapFile=None,iTreeFile=None,iTreeFormat='newick',iDBInfo = None,iStateFile=None):
        super().__init__()
        self.db_info = iDBInfo #Must be given as dictionary that contains information about study/research/data source. Such as author name, article title, publication year, date when data was collected, NGS device information, etc. This information will be used when combining jadval objects and printed into reports during report generation.
        if iStateFile:
            self.load_state(iStateFile)
        if iTreeFile and iTaxaMapFile:
            self.load_taxonomy_database(iTaxaMapFile, iTreeFile,iTreeFormat)
    
    #Following method loads taxonomic database
    def load_taxonomy_database(self, iTaxaMapFile, iTreeFile,iTreeFormat='newick'):       
        self.__load_taxonomy_map(iTaxaMapFile)
        self._load_tree(iTreeFile,iTreeFormat)        
        self.__init_internal_taxonomy()   
        self.reconstruct_internal_lineages()
        return
    
    #Following method loads greengeens taxonomic database
    def __load_taxonomy_map(self,iTaxaMapFile): #iTaxaMapFile - file with ID - Taxonomy data, by default values are separated by \t
        tmp_taxa_map = []
        with open(iTaxaMapFile, 'r') as map_file:
            taxa_map_file = csv.reader(map_file, delimiter='\t')
            for taxon in taxa_map_file:
                tmp_taxa_map.append(taxon)
        self.taxonomy_map = pd.DataFrame(data=[[e[0], e[2]] for e in tmp_taxa_map],index=[e[1] for e in tmp_taxa_map],columns=['lineage','level'])
        return
    
    #This function initiates self.taxonomy_df by breaking down lineages and storing ranks separately (This function is NOT OPTIMIZED. Cython can by considered)
    def __init_internal_taxonomy(self):
        SILVA_RANKS = JADVAL_MAIN_RANKS[:-1]
        tmp_taxonomy_df = pd.DataFrame(columns=(['lineage']+SILVA_RANKS),index=self.taxonomy_map.index,data=None)
        rank = SILVA_RANKS[0]
        tmp_taxonomy_df.loc[:,rank] = self.taxonomy_map.loc[self.taxonomy_map.loc[:,'level']==JADVAL_ITS['r2rank'][rank],:].apply(lambda row: row['lineage'][:-1].split(';')[0],axis=1)
        for rank in SILVA_RANKS[1:]:    
            tmp_taxonomy_df.loc[:,rank] = self.taxonomy_map.loc[self.taxonomy_map.loc[:,'level']==JADVAL_ITS['r2rank'][rank],:].apply(lambda row: row['lineage'][:-1].split(';'),axis=1)
        def reassign(taxons,r):
            new_taxa = {r:taxons.pop()}
            for r_i in SILVA_RANKS[:SILVA_RANKS.index(r)]:
                for t_i in range(len(taxons)):
                    if taxons[t_i] in tmp_taxonomy_df.loc[tmp_taxonomy_df.loc[:,r_i].notna(),r_i].values:
                        new_taxa[r_i] = taxons[t_i]
                        taxons.pop(t_i)
                        break;
            return pd.Series(new_taxa)
        for rank in SILVA_RANKS[1:]:    
            tmp_taxonomy_df.update(tmp_taxonomy_df.loc[tmp_taxonomy_df.loc[:,rank].notna(),rank].apply(reassign,r=rank))
        tmp_taxonomy_df.insert(loc=len(tmp_taxonomy_df.columns),column='s',value=None)
        tmp_taxonomy_df = tmp_taxonomy_df.applymap(lambda x: None if pd.isna(x) else x)
        self.taxonomy_df = tmp_taxonomy_df
        self._fix_missing_taxa()
        return
    
