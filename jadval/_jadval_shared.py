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


import csv
import pandas as pd
import random
import datetime
import pickle
import os

from ._jadval_error import error_handler as jErrorHandler
import jadval._jadval_constants as jConst


def ensure_list(var):
    return var if type(var) is list else [var]

def assert_taxonomic_db(var):
    if type(var) is dict:
        return all([key in var.keys() for key in ['name','type','correlations','tree','details']])
    return False

#This method checks if iOtuTable2 instance contains same taxonomy as iOtuTable1
def is_jadval_taxa_alike(iOtuTable1,iOtuTable2):
    iOtuTable1_lineage_sorted = iOtuTable1.taxonomy.loc[:,'lineage'].sort_values(axis=0).reset_index(drop=True)
    iOtuTable2_lineage_sorted = iOtuTable2.taxonomy.loc[:,'lineage'].sort_values(axis=0).reset_index(drop=True)
    return iOtuTable1_lineage_sorted.equals(iOtuTable2_lineage_sorted)


class JadvalPrime:
    """Jadval is main class of DGRPy and parental class for all others. This class contains definitions, variables, constants and methods shared with other classes."""

    @staticmethod
    @jErrorHandler(1)
    def generate_lineages_from_taxa(iTaxa,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False    
        if iDesiredRanks and not all(e in jConst.MAIN_RANKS for e in iDesiredRanks):
            print('Impossible characters found in iDesiredRanks. Please use: '+','.join(jConst.MAIN_RANKS))
            return False
        iDropRanks = iDropRanks if iDropRanks else []
        add_next = (lambda r,t: ((r+'__'+t+'; ') if t != None else r+'__'+'; ')) if iMissingRank else (lambda r,t: ((r+'__'+t+'; ') if t != None else '')) 
        def make_lineage(taxon):
            tmp_lineage = ''
            for rank in make_ranks:
                tmp_lineage = tmp_lineage + add_next(rank,taxon.loc[rank] if (rank not in iDropRanks) else None)   #Add succeeding rank to lineage
            return tmp_lineage[:-2]
        if not iDesiredRanks:
            make_ranks = jConst.MAIN_RANKS
            new_lineages = iTaxa.apply(make_lineage,axis=1)
        else:
            make_ranks = [rank for rank in jConst.MAIN_RANKS if rank in iDesiredRanks]
            new_lineages = iTaxa.apply(make_lineage,axis=1)
        return new_lineages
    
    #This function find common lineages between iPrimaryLineages and iSecondaryLineages.
    #Result contains dict with {'assigned'} - reserved for lineages that are common to both lineage series
    #{'unassigned'} - reserved for lineages that are not shared between two lineage series
    #{'correlations'} - contains correlation ids between two lineages
    #Both iPrimaryLineages and iSecondaryLineages must be pandas Series with index as taxon ids and values as lineages in format that is common to both of them
    @staticmethod
    @jErrorHandler(2)
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
    @jErrorHandler(3)
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
    @jErrorHandler(4)
    def read_csv(iFilename,iSeparator=',', iQuote='"'):
        file_content= []
        with open(iFilename, 'r') as otu_file:
            otu_file_reader = csv.reader(otu_file, delimiter=iSeparator, quotechar=iQuote)
            for row in otu_file_reader:
                file_content.append(row)
        return file_content

    #Write CSV file and returns its content
    @staticmethod
    @jErrorHandler(17)
    def write_csv(iContent,iFilename,iSeparator=',', iQuote='"'):
        with open(iFilename, 'w') as write_file:
            file_writer = csv.writer(write_file, delimiter=iSeparator, quotechar=iQuote)
            return file_writer.writerows(iContent)
        return False

    @jErrorHandler(5)
    def save_state(self,iFilename=False):
        c_dt = datetime.datetime.now()
        if iFilename:
            file_name = iFilename
        else:
            file_name = self.__class__.__name__+'_'+str(random.randint(1,1000)) + '(' + c_dt.strftime('%Y-%m-%d(%H-%M-%S)') + ')'
        with open(file_name,'wb') as state_file:
            pickle.dump(self.__dict__,state_file,3)
        return

    @jErrorHandler(6)
    def load_state(self,iFilename):
        with open(iFilename,'rb') as state_file:
            tmp_dict = pickle.load(state_file)
        self.__dict__.update(tmp_dict)
        return

    @staticmethod
    @jErrorHandler(52)
    def ensure_new_dir(dir_name):
        new_path_result = False
        if os.path.exists(dir_name):
            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)
                new_path_result = dir_name
            else:
                if '-' in dir_name:
                    sep_i = len(dir_name) - dir_name[::-1].index('-') - 1
                    no_str = dir_name[sep_i + 1:]
                    try:
                        cur_no = int(no_str)
                    except ValueError:
                        cur_no = None
                else:
                    cur_no = None
                dir_name = dir_name+'-' if cur_no is None else dir_name[:sep_i+1]
                n_i = 1 if cur_no is None else cur_no
                while os.path.exists(dir_name+str(n_i)):
                    n_i += 1
                os.mkdir(dir_name+str(n_i))
                new_path_result = dir_name+str(n_i)
        else:
            os.mkdir(dir_name)
            new_path_result = dir_name
        return new_path_result


