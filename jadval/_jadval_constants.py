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


###Following constants should not be modified. JADVAL_ITS  links internal and external representation of values
ITS = {'lineage':'Consensus Lineage',
'r2Rank': {'d':'Domain','k':'Kingdom','p':'Phylum','c':'Class','o':'Order','f':'Family','g':'Genus','s':'Specie'},
'r2rank': {'d':'domain','k':'kingdom','p':'phylum','c':'class','o':'order','f':'family','g':'genus','s':'specie'},
'rank2r': {'domain':'d','kingdom':'k','phylum':'p','class':'c','order':'o','family':'f','genus':'g','specie':'s'},
'rank2Rank': {'domain':'Domain','kingdom':'Kingdom','phylum':'Phylum','class':'Class','order':'Order','family':'Family','genus':'Genus','specie':'Specie'},
'Rank2r': {'Domain':'d','Kingdom':'k','Phylum':'p','Class':'c','Order':'o','Family':'f','Genus':'g','Specie':'s'}}


#JADVAL_MAIN_RANKS list of main taxonomic ranks used in Jadval. Do not modify. Order is important!
MAIN_RANKS = ['d','k','p','c','o','f','g','s']

# ERROR CODES dict of all Jadval Error Code Descriptions
_ERROR_CODES = {
    1:'Error while generating lineages from taxonomy dataframe',#generate_lineages_from_taxa,
    2:'Error while performing primary lineage correlations',#correlate_lineages
    3:'Cutting Lineages Failed',#cut_lineages
    4:'Error While Reading CSV file',#read_csv
    5:'Saving Instance State Failed',#save_state
    6:'Loading Instance State Failed',#load_state
    7:'Loading of Tree File Failed',#load_tree
    8:'Tree Pruning Failed',#prune_tree_for_taxa_ids,
    9:'Can\'t get available ranks for initiated taxonomy',#get_avail_ranks
    10:'Error while fixing missing taxa', #_fix_missing_taxa
    11:'Failed to reconstruct internal lineages',#reconstruct_internal_lineages
    12:'Error during database loading',#load_taxonomy_database
    13:'Failed to initiate taxonomy map',#__load_taxonomy_map
    14:'Failed to initiate internal taxonomy dataframe',#__init_internal_taxonomy
    15:'Error Loading OTU Table',#load_OTU_table
    16:'Error! OTU Table must be initiated first',#load_OTU_lineages
    17:'Failed to Check State of Initiation',#check_init
    18:'Failed to load OTU lineages',#load_OTU_lineages
    19:'Failed to get taxonomy by OTU labels',#get_taxonomy_by_otu_label
    20:'Failed to get taxonomy by OTU ids',#get_taxonomy_by_otu_id
    21:'Failed to reconstruct internal lineages.',#reconstruct_internal_lineages
    22:'Failed to drop OTUs for given IDs.',#drop_otus_by_ids
    23:'Failed to drop OTUs for given labels.',#drop_otus_by_labels
    24:'Failed to drop OTUs without taxonomy.',#drop_otus_without_taxa
    25:'Failed to drop OTUs without given rank(s).',#drop_otus_without_ranks
    26:'Error while merging duplicated taxa.',#merge_duplicated_taxa
    27:'Error while merging OTUs by rank.',#merge_otus_by_rank
    28:'Error loading given taxonomic database.',#load_database
    29:'Failed to generate OTU-Taxon correlations.',#generate_taxa_correlations
    30:'Failed to reconstruct instance for new OTU IDs.',#reconstruct_for_otu_ids
    31:'Failed to retrieve custom OtuTable for required taxonomic database.',#retrieve_custom
    32:'Failed to get active taxonomic database of instance.',#get_active_taxonomic_db
    33:'Error while generating OTU table CSV file.',#save_as_csv
    34:'Error while exporting OtuTable data.',#export_data
    35:'Failed to set new internal IDS.',#__set_internal_ids
    36:'Failed to construct taxonomy from lineages',#__construct_taxonomy_from_lineages
    37:'Failed to derive new taxonomy for given taxonomic database.',#__derive_taxonomy_db
    38:'Failed to record <prior ID - new ID> group relations',#__record_group_relations
    39:'Failed to reassign <OTU ID - OTU Label> associations',#__reassign_label_assoc
    40:'Failed to record removed data to self.lost_otus',#__record_removed_taxa
    41:'Failed to remove OTU IDs and OTU labels self.label_assoc',#__remove_labels_by_ids
    42:'Failed to drop taxa by OTU IDs',#__drop_taxa
    44:'Failed to retrieve complete Jadval table',#jadval
    45:'Failed to retrieve lineages',#lineages
    46:'Failed to retrieve id-label Series',#LabelById
    47:'Failed to retrieve label-id Series',#IdByLabel
    48:'Failed to retrieve info of initial loadings',#loader_info
    49:'Failed to retrieve phylogenetic tree (ete3 Tree Object) ',#phylo_tree
    50:'Failed to retrieve active taxonomic database type',#taxonomy_db_type
    51:'Failed to retrieve active taxonomic database name',#taxonomy_db_name
    52:'Failed to create new directory.',#ensure_new_dir
    53:'Failed to make a copy.',#copy
    54:'Error occured during pattern search.',#search
    55:'Failed to get list of internal OTU ids.',#get_ids
    56:'Failed to get list of OTU labels.',#get_labels
}

_WARNING_CODES = {1:'No Active Taxonomic Database Found.'}
