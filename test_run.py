#THIS IS EXAMPLE TEST RUN
#Note that loading of database files can take some time depending on your RAM and CPU.
#My computer with 16 GB RAM and Intel Core i7-4720HQ 2.60GHz finished in about 30 minutes
#Greengeens database files are downloaded from: ftp://ftp.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
#SILVA database files are downloaded from: https://www.arb-silva.de/no_cache/download/archive/release_132/Exports/taxonomy/
#!!!ENJOY!!!

import jadval

#Loading taxonomic databases
db_gg = jadval.jGreengenes(iTaxaMapFile='97_otu_taxonomy.txt',iTreeFile='97_otus.tree')
db_silva = jadval.jSILVA(iTaxaMapFile='tax_slv_ssu_132.txt',iTreeFile='tax_slv_ssu_132.tre.txt')

#Initiating new Study and loading OTU data
sample = jadval.jOtuTable('Sample')
sample.load_OTU_table('OTUs_sample.csv')
#Filtering OTU data
sample.drop_otus_without_taxa()
sample.merge_duplicated_taxa()
Wolbachia_OTUs = sample.search('Wolbachia')
sample.drop_otus_by_ids(Wolbachia_OTUs)

#Making separate OtuTable objects for separate ranks
sample_genus = sample
sample_family = sample.copy()
sample_order = sample.copy()
sample_genus_TolOrder = sample.copy()

#Preparing OTU data for desired ranks
sample_genus.drop_otus_without_ranks('g')
sample_genus.merge_otus_by_rank('g')

sample_family.drop_otus_without_ranks('f')
sample_family.merge_otus_by_rank('f')

sample_order.drop_otus_without_ranks('o')
sample_order.merge_otus_by_rank('o')

sample_genus_TolOrder.drop_otus_without_ranks('g')
sample_genus_TolOrder.merge_otus_by_rank('g')

#Loading databases for each rank and saving states
sample_genus_gg = sample_genus.load_database('sample_genus_gg',db_gg)
sample_genus_silva = sample_genus.load_database('sample_genus_silva',db_silva)
sample_genus.save_state('sample_genus.state')

sample_family_gg = sample_family.load_database('sample_family_gg',db_gg)
sample_family_silva = sample_family.load_database('sample_family_silva',db_silva)
sample_family.save_state('sample_family.state')

sample_order_gg = sample_order.load_database('sample_order_gg',db_gg)
sample_order_silva = sample_order.load_database('sample_order_silva',db_silva)
sample_order.save_state('sample_order.state')

sample_genus_TolOrder_gg = sample_genus_TolOrder.load_database('sample_genus_TolOrder_gg',db_gg)
sample_genus_TolOrder_silva = sample_genus_TolOrder.load_database('sample_genus_TolOrder_silva',db_silva)
sample_genus_TolOrder.save_state('sample_genus_TolOrder.state')

#Retrieving custom OtuTables for each database
sample_genus_gg = sample_genus.retrieve_custom(sample_genus_gg)
sample_genus_gg.study_name = 'sample_genus_gg'
sample_genus_silva = sample_genus.retrieve_custom(sample_genus_silva)
sample_genus_silva.study_name = 'sample_genus_silva'

sample_family_gg = sample_family.retrieve_custom(sample_family_gg)
sample_family_gg.study_name = 'sample_family_gg'
sample_family_silva = sample_family.retrieve_custom(sample_family_silva)
sample_family_silva.study_name = 'sample_family_silva'

sample_order_gg = sample_order.retrieve_custom(sample_order_gg)
sample_order_gg.study_name = 'sample_order_gg'
sample_order_silva = sample_order.retrieve_custom(sample_order_silva)
sample_order_silva.study_name = 'sample_order_silva'

sample_genus_TolOrder_gg = sample_genus_TolOrder.retrieve_custom(sample_genus_TolOrder_gg)
sample_genus_TolOrder_gg.study_name = 'sample_genus_TolOrder_gg'
sample_genus_TolOrder_silva = sample_genus_TolOrder.retrieve_custom(sample_genus_TolOrder_silva)
sample_genus_TolOrder_silva.study_name = 'sample_genus_TolOrder_silva'

#Exporting data
sample_genus_gg.export_data()
sample_genus_silva.export_data()
sample_family_gg.export_data()
sample_family_silva.export_data()
sample_order_gg.export_data()
sample_order_silva.export_data()
sample_genus_TolOrder_gg.export_data()
sample_genus_TolOrder_silva.export_data()
