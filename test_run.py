#THIS IS EXAMPLE TEST RUN
#Note that loading of database files can take some time depending on your RAM and CPU.
#My computer with 16 GB RAM and Intel Core i7-4720HQ 2.60GHz finished in about 10 minutes
#Greengeens database files are downloaded from: ftp://ftp.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
#SILVA database files are downloaded from: https://www.arb-silva.de/no_cache/download/archive/release_132/Exports/taxonomy/
#!!!ENJOY!!!

import Jadval #Importing DGRPy's Jadval module

jadval_gg = Jadval.JadvalTaxonomyGreengenes('97_otu_taxonomy.txt','97_otus.tree') #Loading Greengeens database
jadval_silva = Jadval.JadvalTaxonomySILVA('tax_slv_ssu_132.txt','tax_slv_ssu_132.tre.txt') #Loading SILVA database

study1 = Jadval.JadvalOTU('Study1') #Creating JadvalOTU object
study1.load_OTU_table('OTUs_sample.csv') #Loading sample OTU table file with OTU ids at first column and consensus lineages at last column
study1.drop_otus_without_taxa() #Removing OTUs without taxonomy
study1.merge_duplicated_taxa() #Merging OTUs with same taxonomy
study1.drop_otus_without_ranks('g') #Removing OTUs that does not have missing genus rank
study1.merge_otus_by_rank('g') #Merging OTUs by with same taxonomy up to genus rank (['d','k','p','c','o','f','g'])
g_gg =  study1.load_database('genus_gg',jadval_gg) #Loading greengeens database files, correlating OTU taxon ids with database taxon ids based on shared lineages, pruning tree for shared taxons and saving correlations and pruned tree into study1.taxonomic_db
g_silva =  study1.load_database('genus_silva',jadval_silva) #Same but for SILVA database

jadval_gg.save_state() #Pickling and saving JadvalTaxonomyGreengenes instance
jadval_silva.save_state() #Pickling and saving JadvalTaxonomySILVA instance
study1.save_state() #Pickling and saving JadvalOTU instance

g_gg['tree'].render('g_gg.pdf') #Rendering pruned greengeens tree and saving as g_gg.pdf
g_silva['tree'].render('g_silva.pdf')#Rendering pruned SILVA tree and saving as g_gg.pdf