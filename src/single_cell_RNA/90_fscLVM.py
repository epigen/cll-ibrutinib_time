#!/usr/bin/python
import os
import sys
import fscLVM
import pdb
import scipy as SP
from fscLVM import plotFactors, plotRelevance, plotLoadings, saveFA, dumpFA


# Data for development:
# dataFile = "/scratch/lab_bock/shared/projects/cll-time_course/results/single_cell_RNA/11_CellTypes/allDataBest_NoDownSampling_noIGH_Monos/fscLVM/matrix.csv.gz"
# outFolder = "/scratch/lab_bock/shared/projects/cll-time_course/results/single_cell_RNA/11_CellTypes/allDataBest_NoDownSampling_noIGH_Monos/fscLVM/fscLVM_files/"
# annoFile = "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt"
# nHidden = int("3")


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)


dataFile = sys.argv[1]
outFolder = sys.argv[2]
annoFile = sys.argv[3]
nHidden = int(sys.argv[4])

print dataFile
print outFolder
print annoFile
print type(nHidden)
print nHidden

if not os.path.exists(dataFile):
    raise Exception("fscLVM dataFile '" + dataFile + "' not found")


if not os.path.exists(outFolder):
    raise Exception("fscLVM outFolder '" + outFolder + "' not found")


if not os.path.exists(annoFile):
    raise Exception("fscLVM '" + annoFile + "' not found")


#Annotation file
# annoFile = os.path.join("/home/nfortelny/resources_labBock/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt") #MSigDB
# annoFile = os.path.join('/home/nfortelny/projects_nfortelny/fscLVM_test/data/h.all.v5.0.symbols.gmt') #MSigDB
annoDB   = 'MSigDB'
if not os.path.exists(annoFile):
    raise Exception("MSigDB annotation file needs to be downloaded manually")
#Note: the license of MSigDB does not permit redistribution of the raw annotation files.
#Please register at http://software.broadinstitute.org/gsea/msigdb and download the
#hallmark gene sets and place the file in data folder.
#annoFile = os.path.join(data_dir,'c2.cp.reactome.v4.0.symbols.gmt.txt') #REACTOME
#annoDB   = 'REACTOME'

# dataFile: csv file with log expresison values
# dataFile = os.path.join("data",'Buettneretal.csv.gz')
data = fscLVM.utils.load_txt(dataFile=dataFile,annoFiles=annoFile,annoDBs=annoDB)

#print statistics for the loaded dataset
print ("Loaded {:d} cells, {:d} genes".format(data['Y'].shape[0],data['Y'].shape[1]))
print ("Annotation: {:d} terms".format(len(data['terms'])))


#I: indicator matrix that assigns genes to pathways
I = data['I'] #if loaded from the hdf file change to I = data['IMSigDB']
#Y: log expresison values
Y = data['Y']
#terms: ther names of the terms
terms = data['terms']

#gene_ids: the ids of the genes in Y
gene_ids = data['genes']

#initialize FA instance, here using a Gaussian noise model and fitting 3 dense hidden factors
FA = fscLVM.initFA(Y, terms,I, gene_ids=gene_ids, noise='gauss', nHidden=3, minGenes=15)


#model training
FA.train()

#print diagnostics
FA.printDiagnostics()


np.savetxt(os.path.join(outFolder, "W.csv"), FA.getW(), delimiter=",")
np.savetxt(os.path.join(outFolder, "X.csv"), FA.getX(), delimiter=",")
np.savetxt(os.path.join(outFolder, "Z.csv"), FA.getZ(), delimiter=",")
np.savetxt(os.path.join(outFolder, "Zchanged.csv"), FA.getZchanged(), delimiter=",")
np.savetxt(os.path.join(outFolder, "Relevance.csv"), FA.getRelevance(), delimiter=",")
np.savetxt(os.path.join(outFolder, "Genes_Indices.csv"), FA.idx_genes, delimiter=",")
np.savetxt(os.path.join(outFolder, "Genes_IDs.csv"), FA.gene_ids, delimiter=",", fmt="%s")
np.savetxt(os.path.join(outFolder, "Terms.csv"), FA.getTerms(), delimiter=",", fmt="%s")
np.savetxt(os.path.join(outFolder, "Annotations.csv"), FA.getAnnotations(), delimiter=",")


# FA.idx_genes # 6635
# FA.gene_ids # 1490
# FA.getTerms() # 47
#
#
# FA.getW() # terms vs genes ? (1490 * 47)
# FA.getX() # Terms associated with cells (182 * 47)
# FA.getZ() # terms vs genes ? (1490 * 47)
# FA.getZchanged() # any genes gained / lost in assignement
# FA.getRelevance() # relevance value - number of terms (47)
# FA.getAnnotations() # TRUE/FALSE - terms vs genes(1490 * 47)
#
#
# ####
# #### For factors just plot FA.getX()
# ####
#
# ####
# ####    Get factor relevance
# ####
#
#
#
# ####
# ####    GET GENE LOADINGS
# ####
#
# term='G2m checkpoint'
# n_genes=20
# Zchanged = FA.getZchanged([term])[:,0]
# W        = FA.getW([term])[:,0]
# Z        = FA.getZ([term])[:,0]
# gene_labels = SP.array(FA.gene_ids)
#
# #plot weights
#
# Wabs = SP.absolute(W)*SP.absolute(Z) # this is the loading?
# gene_index = SP.argsort(-Wabs)[:n_genes]
#
# Igain = (Zchanged[gene_index]==1)
# Ielse = (Zchanged[gene_index]==0)
#
# #fig = plt.figure(figsize=(5,5))
# y = SP.arange(len(gene_index))
# Ielse.any()
# #    plt.plot(abs(W[gene_index][Ielse]*Z[gene_index][Ielse]),y[Ielse],'k.',label='pre annotated')
# Igain.any()
# #    plt.plot(abs(W[gene_index][Igain]*Z[gene_index][Igain]),y[Igain],'r.',label='gains')
#
# Wabs[gene_index]
