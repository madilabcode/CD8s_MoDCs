# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 11:58:47 2024

@author: Kfir Inbal aka. Lionel Christopher Bell
"""
#%%
import os
path = os.path.dirname(os.path.abspath(__file__))
os.chdir(path)

#%%
#Loading libraries
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from biomart import BiomartServer
import statistics
import math
import matplotlib.pyplot as plt
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
from operator import getitem
from sklearn.preprocessing import MinMaxScaler
import networkx as nx

#import pyreadr

#%%

#Function that returns the correlation value of 2 genes from the correlation matrix
def get_correlation_score(row):
    input_gene =  row['Input-node Gene Symbol']
    output_gene = row['Output-node Gene Symbol']
    
    if input_gene in corr_matrix.index and output_gene in corr_matrix.index:
        return corr_matrix.loc[input_gene, output_gene]
    else: #if the pair is not in the correlation matrix
        return pd.NA

#Function that converts human gene symbols to mouse gene symbols. 
#It returns a dataframe where the first column is Human gene symbols and the second column is the corresponding mouse gene symbols.
#If there are more than 1 corresponding mouse gene symbols, then it prioritses the mouse gene symbol that is the exact name of Human gene symbol but with lowercase letters.
#If there is no corresponding mouse gene symbols then it returns None to that human symbol
def human_to_mouse_gene_symbols(human_gene_symbols_df, batch_size=100):
    # Connect to BioMart database
    server = BiomartServer("http://www.ensembl.org/biomart")
    dataset = server.datasets['hsapiens_gene_ensembl']
    
    # Initialize an empty DataFrame to store the results
    result_df = pd.DataFrame(columns=['Human Gene Symbol', 'Mouse Gene Symbol'])
    
    # Iterate over batches of human gene symbols
    for i in range(0, len(human_gene_symbols_df), batch_size):
        # Get a batch of human gene symbols
        batch_df = human_gene_symbols_df.iloc[i:i+batch_size]
        batch_symbols = batch_df.iloc[:, 0].tolist()
        
        # Query BioMart for mouse orthologs of the human gene symbols
        response = dataset.search({
               'filters': {
                   'external_gene_name': batch_symbols
               },
               'attributes': ['external_gene_name', 'mmusculus_homolog_associated_gene_name'],
               'species': 'mmusculus'
        })
        
        # Create a dictionary to store mouse gene symbols for each human gene symbol
        human_to_mouse_mapping = {}
        
        for record in response.iter_lines():
            decoded_record = record.decode('utf-8')
            if decoded_record.strip():  # Skip empty lines
                fields = decoded_record.split('\t')
                human_gene_symbol = fields[0]
                mouse_gene_symbols = fields[1].split(',')
                if len(mouse_gene_symbols) >= 2:
                # Filter mouse gene symbols to prioritize lowercase matches
                    mouse_gene_symbol = next((symbol for symbol in mouse_gene_symbols if symbol.lower() == human_gene_symbol.lower()), None)
                else:
                    mouse_gene_symbol = mouse_gene_symbols[0]
                if mouse_gene_symbol:
                    human_to_mouse_mapping[human_gene_symbol] = mouse_gene_symbol
                    
                    #data.append([human_gene_symbol, mouse_gene_symbol])
        print(human_to_mouse_mapping)
        
        
        # Populate the mouse gene symbols for the current batch
        #for index, row in batch_df.iterrows():
        for human_gene_symbol in list(batch_df.iloc[:, 0]):
            mouse_gene_symbol = human_to_mouse_mapping.get(human_gene_symbol)
            result_df = pd.concat([result_df, pd.DataFrame({'Human Gene Symbol': [human_gene_symbol], 'Mouse Gene Symbol': [mouse_gene_symbol]})], axis=0).reset_index(drop=True)
        print(result_df)
    return result_df


#%%
'''#Converting Human TF table to Mouse TF table'''
'''#LONG RUN TIME#'''
#Reading Human TFs table from Ron and dropping irrelevant columns
TF_Table = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\doro_net.csv").drop(['Unnamed: 0','confidence','mor'],axis = 1)
TF_Table_Mouse = TF_Table.copy()
#Converting the source human gene symbols to corresponding mouse gene symbols
TF_Table_Mouse['Mouse_Source'] = human_to_mouse_gene_symbols(pd.DataFrame(TF_Table['source']))['Mouse Gene Symbol']
#Converting the target human gene symbols to corresponding mouse gene symbols
TF_Table_Mouse['Mouse_Target'] = human_to_mouse_gene_symbols(pd.DataFrame(TF_Table['target']))['Mouse Gene Symbol']
TF_Table_Mouse = TF_Table_Mouse.drop(['source','target'],axis = 1)
#removing Nones
TF_Table_Mouse_noNones = TF_Table_Mouse.dropna()
TF_Table_Mouse = TF_Table_Mouse_noNones.copy()
#Saving in csv file
#TF_Table_Mouse.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Table_Mouse.csv")



#%%

List_Of_Samples = ["T_PCD11_LMoDC_1","T_PCD11_LMoDC_2","T_PCD11_LMoDC_3","T_PCD11_LMoDC_4"]
print("Samples:", List_Of_Samples)
#Loading tables
#########
########Bulk RNA-seq table Raw Data (To check enrichment of TFs):
T_Cells_Priming_Licensing_Raw_Counts = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\T_Cells\Priming_Licensing_Only\T_Cells_Treatments_Only_Raw_Data.csv", index_col = 0)
#Filtering to only T_PcDC_LMoDC samples
T_Cells_Priming_Licensing_Raw_Counts_PcDC_LMoDC = T_Cells_Priming_Licensing_Raw_Counts[List_Of_Samples]
#TRANS
T_Cells_Priming_Licensing_Raw_Counts_PcDC_LMoDC_trans = T_Cells_Priming_Licensing_Raw_Counts_PcDC_LMoDC.T

#########

#directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\directed_mouse_ppi.csv",index_col = 0)
#directed_mouse_PPI = directed_mouse_PPI.drop(['Input-node GeneID','Output-node GeneID','Edge direction score'],axis = 1)
#New PPI (mippie V1_0)
#directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\mippie_ppi_v1_0_symbols.csv",index_col = 0)
directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\mippie_ppi_v1_0_symbols_NEW_getBM.csv", index_col = 0)
directed_mouse_PPI.reset_index(drop=True, inplace=True)
#Try directed_mouse_PPI_2 with reset_index and see what it gives.
#directed_mouse_PPI_2.reset_index(drop=True, inplace=True)
#print((directed_mouse_PPI['GeneNameA'] == directed_mouse_PPI_2['GeneNameA']).all())
#print(directed_mouse_PPI['GeneNameA'].isin(directed_mouse_PPI_2['GeneNameA']).all())
#directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\mippie_ppi_v1_0_symbols_NEW.csv",index_col = 0)
#directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\pickle_ppi_Standard_Gene_Normalized_v3_3_symbols.csv",index_col = 0)
#directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\pickle_ppi_Standard_Gene_Normalized_v3_3_FULL_symbols.csv",index_col = 0)
#directed_mouse_PPI = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\pickle_ppi_Standard_Gene_Normalized_v3_3_FULL_symbols_No_Homology_Conversion.csv",index_col = 0)
directed_mouse_PPI.rename(columns = {"GeneNameA" : 'Input-node Gene Symbol'}, inplace = True)
directed_mouse_PPI.rename(columns = {"GeneNameB" : 'Output-node Gene Symbol'}, inplace = True)



########Bulk RNA-seq table Normalized Data to multiply by capacity later on:
T_Cells_All_Norm_Counts = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\T_Cells\T_Cells_Normalized_Data.csv", index_col = 0)

#Apply MinMax Normalization:
scaler = MinMaxScaler(feature_range=(0,max(list(T_Cells_All_Norm_Counts.max()))))
#scaler = MinMaxScaler(feature_range=(0,1))
T_Cells_Norm_Counts_MinMax = pd.DataFrame(scaler.fit_transform(T_Cells_All_Norm_Counts),columns = T_Cells_All_Norm_Counts.columns, index = T_Cells_All_Norm_Counts.index)
#Filtering to only T_PcDC_LMoDC samples
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax = T_Cells_Norm_Counts_MinMax[List_Of_Samples]
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax["MeanExp"] = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.mean(axis=1)
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax[["MeanExp"]]

'''
#Filtering to only T_PcDC_LMoDC samples
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC = T_Cells_All_Norm_Counts[List_Of_Samples]
#Apply MinMax Normalization:
scaler = MinMaxScaler(feature_range=(0,max(list(T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC.max()))))
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax = pd.DataFrame(scaler.fit_transform(T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC),columns = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC.columns, index = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC.index)
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax["MeanExp"] = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.mean(axis=1)
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax[["MeanExp"]]
'''
########Bulk RNA-seq table Normalized Data (To create correlation matrix):
#All T cells samples
T_Cells_All_Norm_Counts_trans = T_Cells_All_Norm_Counts.T
#Pearson correlation matrix:
corr_matrix = T_Cells_All_Norm_Counts_trans.corr(method = "pearson") #Pearson is default


plt.hist(abs(corr_matrix.loc["Csf2rb"]))
plt.show()

#1) different correlation method, look at Csf2rb.
#2) add all samples of T cells.


#Renaming Kars1 to its synonym Kars because it skips a receptor because of it, the correlation value exists for Kars name and not Kars1 name
directed_mouse_PPI.replace("Kars1", "Kars", inplace=True)



#Adding a column of correlation_score for each pair 
directed_mouse_PPI['correlation_score'] = directed_mouse_PPI.apply(get_correlation_score, axis=1)
nan_values = directed_mouse_PPI[directed_mouse_PPI.isnull().any(axis=1)]
#Removing NAs
directed_mouse_PPI = directed_mouse_PPI.dropna()
#Removing rows with correlation score of 1 #NOT A MUST
directed_mouse_PPI = directed_mouse_PPI[directed_mouse_PPI.correlation_score != 1]
#Resetting indexes
directed_mouse_PPI = directed_mouse_PPI.reset_index(drop=True)
plt.hist(abs(directed_mouse_PPI["correlation_score"]))
plt.show()

#Icam1 to Icam1_R
directed_mouse_PPI.replace("Icam1", "Icam1_R", inplace=True)
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.rename(index={'Icam1': 'Icam1_R'})

#Atp5f1b to Atp5b
directed_mouse_PPI.replace("Atp5f1b", "Atp5b", inplace=True)
T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax = T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.rename(index={'Atp5f1b': 'Atp5b'})



#%%

#ENRICHMNET TEST HYPERGEOMETRIC DISTRIBUTION

rna_seq_df = T_Cells_Priming_Licensing_Raw_Counts_PcDC_LMoDC_trans.copy()
#Reading Mouse TFs table
TF_Table_Mouse = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Table_Mouse.csv", index_col = 0)

# Dictionary to store enrichment results for each TF
enrichment_results = {}

# Total number of genes in the expression table
total_genes = len(rna_seq_df.columns)

# Total number of differentially expressed genes across all samples

#DEG = for raw counts, expression > 10 in all samples
#.all() makes sure its greater than 10 in all samples(rows).
#sum() sums how many genes have expression greater than 10 in all samples.
total_diff_expr_genes = sum((rna_seq_df > 10).all())

for tf in TF_Table_Mouse['Mouse_Source'].unique():
    # Extract target genes for the current TF
    target_genes = TF_Table_Mouse[TF_Table_Mouse['Mouse_Source'] == tf]['Mouse_Target'].tolist()
    
    # Create a new list of column names that exist in the dataframe
    valid_columns = [col for col in target_genes if col in rna_seq_df.columns]
    
    # Count the number of differentially expressed target genes
    #.all() makes sure its greater than 10 in all samples(rows).
    #sum() sums how many genes have expression greater than 10 in all samples.
    diff_expr_target_genes = sum((rna_seq_df[valid_columns] > 10).all())
  
    # Number of target genes for the current TF
    num_target_genes = len(target_genes)
    
    # Perform hypergeometric test
        #M = The population size, total number of genes in the bulkRNAseq expression table.
        #n = The number of successes in the population, total number of genes that are differentially expressed in the table.
        #N = The sample size. Number of target genes we check for the particular TF.
        #k = The number of drawn successes, target genes that are differentially expressed.
    p_value = hypergeom.sf(M=total_genes, n=total_diff_expr_genes, N=num_target_genes,k=diff_expr_target_genes - 1)
    
    
    enrichment_results[tf] =   {
         'Target Gene': target_genes,
         'Number of Differentially Expressed Target Genes': diff_expr_target_genes,
         'Total Number of Target Genes': num_target_genes,
         'Hypergeometric p-value': p_value
     }

    
# Sort keys based on the pvalues of the nested dictionaries
sorted_enrichment_results = dict(sorted(enrichment_results.items(), key=lambda x: x[1]['Hypergeometric p-value']))

# Correct for multiple testing
corrected_p_values = multipletests([value['Hypergeometric p-value'] for value in sorted_enrichment_results.values()], method='bonferroni')[1]

#Adding the adjusted pvalue for each TF in the nested dict
i = 0
for key in sorted_enrichment_results.keys():
    if i == len(corrected_p_values):
        break
    else:
        sorted_enrichment_results[key]['Adjusted p-value'] = corrected_p_values[i]
        i+=1

#Creating Dataframe that would store the results from the nested dictionary
TF_DF_results = pd.DataFrame(columns = ['TF','Total_Num_Targets','Num_Targets_DE','p-value','Adj. p-value'])

#Filling the dataframe
for key in sorted_enrichment_results.keys():
    TF_DF_results.loc[len(TF_DF_results)] = [key,sorted_enrichment_results[key]['Total Number of Target Genes'], sorted_enrichment_results[key]['Number of Differentially Expressed Target Genes'] ,sorted_enrichment_results[key]['Hypergeometric p-value'],
                                           sorted_enrichment_results[key]['Adjusted p-value']]

#Saving to csv file
TF_DF_results.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_TFs_In_T_Cell_PcDC_LMoDC_Results.csv")

#Filtering by pvalue < 0.05
filt_pval_TF_DF_results005 = TF_DF_results[TF_DF_results['p-value'] < 0.05]
#Filtering by adjusted p val < 0.05
filt_FDR_TF_DF_results005 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.05]
#Filtering by adjusted p val < 0.01
filt_FDR_TF_DF_results001 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.01]

#Saving to csv files
filt_pval_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_TFs_In_T_Cell_PcDC_LMoDC_Results_Pval_lt_0.05.csv")
filt_FDR_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_TFs_In_T_Cell_PcDC_LMoDC_Results_FDR_lt_0.05.csv")
filt_FDR_TF_DF_results001.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_TFs_In_T_Cell_PcDC_LMoDC_Results_FDR_lt_0.01.csv")

#%%

#TF lists
list_TFs_FDR005 = filt_FDR_TF_DF_results005['TF'].tolist()
list_TFs_FDR005_Exist_In_PPI = []

#Checking if the TF is either in Input-node or in Output-node then add to list list_TFs_FDR005_Exist_In_PPI; If in neither of them then I give up on that TF
for TF in list_TFs_FDR005:
    if (TF in directed_mouse_PPI['Input-node Gene Symbol'].tolist()) or (TF in directed_mouse_PPI['Output-node Gene Symbol'].tolist()):
        list_TFs_FDR005_Exist_In_PPI.append(TF)

directed_mouse_PPI_Virtual_Sink_abs = directed_mouse_PPI.copy()
directed_mouse_PPI_Virtual_Sink_MinMax = directed_mouse_PPI.copy()


#MinMax scaling with range 0 to 1 to change negative values to positive ones.
scaler = MinMaxScaler(feature_range=(0, 1))
directed_mouse_PPI_Virtual_Sink_MinMax[['correlation_score']] = scaler.fit_transform(directed_mouse_PPI_Virtual_Sink_MinMax[['correlation_score']])

#Alternatively abs to change negative values to positive ones.
directed_mouse_PPI_Virtual_Sink_abs['correlation_score'] = abs(directed_mouse_PPI_Virtual_Sink_abs['correlation_score'])



#Adding rows to directed_mouse_PPI so that each enriched TF would be connected to virtual sink with score of infinity, float('inf')
for TF in list_TFs_FDR005_Exist_In_PPI:
    directed_mouse_PPI_Virtual_Sink_MinMax.loc[len(directed_mouse_PPI_Virtual_Sink_MinMax)] = [TF,'V_S',math.inf]


directed_mouse_PPI_Virtual_Sink_MinMax = directed_mouse_PPI_Virtual_Sink_MinMax.rename(columns={"Input-node Gene Symbol": "Source", "Output-node Gene Symbol": "Target", "correlation_score":"capacity"})

#Saving to csv file
directed_mouse_PPI_Virtual_Sink_MinMax.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\directed_mouse_PPI_MinMax_Weighted_Flow_Network_EnrichedTFsInf_FDR005_In_T_Cell_PcDC_LMoDC.csv")


#Adding rows to directed_mouse_PPI so that each enriched TF would be connected to virtual sink with score of infinity, float('inf')
for TF in list_TFs_FDR005_Exist_In_PPI:
    directed_mouse_PPI_Virtual_Sink_abs.loc[len(directed_mouse_PPI_Virtual_Sink_abs)] = [TF,'V_S',math.inf]


directed_mouse_PPI_Virtual_Sink_abs = directed_mouse_PPI_Virtual_Sink_abs.rename(columns={"Input-node Gene Symbol": "Source", "Output-node Gene Symbol": "Target", "correlation_score":"capacity"})

#Saving to csv file
directed_mouse_PPI_Virtual_Sink_abs.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\directed_mouse_PPI_abs_Weighted_Flow_Network_EnrichedTFsInf_FDR005_In_T_Cell_PcDC_LMoDC.csv")

#%%
#Multiply capacity by meanExp of LMoDC

for s,t in zip(directed_mouse_PPI_Virtual_Sink_abs.Source, directed_mouse_PPI_Virtual_Sink_abs.Target):
    #if s in T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.index:
    s_Exp = float(T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.loc[s])
    
    #if t in T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.index:
    if t != "V_S":
        t_Exp = float(T_Cells_Priming_Licensing_Norm_Counts_PcDC_LMoDC_MinMax.loc[t])
    else:
        t_Exp = math.inf
        
    Old_Cap = float(directed_mouse_PPI_Virtual_Sink_abs.query('Source==@s and Target==@t')['capacity'])
 
    multiply  = (statistics.mean([s_Exp,t_Exp])) #min([s_Exp,t_Exp])#  #max([s_Exp,t_Exp]) #
    #print(multiply)
    new_capacity = Old_Cap * multiply
    row_index = directed_mouse_PPI_Virtual_Sink_abs.query('Source==@s and Target==@t')['capacity'].index[0]
    directed_mouse_PPI_Virtual_Sink_abs.at[row_index, 'capacity'] = new_capacity
    #print("For Source",s,"and Target",t,"the previous capicity was:",Old_Cap,"; And the new is:",new_capacity)

#%%

plt.hist((directed_mouse_PPI_Virtual_Sink_abs[directed_mouse_PPI_Virtual_Sink_abs["capacity"] != math.inf]["capacity"]))
plt.show()

#def returnMaxFlowPath(receptor, target, flow_dict):
    #for key in flow_dic
    #

def listOfTuples(l1, l2):
    return list(map(lambda x, y:(x,y), l1, l2))

def extract_flow_path(flow_dict, start, end):
    path = [start]
    flow_path_ = []
    current_node = start
    while current_node != end:
        next_nodes = flow_dict[current_node].copy()
        next_node = max(next_nodes, key=next_nodes.get)
        flow_= next_nodes[next_node]
        while(next_node in path and next_nodes):
            next_nodes.pop(next_node)
            if next_nodes:
                next_node = max(next_nodes, key=next_nodes.get)
                flow_ = next_nodes[next_node]
       
        if next_nodes is {}:
            raise ValueError(f"No path found from {current_node} to {end}")
        path.append(next_node)
        flow_path_.append(flow_)
        current_node = next_node
    return listOfTuples(path, flow_path_)
    #return path + flow_path_


#%%
#Explanation with example:
    #flow_value, flow_dict = nx.maximum_flow(FN_abs, _s = "Ifngr1", _t = "V_S", capacity = "capacity", flow_func= nx.algorithms.flow.dinitz)
    #The flow_dict["Ifngr1"] returns edges that goes out of Ifngr1, and the flow values you should flow through each edge.
    #And every other gene in the network is having a different flow in order to make the maximum flow of the whole network where the receptor is the source and the V_S is the target.
    #The max flow(flow_value) is the sum of these flow values and it also equals to the sum of the flow values that get into the target V_S at the end of the network.
#Without permutation test, raw flow value

List_Receptors_MoDC6H_Uniq___T_PcDC_LMoDC_T_PcDC = ["Il2rg","Csf1r","Icam1_R","Vcam1","Pde1b",
                                                    "Itga5","Sdc3","Ogfr","Plxnb2","Kdr","Tnfrsf1b",
                                                    "Traf2","Cxcr4","Robo3","Rpsa","Csf2rb","Tnfrsf12a",
                                                    "F11r","Anxa6","Hcar2","Ptafr","P2ry2","Atp5b"]
paths_max_flow = {}
paths_max_flow_directed = {}
directed_flag = False

#Create flow network
FN_abs = nx.from_pandas_edgelist(directed_mouse_PPI_Virtual_Sink_abs, "Source", "Target", edge_attr="capacity")#, create_using = nx.DiGraph)


flow_values_dict = {}
normalized_flow_values_dict = {}
receptors_left_out = []
receptors_paths_original = {}
flows_outdegree = {}

for receptor in List_Receptors_MoDC6H_Uniq___T_PcDC_LMoDC_T_PcDC:
    #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
    if (receptor in directed_mouse_PPI_Virtual_Sink_abs["Source"].tolist()): #or (receptor in directed_mouse_PPI_Virtual_Sink_abs["Target"].tolist()):
    #if (receptor in directed_mouse_PPI_Virtual_Sink_abs["Source"].tolist()) and (receptor in directed_mouse_PPI_Virtual_Sink_abs["Target"].tolist()):

        flow_value, flow_dict = nx.maximum_flow(FN_abs, _s = receptor, _t = "V_S", capacity = "capacity", flow_func= nx.algorithms.flow.dinitz)
        print(receptor.upper(), ":")
        print("Flow value =", flow_value)
        degree_centrality_dict = nx.degree_centrality(FN_abs)
        Receptor_dc = degree_centrality_dict[receptor]
        print("DC value =", Receptor_dc) #Degree Centrality
        normalized_flow_value = flow_value/Receptor_dc
        print("Normalized flow value=", normalized_flow_value, "\n\n")
        normalized_flow_values_dict[receptor] = normalized_flow_value
        
        flows_outdegree[receptor] = flow_dict[receptor]
      
        flow_values_dict[receptor] = flow_value
        if flow_value > 0:
            if directed_flag:
                paths_max_flow_directed[receptor] = extract_flow_path(flow_dict, receptor, "V_S")
            else:
                paths_max_flow[receptor] = extract_flow_path(flow_dict, receptor, "V_S")
    else:
        receptors_left_out.append(receptor)

# Create a list to store the data for the DataFrame
df_data = []

# Extract the receptor and the last element of each list
for receptor, values in paths_max_flow.items():
    last_tf = values[-1][0]  # Get the last element's first value (TF)
    df_data.append((receptor, last_tf))

# Create the DataFrame
df = pd.DataFrame(df_data, columns=["Receptor", "TF"])  

#%% 
'''Long runtime'''    
#Multi source max flow network. Source –inf-> receptors -> TFs --inf>Sink
#If the max flow value for the network is greater with a specific TF than without it, then the TF is important
directed_mouse_PPI_Virtual_Sink_abs_MultiSource = directed_mouse_PPI_Virtual_Sink_abs.copy()
list_of_V_Source_edges = []
TF_flow_WO_DF = pd.DataFrame({
                    "TF" : [],
                    "Flow_With_TF" : [],
                    "Flow_Without_TF" : [],
                    "Delta_With-Without" : []})
#CHECK
all(elem in list(directed_mouse_PPI_Virtual_Sink_abs_MultiSource['Source']) for elem in list_TFs_FDR005_Exist_In_PPI)


for receptor in flow_values_dict.keys():
    list_of_V_Source_edges.append(["V_Source", receptor, math.inf])

for V_Source_edge in list_of_V_Source_edges:
    directed_mouse_PPI_Virtual_Sink_abs_MultiSource.loc[len(directed_mouse_PPI_Virtual_Sink_abs_MultiSource)] = V_Source_edge

for TF in list_TFs_FDR005_Exist_In_PPI:
    FN_abs = nx.from_pandas_edgelist(directed_mouse_PPI_Virtual_Sink_abs_MultiSource, "Source", "Target", edge_attr="capacity")#, create_using = nx.DiGraph)
    With_flow_value, With_flow_dict = nx.maximum_flow(FN_abs, _s = "V_Source", _t = "V_S", capacity = "capacity", flow_func= nx.algorithms.flow.dinitz)
    
    directed_mouse_PPI_Virtual_Sink_abs_MultiSource_without_TF = directed_mouse_PPI_Virtual_Sink_abs_MultiSource[directed_mouse_PPI_Virtual_Sink_abs_MultiSource['Source'] != TF]
    FN_abs = nx.from_pandas_edgelist(directed_mouse_PPI_Virtual_Sink_abs_MultiSource_without_TF, "Source", "Target", edge_attr="capacity")#, create_using = nx.DiGraph)
    Without_flow_value, Without_flow_dict = nx.maximum_flow(FN_abs, _s = "V_Source", _t = "V_S", capacity = "capacity", flow_func= nx.algorithms.flow.dinitz)
    
    TF_flow_WO_DF.loc[len(TF_flow_WO_DF)] = [TF, With_flow_value, Without_flow_value, (With_flow_value-Without_flow_value)]
    
TF_flow_WO_DF = TF_flow_WO_DF.sort_values(by=['Delta_With-Without'], ascending=False)
TF_flow_WO_DF.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Multisource_Delta_Flow_Values.csv")

#%%
#Multisource Plots
#Bar plots:
TF_flow_WO_DF = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Multisource_Delta_Flow_Values.csv", index_col = 0)
# Create the bar plot with wider x-axis ticks
plt.figure(figsize=(12, 6))
plt.bar(TF_flow_WO_DF[TF_flow_WO_DF['Delta_With-Without'] > 1]['TF'], TF_flow_WO_DF[TF_flow_WO_DF['Delta_With-Without'] > 1]['Delta_With-Without'], color='skyblue')
plt.xlabel('TF')
plt.ylabel('Delta With-Without')
plt.title('Flow Value Delta With-Without for Each TF')
plt.xticks(rotation=90, fontsize=10)  # Rotate x-axis labels for better readability and increase font size
plt.tight_layout()  # Adjust layout to make room for x-axis labels
#plt.show()
plt.savefig('FlowValueDelta_With-Without_for_Each_TF.png', dpi = 300)

# Plotting the bar plot
plt.figure(figsize=(12, 6))

plt.bar(T_Cells_All_Norm_Counts.columns, T_Cells_All_Norm_Counts.loc["Rara"], color='skyblue')
plt.xticks(rotation=90, fontsize=10)  # Rotate x-axis labels for better readability and increase font size
# Adding labels and title
plt.xlabel('Conditions')
plt.ylabel('Rara Normalized Expression')
plt.title('Normalized Expression of Rara Gene Across Different CD8+ T Cells Conditions')
plt.tight_layout()  # Adjust layout to make room for x-axis labels
# Display the plot
#plt.show()
plt.savefig('Normalized_Expression_Rara_Gene_Across_CD8+_T_Cells_Conditions.png', dpi = 300)
#%%

#ENRICHMNET TEST HYPERGEOMETRIC DISTRIBUTION Subsetting table for HIGH FLOW TFs
filt_FDR_TF_DF_results005_HighFlow_TFs_Subset = filt_FDR_TF_DF_results005[filt_FDR_TF_DF_results005['TF'].isin(TF_flow_WO_DF[TF_flow_WO_DF['Delta_With-Without'] > 1]['TF'].tolist())]
filt_FDR_TF_DF_results005_HighFlow_TFs_Subset.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\\filt_FDR_TF_DF_results005_HighFlow_TFs_Subset.csv")
#%%

#ENRICHMNET TEST HYPERGEOMETRIC DISTRIBUTION IN ALL SAMPLES OF HIGH FLOW TFs

########Bulk RNA-seq table Raw Data ALL SAMPLES (To check enrichment of TFs):
T_Cells_Raw_Counts_All = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\T_Cells\T_Cells_Raw_Data.csv", index_col = 0)
#TRANS
T_Cells_Raw_Counts_All_trans = T_Cells_Raw_Counts_All.T

rna_seq_df = T_Cells_Raw_Counts_All_trans.copy()
#Reading Mouse TFs table
TF_Table_Mouse = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Table_Mouse.csv", index_col = 0)
TF_flow_WO_DF = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Multisource_Delta_Flow_Values.csv", index_col = 0)
TF_flow_WO_DF_high_flow = TF_flow_WO_DF[TF_flow_WO_DF['Delta_With-Without'] > 1]['TF'].tolist()
# Dictionary to store enrichment results for each TF
enrichment_results = {}

# Total number of genes in the expression table
total_genes = len(rna_seq_df.columns)

# Total number of differentially expressed genes across all samples

#DEG = for raw counts, expression > 10 in all samples
#.all() makes sure its greater than 10 in all samples(rows).
#sum() sums how many genes have expression greater than 10 in all samples.
total_diff_expr_genes = sum((rna_seq_df > 500).all())

for tf in TF_flow_WO_DF_high_flow:
    # Extract target genes for the current TF
    target_genes = TF_Table_Mouse[TF_Table_Mouse['Mouse_Source'] == tf]['Mouse_Target'].tolist()
    
    # Create a new list of column names that exist in the dataframe
    valid_columns = [col for col in target_genes if col in rna_seq_df.columns]
    
    # Count the number of differentially expressed target genes
    #.all() makes sure its greater than 10 in all samples(rows).
    #sum() sums how many genes have expression greater than 10 in all samples.
    diff_expr_target_genes = sum((rna_seq_df[valid_columns] > 500).all())
  
    # Number of target genes for the current TF
    num_target_genes = len(target_genes)
    
    # Perform hypergeometric test
        #M = The population size, total number of genes in the bulkRNAseq expression table.
        #n = The number of successes in the population, total number of genes that are differentially expressed in the table.
        #N = The sample size. Number of target genes we check for the particular TF.
        #k = The number of drawn successes, target genes that are differentially expressed.
    p_value = hypergeom.sf(M=total_genes, n=total_diff_expr_genes, N=num_target_genes,k=diff_expr_target_genes - 1)
    
    
    enrichment_results[tf] =   {
         'Target Gene': target_genes,
         'Number of Differentially Expressed Target Genes': diff_expr_target_genes,
         'Total Number of Target Genes': num_target_genes,
         'Hypergeometric p-value': p_value
     }

    
# Sort keys based on the pvalues of the nested dictionaries
sorted_enrichment_results = dict(sorted(enrichment_results.items(), key=lambda x: x[1]['Hypergeometric p-value']))

# Correct for multiple testing
corrected_p_values = multipletests([value['Hypergeometric p-value'] for value in sorted_enrichment_results.values()], method='bonferroni')[1]

#Adding the adjusted pvalue for each TF in the nested dict
i = 0
for key in sorted_enrichment_results.keys():
    if i == len(corrected_p_values):
        break
    else:
        sorted_enrichment_results[key]['Adjusted p-value'] = corrected_p_values[i]
        i+=1

#Creating Dataframe that would store the results from the nested dictionary
TF_DF_results = pd.DataFrame(columns = ['TF','Total_Num_Targets','Num_Targets_DE','p-value','Adj. p-value'])

#Filling the dataframe
for key in sorted_enrichment_results.keys():
    TF_DF_results.loc[len(TF_DF_results)] = [key,sorted_enrichment_results[key]['Total Number of Target Genes'], sorted_enrichment_results[key]['Number of Differentially Expressed Target Genes'] ,sorted_enrichment_results[key]['Hypergeometric p-value'],
                                           sorted_enrichment_results[key]['Adjusted p-value']]

#Saving to csv file
TF_DF_results.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_All_Samples_Results.csv")

#Filtering by pvalue < 0.05
filt_pval_TF_DF_results005 = TF_DF_results[TF_DF_results['p-value'] < 0.05]
#Filtering by adjusted p val < 0.05
filt_FDR_TF_DF_results005 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.05]
#Filtering by adjusted p val < 0.01
filt_FDR_TF_DF_results001 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.01]

#Saving to csv files
filt_pval_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_All_Samples_Results_Pval_lt_0.05.csv")
filt_FDR_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_All_Samples_Results_FDR_lt_0.05.csv")
filt_FDR_TF_DF_results001.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_All_Samples_Results_FDR_lt_0.01.csv")






#%%

#ENRICHMNET TEST HYPERGEOMETRIC DISTRIBUTION IN NAIVE SAMPLES OF HIGH FLOW TFs

########Bulk RNA-seq table Raw Data ALL SAMPLES (To check enrichment of TFs):
T_Cells_Raw_Counts_All = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\T_Cells\T_Cells_Raw_Data.csv", index_col = 0)
T_Cells_Raw_Counts_Naive = T_Cells_Raw_Counts_All[["Naive_T_cell_1", "Naive_T_cell_2", "Naive_T_cell_3", "Naive_T_cell_4"]]
#TRANS
T_Cells_Raw_Counts_All_trans = T_Cells_Raw_Counts_Naive.T

rna_seq_df = T_Cells_Raw_Counts_All_trans.copy()
#Reading Mouse TFs table
TF_Table_Mouse = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Table_Mouse.csv", index_col = 0)
TF_flow_WO_DF = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Multisource_Delta_Flow_Values.csv", index_col = 0)
TF_flow_WO_DF_high_flow = TF_flow_WO_DF[TF_flow_WO_DF['Delta_With-Without'] > 1]['TF'].tolist()
# Dictionary to store enrichment results for each TF
enrichment_results = {}

# Total number of genes in the expression table
total_genes = len(rna_seq_df.columns)

# Total number of differentially expressed genes across all samples

#DEG = for raw counts, expression > 10 in all samples
#.all() makes sure its greater than 10 in all samples(rows).
#sum() sums how many genes have expression greater than 10 in all samples.
total_diff_expr_genes = sum((rna_seq_df > 500).all())

for tf in TF_flow_WO_DF_high_flow:
    # Extract target genes for the current TF
    target_genes = TF_Table_Mouse[TF_Table_Mouse['Mouse_Source'] == tf]['Mouse_Target'].tolist()
    
    # Create a new list of column names that exist in the dataframe
    valid_columns = [col for col in target_genes if col in rna_seq_df.columns]
    
    # Count the number of differentially expressed target genes
    #.all() makes sure its greater than 10 in all samples(rows).
    #sum() sums how many genes have expression greater than 10 in all samples.
    diff_expr_target_genes = sum((rna_seq_df[valid_columns] > 500).all())
  
    # Number of target genes for the current TF
    num_target_genes = len(target_genes)
    
    # Perform hypergeometric test
        #M = The population size, total number of genes in the bulkRNAseq expression table.
        #n = The number of successes in the population, total number of genes that are differentially expressed in the table.
        #N = The sample size. Number of target genes we check for the particular TF.
        #k = The number of drawn successes, target genes that are differentially expressed.
    p_value = hypergeom.sf(M=total_genes, n=total_diff_expr_genes, N=num_target_genes,k=diff_expr_target_genes - 1)
    
    
    enrichment_results[tf] =   {
         'Target Gene': target_genes,
         'Number of Differentially Expressed Target Genes': diff_expr_target_genes,
         'Total Number of Target Genes': num_target_genes,
         'Hypergeometric p-value': p_value
     }

    
# Sort keys based on the pvalues of the nested dictionaries
sorted_enrichment_results = dict(sorted(enrichment_results.items(), key=lambda x: x[1]['Hypergeometric p-value']))

# Correct for multiple testing
corrected_p_values = multipletests([value['Hypergeometric p-value'] for value in sorted_enrichment_results.values()], method='bonferroni')[1]

#Adding the adjusted pvalue for each TF in the nested dict
i = 0
for key in sorted_enrichment_results.keys():
    if i == len(corrected_p_values):
        break
    else:
        sorted_enrichment_results[key]['Adjusted p-value'] = corrected_p_values[i]
        i+=1

#Creating Dataframe that would store the results from the nested dictionary
TF_DF_results = pd.DataFrame(columns = ['TF','Total_Num_Targets','Num_Targets_DE','p-value','Adj. p-value'])

#Filling the dataframe
for key in sorted_enrichment_results.keys():
    TF_DF_results.loc[len(TF_DF_results)] = [key,sorted_enrichment_results[key]['Total Number of Target Genes'], sorted_enrichment_results[key]['Number of Differentially Expressed Target Genes'] ,sorted_enrichment_results[key]['Hypergeometric p-value'],
                                           sorted_enrichment_results[key]['Adjusted p-value']]

#Saving to csv file
TF_DF_results.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_Naive_Samples_Results.csv")

#Filtering by pvalue < 0.05
filt_pval_TF_DF_results005 = TF_DF_results[TF_DF_results['p-value'] < 0.05]
#Filtering by adjusted p val < 0.05
filt_FDR_TF_DF_results005 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.05]
#Filtering by adjusted p val < 0.01
filt_FDR_TF_DF_results001 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.01]

#Saving to csv files
filt_pval_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_Naive_Samples_Results_Pval_lt_0.05.csv")
filt_FDR_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_Naive_Samples_Results_FDR_lt_0.05.csv")
filt_FDR_TF_DF_results001.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_Naive_Samples_Results_FDR_lt_0.01.csv")
#%%

#ENRICHMNET TEST HYPERGEOMETRIC DISTRIBUTION IN PcDC_LMoDC SAMPLES OF HIGH FLOW TFs

########Bulk RNA-seq table Raw Data ALL SAMPLES (To check enrichment of TFs):
T_Cells_Raw_Counts_All = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\T_Cells\T_Cells_Raw_Data.csv", index_col = 0)
T_Cells_Raw_Counts_PcDC_LMoDC = T_Cells_Raw_Counts_All[List_Of_Samples]
#TRANS
T_Cells_Raw_Counts_All_trans = T_Cells_Raw_Counts_PcDC_LMoDC.T

rna_seq_df = T_Cells_Raw_Counts_All_trans.copy()
#Reading Mouse TFs table
TF_Table_Mouse = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Table_Mouse.csv", index_col = 0)
TF_flow_WO_DF = pd.read_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\TF_Multisource_Delta_Flow_Values.csv", index_col = 0)
TF_flow_WO_DF_high_flow = TF_flow_WO_DF[TF_flow_WO_DF['Delta_With-Without'] > 1]['TF'].tolist()
# Dictionary to store enrichment results for each TF
enrichment_results = {}

# Total number of genes in the expression table
total_genes = len(rna_seq_df.columns)

# Total number of differentially expressed genes across all samples

#DEG = for raw counts, expression > 10 in all samples
#.all() makes sure its greater than 10 in all samples(rows).
#sum() sums how many genes have expression greater than 10 in all samples.
total_diff_expr_genes = sum((rna_seq_df > 500).all())

for tf in TF_flow_WO_DF_high_flow:
    # Extract target genes for the current TF
    target_genes = TF_Table_Mouse[TF_Table_Mouse['Mouse_Source'] == tf]['Mouse_Target'].tolist()
    
    # Create a new list of column names that exist in the dataframe
    valid_columns = [col for col in target_genes if col in rna_seq_df.columns]
    
    # Count the number of differentially expressed target genes
    #.all() makes sure its greater than 10 in all samples(rows).
    #sum() sums how many genes have expression greater than 10 in all samples.
    diff_expr_target_genes = sum((rna_seq_df[valid_columns] > 500).all())
  
    # Number of target genes for the current TF
    num_target_genes = len(target_genes)
    
    # Perform hypergeometric test
        #M = The population size, total number of genes in the bulkRNAseq expression table.
        #n = The number of successes in the population, total number of genes that are differentially expressed in the table.
        #N = The sample size. Number of target genes we check for the particular TF.
        #k = The number of drawn successes, target genes that are differentially expressed.
    p_value = hypergeom.sf(M=total_genes, n=total_diff_expr_genes, N=num_target_genes,k=diff_expr_target_genes - 1)
    
    
    enrichment_results[tf] =   {
         'Target Gene': target_genes,
         'Number of Differentially Expressed Target Genes': diff_expr_target_genes,
         'Total Number of Target Genes': num_target_genes,
         'Hypergeometric p-value': p_value
     }

    
# Sort keys based on the pvalues of the nested dictionaries
sorted_enrichment_results = dict(sorted(enrichment_results.items(), key=lambda x: x[1]['Hypergeometric p-value']))

# Correct for multiple testing
corrected_p_values = multipletests([value['Hypergeometric p-value'] for value in sorted_enrichment_results.values()], method='bonferroni')[1]

#Adding the adjusted pvalue for each TF in the nested dict
i = 0
for key in sorted_enrichment_results.keys():
    if i == len(corrected_p_values):
        break
    else:
        sorted_enrichment_results[key]['Adjusted p-value'] = corrected_p_values[i]
        i+=1

#Creating Dataframe that would store the results from the nested dictionary
TF_DF_results = pd.DataFrame(columns = ['TF','Total_Num_Targets','Num_Targets_DE','p-value','Adj. p-value'])

#Filling the dataframe
for key in sorted_enrichment_results.keys():
    TF_DF_results.loc[len(TF_DF_results)] = [key,sorted_enrichment_results[key]['Total Number of Target Genes'], sorted_enrichment_results[key]['Number of Differentially Expressed Target Genes'] ,sorted_enrichment_results[key]['Hypergeometric p-value'],
                                           sorted_enrichment_results[key]['Adjusted p-value']]

#Saving to csv file
TF_DF_results.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_PcDC_LMoDC_Samples_Results.csv")

#Filtering by pvalue < 0.05
filt_pval_TF_DF_results005 = TF_DF_results[TF_DF_results['p-value'] < 0.05]
#Filtering by adjusted p val < 0.05
filt_FDR_TF_DF_results005 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.05]
#Filtering by adjusted p val < 0.01
filt_FDR_TF_DF_results001 = TF_DF_results[TF_DF_results['Adj. p-value'] < 0.01]

#Saving to csv files
filt_pval_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_PcDC_LMoDC_Samples_Results_Pval_lt_0.05.csv")
filt_FDR_TF_DF_results005.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_PcDC_LMoDC_Samples_Results_FDR_lt_0.05.csv")
filt_FDR_TF_DF_results001.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Network_Calc\Enriched_High_Flow_TFs_In_T_Cell_Raw_PcDC_LMoDC_Samples_Results_FDR_lt_0.01.csv")








#%%

#FIXED that the TF->V_S are fixed and not permutated

#With permutation test:

import numpy.random
#Table of receptor- numpy array of (flow_value/Receptor_dc) values and then compare how many values are bigger than the (flow_value/Receptor_dc) value in normalized_flow_values_dict and thus calculate p value from it. (empirical)
receptors_left_out = []
receptors_existed = []
directed_mouse_PPI_Virtual_Sink_abs_perm = directed_mouse_PPI_Virtual_Sink_abs.copy()
flow_values_dict_permutated = {}
for receptor in List_Receptors_MoDC6H_Uniq___T_PcDC_LMoDC_T_PcDC:
    #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
    if (receptor in directed_mouse_PPI_Virtual_Sink_abs_perm["Source"].tolist()): #or (receptor in directed_mouse_PPI_Virtual_Sink_abs_perm["Target"].tolist()):
            flow_values_dict_permutated[receptor] = np.array([])
            receptors_existed.append(receptor)
    else:
        receptors_left_out.append(receptor)




np.random.seed(5)
#np.random.seed(42)
#np.random.seed(1)

# Get the indices of rows where the target is not "V_S"
non_vs_indices = directed_mouse_PPI_Virtual_Sink_abs_perm[directed_mouse_PPI_Virtual_Sink_abs_perm['Target'] != 'V_S'].index

#For random.shuffle:
#non_vs_indices = list(non_vs_indices)

# Permute only the rows where the target is not "V_S"
#permuted_indices = np.random.permutation(non_vs_indices)
# Update the "Target" and "capacity" columns with the permuted values
#directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'Target'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'Target'].values
#directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'capacity'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'capacity'].values


permuted_indices = np.random.permutation(non_vs_indices)
#import random

count_perm = 0
while count_perm < 200:
    print(count_perm)
    # Permute only the rows where the target is not "V_S"
    permuted_indices = np.random.permutation(permuted_indices) # np.random.permutation(non_vs_indices)
    #permuted_indices = non_vs_indices.copy()
    #random.shuffle(permuted_indices)
    
    
    # Update the "Target" and "capacity" columns with the permuted values
    directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'Target'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'Target'].values
    permuted_indices = np.random.permutation(permuted_indices) # np.random.permutation(non_vs_indices)
    directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'capacity'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'capacity'].values

    FN_abs = nx.from_pandas_edgelist(directed_mouse_PPI_Virtual_Sink_abs_perm, "Source", "Target", edge_attr="capacity")#, create_using = nx.DiGraph)
    
    for receptor in receptors_existed:
        #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
        
        flow_value, flow_dict = nx.maximum_flow(FN_abs, _s = receptor, _t = "V_S", capacity = "capacity", flow_func=nx.algorithms.flow.dinitz)
        print(receptor, "-", flow_value)
        flow_values_dict_permutated[receptor] = np.append(flow_values_dict_permutated[receptor],flow_value) #flow_value
        
    count_perm += 1




    
# Dictionary to store p-values for each receptor
p_values = {}

for receptor in flow_values_dict.keys():
    print(receptor)
    original_value = flow_values_dict[receptor]
    print(original_value)
    permuted_values = flow_values_dict_permutated[receptor]
    
    # Count how many permuted values are greater than or equal to the original value
    count_greater_or_equal = np.sum(original_value <= permuted_values)
    
    # Calculate p-value as the fraction of permuted values that are greater than or equal to the original value
    p_value = count_greater_or_equal / len(permuted_values)
    
    # Store p-value in the dictionary
    p_values[receptor] = p_value
    
from scipy import stats
p_values_1samp_ttest = {}

for receptor in flow_values_dict.keys():
    print(receptor)
    original_value = flow_values_dict[receptor]
    print(original_value)
    permuted_values = flow_values_dict_permutated[receptor]
    
    # Store p-value in the dictionary
    p_values_1samp_ttest[receptor] = stats.ttest_1samp(flow_values_dict_permutated[receptor], flow_values_dict[receptor])[1]








#%%
#Preparing dataframe for circos:
#Genes with p value above 0.05 would be gray
#Genes with p value below 0.05 would be shades of purple


receptor_names = list(p_values_1samp_ttest.keys())
receptors_Pvals = list(p_values_1samp_ttest.values())
fdr_corrected_p_values = multipletests(receptors_Pvals, method='bonferroni')[1]

# Update the dictionary with FDR values
for gene, fdr_corrected in zip(receptor_names, fdr_corrected_p_values):
    p_values_1samp_ttest[gene] = {
        'DSA': normalized_flow_values_dict[gene],
        'p_value': p_values_1samp_ttest[gene],
        'fdr_value': fdr_corrected,
        'n': len(flows_outdegree[gene])
    }
    

final_DSA_DF = pd.DataFrame([
    {'Receptor': gene, 'DSA': values['DSA'],'p_value': values['p_value'], 'fdr_value': values['fdr_value'], 'n': values['n']}
    for gene, values in p_values_1samp_ttest.items()
])

final_DSA_DF.set_index(final_DSA_DF.Receptor, inplace=True, drop =False)
final_DSA_DF.index.names = [None]
#final_DSA_DF["n"] = 3
final_DSA_DF.to_csv("C:\Kfir_Thesis_Asaf_Madi\Project_2 - BulkRNAseq\Circos_Code\DSA_Flow_Calculation_DF.csv")
#%%






































































































































#%%

#WITH ERROR HANDLING
#FIX that the TF->V_S are fixed and not permutated

#With permutation test:

import numpy.random
#Table of receptor- numpy array of (flow_value/Receptor_dc) values and then compare how many values are bigger than the (flow_value/Receptor_dc) value in normalized_flow_values_dict and thus calculate p value from it. (empirical)
receptors_left_out = []
receptors_existed = []
directed_mouse_PPI_Virtual_Sink_abs_perm = directed_mouse_PPI_Virtual_Sink_abs.copy()
flow_values_dict_permutated = {}
for receptor in List_Receptors_MoDC6H_Uniq___T_PcDC_LMoDC_T_PcDC:
    #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
    if (receptor in directed_mouse_PPI_Virtual_Sink_abs_perm["Source"].tolist()) or (receptor in directed_mouse_PPI_Virtual_Sink_abs_perm["Target"].tolist()):
            flow_values_dict_permutated[receptor] = np.array([])
            receptors_existed.append(receptor)
    else:
        receptors_left_out.append(receptor)
count_perm = 0
flow_dict_perm_temp = {}
perm_succesful = False


#np.random.seed(42)
np.random.seed(1)

# Get the indices of rows where the target is not "V_S"
non_vs_indices = directed_mouse_PPI_Virtual_Sink_abs_perm[directed_mouse_PPI_Virtual_Sink_abs_perm['Target'] != 'V_S'].index
# Permute only the rows where the target is not "V_S"
permuted_indices = np.random.permutation(non_vs_indices)
# Update the "Target" and "capacity" columns with the permuted values
directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'Target'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'Target'].values
directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'capacity'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'capacity'].values

while count_perm < 100:
    flow_dict_perm_temp = {}
    print(count_perm)
    #FIX - first check if you can save flow values for all receptors without exception and if you can then move it to the dictionary
    try:
        # Permute only the rows where the target is not "V_S"
        permuted_indices = np.random.permutation(non_vs_indices)
        # Update the "Target" and "capacity" columns with the permuted values
        directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'Target'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'Target'].values
        directed_mouse_PPI_Virtual_Sink_abs_perm.loc[non_vs_indices, 'capacity'] = directed_mouse_PPI_Virtual_Sink_abs_perm.loc[permuted_indices, 'capacity'].values

        FN_abs = nx.from_pandas_edgelist(directed_mouse_PPI_Virtual_Sink_abs_perm, "Source", "Target", edge_attr="capacity")
        
        for receptor in receptors_existed:
            #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
            flow_value, flow_dict = nx.maximum_flow(FN_abs, _s = receptor, _t = "V_S", capacity = "capacity", flow_func=nx.algorithms.flow.dinitz)
            flow_dict_perm_temp[receptor] = flow_value
                    
        perm_succesful = True  
        if perm_succesful:
            
            for key in flow_dict_perm_temp.keys():
                print(flow_dict_perm_temp[key])
                flow_values_dict_permutated[key] = np.append(flow_values_dict_permutated[key],flow_dict_perm_temp[key])
   
        count_perm += 1
        
    except nx.NetworkXUnbounded:
        perm_succesful = False
        continue


# Dictionary to store p-values for each receptor
p_values = {}

for receptor in flow_values_dict.keys():
    print(receptor)
    original_value = flow_values_dict[receptor]
    print(original_value)
    permuted_values = flow_values_dict_permutated[receptor]
    
    # Count how many permuted values are greater than or equal to the original value
    count_greater_or_equal = np.sum(original_value <= permuted_values)
    
    # Calculate p-value as the fraction of permuted values that are greater than or equal to the original value
    p_value = count_greater_or_equal / len(permuted_values)
    
    # Store p-value in the dictionary
    p_values[receptor] = p_value
    
    


#%%
#With permutation test:

import numpy.random
#Table of receptor- numpy array of (flow_value/Receptor_dc) values and then compare how many values are bigger than the (flow_value/Receptor_dc) value in normalized_flow_values_dict and thus calculate p value from it. (empirical)
receptors_left_out = []
receptors_existed = []
directed_mouse_PPI_Virtual_Sink_abs_perm = directed_mouse_PPI_Virtual_Sink_abs.copy()
flow_values_dict_permutated = {}
for receptor in List_Receptors_MoDC6H_Uniq___T_PcDC_LMoDC_T_PcDC:
    #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
    if (receptor in directed_mouse_PPI_Virtual_Sink_abs_perm["Source"].tolist()) or (receptor in directed_mouse_PPI_Virtual_Sink_abs_perm["Target"].tolist()):
            flow_values_dict_permutated[receptor] = np.array([])
            receptors_existed.append(receptor)
    else:
        receptors_left_out.append(receptor)
count_perm = 0
flow_dict_perm_temp = {}
perm_succesful = False


#np.random.seed(42)
np.random.seed(1)

permuted_indices = numpy.random.permutation(len(directed_mouse_PPI_Virtual_Sink_abs_perm))
directed_mouse_PPI_Virtual_Sink_abs_perm["Target"] = directed_mouse_PPI_Virtual_Sink_abs_perm["Target"].values[permuted_indices]
directed_mouse_PPI_Virtual_Sink_abs_perm["capacity"] = directed_mouse_PPI_Virtual_Sink_abs_perm["capacity"].values[permuted_indices]


while count_perm < 400:
    flow_dict_perm_temp = {}
    print(count_perm)
    #FIX - first check if you can save flow values for all receptors without exception and if you can then move it to the dictionary
    try:
        # Permute only the "Target" and "Score" columns
        #permuted_indices = numpy.random.permutation(len(directed_mouse_PPI_Virtual_Sink_abs))
        permuted_indices = numpy.random.permutation(len(directed_mouse_PPI_Virtual_Sink_abs_perm))
        directed_mouse_PPI_Virtual_Sink_abs_perm["Target"] = directed_mouse_PPI_Virtual_Sink_abs_perm["Target"].values[permuted_indices]
        directed_mouse_PPI_Virtual_Sink_abs_perm["capacity"] = directed_mouse_PPI_Virtual_Sink_abs_perm["capacity"].values[permuted_indices]
        
        FN_abs = nx.from_pandas_edgelist(directed_mouse_PPI_Virtual_Sink_abs_perm, "Source", "Target", edge_attr="capacity")
        
        for receptor in receptors_existed:
            #Chosen receptor from Circos plot MoDC_6H unique vs T_PcDC_LMoDC/T_PcDC
            flow_value, flow_dict = nx.maximum_flow(FN_abs, _s = receptor, _t = "V_S", capacity = "capacity", flow_func=nx.algorithms.flow.dinitz)
            flow_dict_perm_temp[receptor] = flow_value
                    
        perm_succesful = True  
        if perm_succesful:
            
            for key in flow_dict_perm_temp.keys():
                print(flow_dict_perm_temp[key])
                flow_values_dict_permutated[key] = np.append(flow_values_dict_permutated[key],flow_dict_perm_temp[key])
   
        count_perm += 1
        
    except nx.NetworkXUnbounded:
        perm_succesful = False
        continue


# Dictionary to store p-values for each receptor
p_values = {}

for receptor in flow_values_dict.keys():
    print(receptor)
    original_value = flow_values_dict[receptor]
    print(original_value)
    permuted_values = flow_values_dict_permutated[receptor]
    
    # Count how many permuted values are greater than or equal to the original value
    count_greater_or_equal = np.sum(original_value >= permuted_values)
    
    # Calculate p-value as the fraction of permuted values that are greater than or equal to the original value
    p_value = count_greater_or_equal / len(permuted_values)
    
    # Store p-value in the dictionary
    p_values[receptor] = p_value
    
    


