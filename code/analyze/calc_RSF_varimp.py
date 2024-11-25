from sksurv.ensemble import RandomSurvivalForest
import _pickle as cPickle
from sklearn.inspection import permutation_importance
import logging
import os
import copy
import pandas as pd
import numpy as np


top_feature_num = 20

n_rep = 10

path = '/work/kawakami/Project/Kimitsu_Diabetes/analysis'

data_all = pd.read_csv("before_interpolate.csv")
DME_cluster = pd.read_csv("RSF_maculopathy_risk_score_DME_cluster.txt",sep="\t")

cluster1_ID = list(DME_cluster.loc[DME_cluster["cluster"]==1,"Patient_ID"])
cluster2_ID = list(DME_cluster.loc[DME_cluster["cluster"]==2,"Patient_ID"])
cluster3_ID = list(DME_cluster.loc[DME_cluster["cluster"]==3,"Patient_ID"])
nonDME_ID = list(set(data_all.loc[data_all["is_mac"]==0,"加入者id"]))

train_data = pd.read_csv('../data/train_data.csv')
test_data = pd.read_csv('../data/test_data.csv')

mapping = {1: True, 0: False}
train_data['is_mac'] = train_data['is_mac'].replace(mapping).astype(bool)
test_data['is_mac'] = test_data['is_mac'].replace(mapping).astype(bool)

test_data_cluster1 = test_data.loc[test_data["加入者id"].isin(cluster1_ID),:]
test_data_cluster2 = test_data.loc[test_data["加入者id"].isin(cluster2_ID),:]
test_data_cluster3 = test_data.loc[test_data["加入者id"].isin(cluster3_ID),:]
test_data_nonDME = test_data.loc[test_data["加入者id"].isin(nonDME_ID),:]

###### train ######
X_train = train_data.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_train = train_data[["is_mac","month_maculopathy"]]
yy_train = y_train.to_records(index=False)

###### test ######
X_test = test_data.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_test = test_data[["is_mac","month_maculopathy"]]
yy_test = y_test.to_records(index=False)

X_test_cluster1 = test_data_cluster1.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_test_cluster1 = test_data_cluster1[["is_mac","month_maculopathy"]]
yy_test_cluster1 = y_test_cluster1.to_records(index=False)

X_test_cluster2 = test_data_cluster2.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_test_cluster2 = test_data_cluster2[["is_mac","month_maculopathy"]]
yy_test_cluster2 = y_test_cluster2.to_records(index=False)

X_test_cluster3 = test_data_cluster3.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_test_cluster3 = test_data_cluster3[["is_mac","month_maculopathy"]]
yy_test_cluster3 = y_test_cluster3.to_records(index=False)

# X_test_nonDME = test_data_nonDME.drop(columns=['加入者id','is_mac','month_maculopathy'])
# y_test_nonDME = test_data_nonDME[["is_mac","month_maculopathy"]]
# yy_test_nonDME = y_test_nonDME.to_records(index=False)

imp_summary = pd.DataFrame([])
imp_summary_cluster1 = pd.DataFrame([])
imp_summary_cluster2 = pd.DataFrame([])
imp_summary_cluster3 = pd.DataFrame([])
# imp_summary_nonDME = pd.DataFrame([])

feature_names = X_train.columns.tolist()

for i in range(n_rep):
	print(str(i+1))

	with open(os.path.join(f"{path}","RSF_models","RSF_maculopathy_risk_rep"+str(i+1)+".sav"), 'rb') as f:
		rsf = cPickle.load(f)
	
	print("All")
	perm = permutation_importance(rsf,X_test, yy_test, n_repeats=5, random_state=i)
	imp_summary = pd.concat([imp_summary, pd.DataFrame(perm.importances_mean)], axis=1)
	print("cluster1")	
	perm = permutation_importance(rsf,X_test_cluster1, yy_test_cluster1, n_repeats=5, random_state=i)
	imp_summary_cluster1 = pd.concat([imp_summary_cluster1, pd.DataFrame(perm.importances_mean)], axis=1)
	print("cluster2")
	perm = permutation_importance(rsf,X_test_cluster2, yy_test_cluster2, n_repeats=5, random_state=i)
	imp_summary_cluster2 = pd.concat([imp_summary_cluster2, pd.DataFrame(perm.importances_mean)], axis=1)
	print("cluster3")
	perm = permutation_importance(rsf,X_test_cluster3, yy_test_cluster3, n_repeats=5, random_state=i)
	imp_summary_cluster3 = pd.concat([imp_summary_cluster3, pd.DataFrame(perm.importances_mean)], axis=1)
	# print("nonDME")
	# perm = permutation_importance(rsf,X_test_nonDME, yy_test_nonDME, n_repeats=5, random_state=i)
	# imp_summary_nonDME = pd.concat([imp_summary_nonDME, pd.DataFrame(perm.importances_mean)], axis=1)
	
# all
imp_summary.columns = [ "rep_"+str(i+1) for i in range(n_rep)]
imp_summary.index = X_train.columns
imp_summary = imp_summary.reset_index()
imp_summary_tops = copy.deepcopy(imp_summary).set_index('index')
imp_summary_tops["median"] = imp_summary_tops.median(axis="columns")
imp_summary_tops = imp_summary_tops.sort_values("median",ascending=False)
imp_summary_tops.to_csv(os.path.join(f"{path}","RSF_maculopathy_risk_varimp_descending.txt"),sep="\t")

# cluster1
imp_summary_cluster1.columns = [ "rep_"+str(i+1) for i in range(n_rep)]
imp_summary_cluster1.index = X_train.columns
imp_summary_cluster1 = imp_summary_cluster1.reset_index()
imp_summary_cluster1_tops = copy.deepcopy(imp_summary_cluster1).set_index('index')
imp_summary_cluster1_tops["median"] = imp_summary_cluster1_tops.median(axis="columns")
imp_summary_cluster1_tops = imp_summary_cluster1_tops.sort_values("median",ascending=False)
imp_summary_cluster1_tops.to_csv(os.path.join(f"{path}","RSF_maculopathy_risk_varimp_cluster1_descending.txt"),sep="\t")

# cluster2
imp_summary_cluster2.columns = [ "rep_"+str(i+1) for i in range(n_rep)]
imp_summary_cluster2.index = X_train.columns
imp_summary_cluster2 = imp_summary_cluster2.reset_index()
imp_summary_cluster2_tops = copy.deepcopy(imp_summary_cluster2).set_index('index')
imp_summary_cluster2_tops["median"] = imp_summary_cluster2_tops.median(axis="columns")
imp_summary_cluster2_tops = imp_summary_cluster2_tops.sort_values("median",ascending=False)
imp_summary_cluster2_tops.to_csv(os.path.join(f"{path}","RSF_maculopathy_risk_varimp_cluster2_descending.txt"),sep="\t")

# cluster3
imp_summary_cluster3.columns = [ "rep_"+str(i+1) for i in range(n_rep)]
imp_summary_cluster3.index = X_train.columns
imp_summary_cluster3 = imp_summary_cluster3.reset_index()
imp_summary_cluster3_tops = copy.deepcopy(imp_summary_cluster3).set_index('index')
imp_summary_cluster3_tops["median"] = imp_summary_cluster3_tops.iloc[:,1:].median(axis="columns")
imp_summary_cluster3_tops = imp_summary_cluster3_tops.sort_values("median",ascending=False)
imp_summary_cluster3_tops.to_csv(os.path.join(f"{path}","RSF_maculopathy_risk_varimp_cluster3_descending.txt"),sep="\t")

# # nonDME
# imp_summary_nonDME.columns = [ "rep_"+str(i+1) for i in range(n_rep)]
# imp_summary_nonDME.index = X_train.columns
# imp_summary_nonDME = imp_summary_nonDME.reset_index()
# imp_summary_nonDME_tops = copy.deepcopy(imp_summary_nonDME).set_index('index')
# imp_summary_nonDME_tops["median"] = imp_summary_nonDME_tops.median(axis="columns")
# imp_summary_nonDME_tops = imp_summary_nonDME_tops.sort_values("median",ascending=False)
# imp_summary_nonDME_tops.to_csv(os.path.join(f"{path}","RSF_models","RSF_maculopathy_risk_varimp_nonDME_descending.txt",sep="\t"))