import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split


df = pd.read_csv('../data/after_interpolate.csv')

# 同じ人はtrain/test内でセットになるように分割
unique_df = df[['加入者id','is_mac']].drop_duplicates()

# 加入者idを4:1に分割してtrain_idsとtest_idsを作成
X = unique_df['加入者id']
y = unique_df['is_mac']
train_ids, test_ids, y_train, y_test = train_test_split(X, y,test_size=0.2, random_state=42,stratify=y)

# 元のデータをtrain_idsとtest_idsに基づいて分割
train_data = df[df['加入者id'].isin(train_ids)]
test_data = df[df['加入者id'].isin(test_ids)]

train_data.to_csv('../data/train_data.csv',index=False)
test_data.to_csv('../data/test_data.csv',index=False)

# is_mac列を置換
mapping = {1: True, 0: False}
train_data['is_mac'] = train_data['is_mac'].replace(mapping).astype(bool)
test_data['is_mac'] = test_data['is_mac'].replace(mapping).astype(bool)

###### train ######
X_train = train_data.drop(columns=['加入者id','is_mac','month_maculopathy']) # 不要な列を削除
y_train = train_data[["is_mac","month_maculopathy"]]
yy_train = y_train.to_records(index=False)

###### test ######
X_test = test_data.drop(columns=['加入者id','is_mac','month_maculopathy']) # 不要な列を削除
y_test = test_data[["is_mac","month_maculopathy"]]
yy_test = y_test.to_records(index=False)

from sksurv.ensemble import RandomSurvivalForest
import _pickle as cPickle
from sksurv.metrics import integrated_brier_score, cumulative_dynamic_auc
import logging
import os
import copy
import numpy as np

n_rep = 10
path = '/work/kawakami/Project/Kimitsu_Diabetes/analysis'

array_risk_scores = []
array_metric = []

for i in range(n_rep):

	print(i)
	
	with open(os.path.join(f"{path}","RSF_models","RSF_maculopathy_risk_rep"+str(i+1)+".sav"), 'rb') as f:
		rsf = cPickle.load(f)

	rsf_risk_scores = rsf.predict(X_test)

	# Survival probability
	df_risk_scores = test_data[['加入者id','is_mac','month_maculopathy']].reset_index(drop=True).copy()
	df_risk_scores["risk_score"] = rsf_risk_scores
	array_risk_scores.append(df_risk_scores)

	array_cindex = []
	array_IBS = []
	array_mean_auc = []

	for t in range(12, 144, 12):
		test_data_int = test_data.loc[test_data["month_maculopathy"] >= (t-12),:]
		X_test_int = test_data_int.drop(columns=['加入者id','is_mac','month_maculopathy']) # 不要な列を削除
		y_test_int = test_data_int[["is_mac","month_maculopathy"]]
		times = np.arange(min(test_data_int["month_maculopathy"]), max(test_data_int["month_maculopathy"]))
		yy_test_int = y_test_int.to_records(index=False)
		array_cindex.append(rsf.score(X_test_int, yy_test_int))
		survs_int = rsf.predict_survival_function(X_test_int)
		risk_scores_int = rsf.predict(X_test_int)
		preds_int = np.asarray([[fn(t) for t in times] for fn in survs_int])
		array_IBS.append(integrated_brier_score(yy_train, yy_test_int, preds_int, times))
		auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test_int, risk_scores_int, times)
		array_mean_auc.append(mean_auc)
	
	# test_data_int = test_data.loc[test_data["month_maculopathy"] >= 120,:]
	# X_test_int = test_data_int.drop(columns=['加入者id','is_mac','month_maculopathy']) # 不要な列を削除
	# y_test_int = test_data_int[["is_mac","month_maculopathy"]]
	# yy_test_int = y_test_int.to_records(index=False)
	# times = np.arange(min(test_data_int["month_maculopathy"]), max(test_data_int["month_maculopathy"]))
	# array_cindex.append(rsf.score(X_test_int, yy_test_int))
	# survs_int = rsf.predict_survival_function(X_test_int)
	# risk_scores_int = rsf.predict(X_test_int)
	# preds_int = np.asarray([[fn(t) for t in times] for fn in survs_int])
	# array_IBS.append(integrated_brier_score(yy_train, yy_test_int, preds_int, times))
	# auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test_int, risk_scores_int, times)
	# array_mean_auc.append(mean_auc)

	df_metric = pd.DataFrame()
	df_metric["time"] = list(range(0, 132, 12))
	df_metric["cindex"] = array_cindex
	df_metric["IBS"] = array_IBS
	df_metric["meanAUC"] = array_mean_auc

	array_metric.append(df_metric)

df_psfs = pd.concat(array_risk_scores)
df_psfs_mean = df_psfs.groupby(['加入者id','month_maculopathy']).mean().reset_index()
df_psfs_mean.to_csv(os.path.join(f"{path}","RSF_models_maculopathy_risk_score_mean.csv"),index=False,sep="\t")

df_metric = pd.concat(array_metric)
df_metric_mean = df_metric.groupby(['time']).mean().reset_index()
df_metric_mean.to_csv(os.path.join(f"{path}","RSF_models_maculopathy_interval_metrics.csv"),index=False,sep="\t")
