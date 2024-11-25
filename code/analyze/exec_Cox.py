from sksurv.linear_model import CoxPHSurvivalAnalysis
import _pickle as cPickle
from sksurv.metrics import integrated_brier_score, cumulative_dynamic_auc
import logging
import os
import copy
import pandas as pd
import numpy as np


n_rep = 10

train_data = pd.read_csv('train_data.txt', sep='\t')
test_data = pd.read_csv('test_data.txt', sep='\t')

mapping = {1: True, 0: False}
train_data['is_mac'] = train_data['is_mac'].replace(mapping).astype(bool)
test_data['is_mac'] = test_data['is_mac'].replace(mapping).astype(bool)

###### train ######
X_train = train_data.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_train = train_data[["is_mac","month_maculopathy"]]
yy_train = y_train.to_records(index=False)

###### test ######
X_test = test_data.drop(columns=['加入者id','is_mac','month_maculopathy'])
y_test = test_data[["is_mac","month_maculopathy"]]
yy_test = y_test.to_records(index=False)

cindex_test = []
IBS_test = []
cAUC_test = []

feature_names = X_train.columns.tolist()

for i in range(n_rep):
	print('CoxPH start',i)

	cph = CoxPHSurvivalAnalysis(alpha=0)
	cph.fit(X_train, yy_train)
	with open("./Cox_models/CoxPH_maculopathy_risk_rep"+str(i+1)+".sav", 'wb') as f:
		cPickle.dump(cph, f)

	cindex_test.append(cph.score(X_test, yy_test))
	
	survs = cph.predict_survival_function(X_test)
	times = np.arange(12, 120)
	preds = np.asarray([[fn(t) for t in times] for fn in survs])
	IBS_test.append(integrated_brier_score(yy_train, yy_test, preds, times))

	auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test, preds, times)

	cAUC_test.append(mean_auc)

metric_summary = pd.DataFrame([])
metric_summary["cindex"] = cindex_test
metric_summary["IBS"] = IBS_test
metric_summary["cAUC"] = cAUC_test
metric_summary.to_csv("CoxPH_maculopathy_risk_performance_metrics.txt", sep='\t',index=False)
