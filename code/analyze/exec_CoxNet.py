from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.linear_model import CoxPHSurvivalAnalysis
import _pickle as cPickle
from sksurv.metrics import integrated_brier_score, cumulative_dynamic_auc
import logging
import os
import copy
import pandas as pd
import numpy as np

import warnings

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.exceptions import FitFailedWarning

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

cindex_LASSO_test = []
cindex_RIDGE_test = []
cindex_CoxNet_test = []

IBS_LASSO_test = []
IBS_RIDGE_test = []
IBS_CoxNet_test = []

cAUC_LASSO_test = []
cAUC_RIDGE_test = []
cAUC_CoxNet_test = []

feature_names = X_train.columns.tolist()

for i in range(n_rep):
	print('CoxNet start',i)

	coxnet_pipe = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=1, alpha_min_ratio=0.01, max_iter=100))
	warnings.simplefilter("ignore", UserWarning)
	warnings.simplefilter("ignore", FitFailedWarning)
	coxnet_pipe.fit(X_train, yy_train)

	alphas = 10.0 ** np.linspace(-4, 4, 50)
	cv = KFold(n_splits=5, shuffle=True, random_state=0)
	gcv = GridSearchCV(
	    make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=1)),
	    param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in alphas]},
	    cv=cv,
	    error_score=0.5,
	    n_jobs=16,
	).fit(X_train, yy_train)

	cnet = CoxnetSurvivalAnalysis(l1_ratio=1,max_iter=100,alphas=[gcv.best_params_["coxnetsurvivalanalysis__alphas"][0]],fit_baseline_model=True)
	cnet.fit(X_train, yy_train)
	with open("./Cox_models/CoxLASSO_maculopathy_risk_rep"+str(i+1)+".sav", 'wb') as f:
		cPickle.dump(cnet, f)

	cindex_LASSO_test.append(cnet.score(X_test, yy_test))

	survs = cnet.predict_survival_function(X_test)
	risk_scores = cnet.predict(X_test)
	times = np.arange(12, 120)
	preds = np.asarray([[fn(t) for t in times] for fn in survs])
	IBS_LASSO_test.append(integrated_brier_score(yy_train, yy_test, preds, times))

	auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test, risk_scores, times)

	cAUC_LASSO_test.append(mean_auc)

	coxnet_pipe = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=1e-16, alpha_min_ratio=0.01, max_iter=100))
	warnings.simplefilter("ignore", UserWarning)
	warnings.simplefilter("ignore", FitFailedWarning)
	# coxnet_pipe.fit(X_train, yy_train)

	cv = KFold(n_splits=5, shuffle=True, random_state=0)
	gcv = GridSearchCV(
	    make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=1e-16)),
	    param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in alphas]},
	    cv=cv,
	    error_score=0.5,
	    n_jobs=16,
	).fit(X_train, yy_train)

	cnet = CoxPHSurvivalAnalysis(alpha=gcv.best_params_["coxnetsurvivalanalysis__alphas"][0])
	cnet.fit(X_train, yy_train)
	with open("./Cox_models/CoxRIDGE_maculopathy_risk_rep"+str(i+1)+".sav", 'wb') as f:
		cPickle.dump(cnet, f)

	cindex_RIDGE_test.append(cnet.score(X_test, yy_test))

	survs = cnet.predict_survival_function(X_test)
	risk_scores = cnet.predict(X_test)
	times = np.arange(12, 120)
	preds = np.asarray([[fn(t) for t in times] for fn in survs])
	IBS_RIDGE_test.append(integrated_brier_score(yy_train, yy_test, preds, times))

	auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test, risk_scores, times)

	cAUC_RIDGE_test.append(mean_auc)

	coxnet_pipe = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.5, alpha_min_ratio=0.01, max_iter=100))
	warnings.simplefilter("ignore", UserWarning)
	warnings.simplefilter("ignore", FitFailedWarning)
	coxnet_pipe.fit(X_train, yy_train)

	estimated_alphas = coxnet_pipe.named_steps["coxnetsurvivalanalysis"].alphas_
	cv = KFold(n_splits=5, shuffle=True, random_state=0)
	gcv = GridSearchCV(
	    make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.5)),
	    param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in alphas]},
	    cv=cv,
	    error_score=0.5,
	    n_jobs=16,
	).fit(X_train, yy_train)

	cnet = CoxnetSurvivalAnalysis(l1_ratio=0.5,max_iter=100,alphas=[gcv.best_params_["coxnetsurvivalanalysis__alphas"][0]],fit_baseline_model=True)
	cnet.fit(X_train, yy_train)
	with open("./Cox_models/CoxNet0.5_maculopathy_risk_rep"+str(i+1)+".sav", 'wb') as f:
		cPickle.dump(cnet, f)

	cindex_CoxNet_test.append(cnet.score(X_test, yy_test))

	survs = cnet.predict_survival_function(X_test)
	risk_scores = cnet.predict(X_test)
	times = np.arange(12, 120)
	preds = np.asarray([[fn(t) for t in times] for fn in survs])
	IBS_CoxNet_test.append(integrated_brier_score(yy_train, yy_test, preds, times))

	auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test, risk_scores, times)

	cAUC_CoxNet_test.append(mean_auc)


metric_LASSO_summary = pd.DataFrame([])
metric_LASSO_summary["cindex"] = cindex_LASSO_test
metric_LASSO_summary["IBS"] = IBS_LASSO_test
metric_LASSO_summary["cAUC"] = cAUC_LASSO_test
metric_LASSO_summary.to_csv("CoxLASSO_maculopathy_risk_performance_metrics.txt", sep='\t',index=False)

metric_RIDGE_summary = pd.DataFrame([])
metric_RIDGE_summary["cindex"] = cindex_RIDGE_test
metric_RIDGE_summary["IBS"] = IBS_RIDGE_test
metric_RIDGE_summary["cAUC"] = cAUC_RIDGE_test
metric_RIDGE_summary.to_csv("CoxRIDGE_maculopathy_risk_performance_metrics.txt", sep='\t',index=False)

metric_CoxNet_summary = pd.DataFrame([])
metric_CoxNet_summary["cindex"] = cindex_CoxNet_test
metric_CoxNet_summary["IBS"] = IBS_CoxNet_test
metric_CoxNet_summary["cAUC"] = cAUC_CoxNet_test
metric_CoxNet_summary.to_csv("CoxNet0.5_maculopathy_risk_performance_metrics.txt", sep='\t',index=False)
