import pandas as pd
from sklearn.model_selection import train_test_split


df = pd.read_csv('../data/after_interpolate.csv')

unique_df = df[['加入者id','is_mac']].drop_duplicates()
X = unique_df['加入者id']
y = unique_df['is_mac']
train_ids, test_ids, y_train, y_test = train_test_split(X, y,test_size=0.2, random_state=42,stratify=y)

train_data = df[df['加入者id'].isin(train_ids)]
test_data = df[df['加入者id'].isin(test_ids)]

train_data.to_csv('../data/train_data.csv',index=False)
test_data.to_csv('../data/test_data.csv',index=False)

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

from sksurv.ensemble import RandomSurvivalForest
import _pickle as cPickle
from sklearn.inspection import permutation_importance
from sksurv.metrics import integrated_brier_score, cumulative_dynamic_auc
import logging
import os
import copy
import numpy as np

# ログの設定
logging.basicConfig(filename='log.txt', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ログメッセージの出力
logging.info('====================== start ======================')

n_rep = 10
path = '/work/kawakami/Project/Kimitsu_Diabetes/analysis'

# imp_summary = pd.read_csv(os.path.join(path, 'RSF_maculopathy_risk_varimp.csv'))
# imp_summary['median']=imp_summary.filter(like='rep').apply(lambda x: x.median(), axis=1)

cindex_summary = pd.DataFrame([])

# imp_features = imp_summary.sort_values('median', ascending=False).head(53)['Unnamed: 0'].values
# X_train_selected = X_train[imp_features]
# X_test_selected = X_test[imp_features]

imp_summary = pd.DataFrame([])
cindex_test = []
IBS_test = []
cAUC_test = []
array_psf = []

feature_names = X_train.columns.tolist()

for i in range(n_rep):
	print('RandomSurvivalForest start',i)
	logging.info("RandomSurvivalForest(n_estimators=200,min_samples_split=10,min_samples_leaf=15,max_features='sqrt',oob_score=True,n_jobs=31,random_state=i=%d)",i)

	rsf = RandomSurvivalForest(n_estimators=200,min_samples_split=10,min_samples_leaf=15,max_features="sqrt",oob_score=True,n_jobs=31,random_state=i)
	rsf.fit(X_train, yy_train)
	with open(os.path.join(f"{path}","RSF_models","RSF_maculopathy_risk_rep"+str(i+1)+".sav"), 'wb') as f:
		cPickle.dump(rsf, f)

	cindex_test.append(rsf.score(X_test, yy_test))
	survs = rsf.predict_survival_function(X_test)
	risk_scores = rsf.predict(X_test)
	times = np.arange(12, 120)
	preds = np.asarray([[fn(t) for t in times] for fn in survs])
	IBS_test.append(integrated_brier_score(yy_train, yy_test, preds, times))

	auc, mean_auc = cumulative_dynamic_auc(yy_train, yy_test, risk_scores, times)

	cAUC_test.append(mean_auc)

	# Survival probability
	logging.info('predict_survival_function')
	surv = rsf.predict_survival_function(X_test, return_array=True)
	df_psf = pd.DataFrame(surv,columns=rsf.unique_times_)
	df_psf = pd.concat([test_data[['加入者id','is_mac','month_maculopathy']].reset_index(drop=True),df_psf], axis=1)
	array_psf.append(df_psf)

	# # Permutation importance
	# logging.info('permutation_importance')
	# perm = permutation_importance(rsf,X_train, yy_train, n_repeats=5, random_state=i)

	# imp_summary = pd.concat([imp_summary, pd.DataFrame(perm.importances_mean)], axis=1)

# imp_summary.to_csv(os.path.join(f"{path}","RSF_maculopathy_risk_impsummary.csv"),index=False)

# df_psfs = pd.concat(array_psf)
# #df_psfs.to_csv(os.path.join(f"{path}","analyze","RSF_models","RSF_models_maculopathy_risk_permutation_importance_all.csv"),index=False)
# df_psfs_mean = df_psfs.groupby(['加入者id','month_maculopathy']).mean().reset_index()
# df_psfs_mean.to_csv(os.path.join(f"{path}","RSF_models","RSF_models_maculopathy_risk_survival_func_mean.csv"),index=False)

metric_summary = pd.DataFrame([])
metric_summary["cindex"] = cindex_test
metric_summary["IBS"] = IBS_test
metric_summary["cAUC"] = cAUC_test
metric_summary.to_csv("RSF_maculopathy_risk_performance_metrics.txt", sep='\t',index=False)

# imp_summary.columns = [ "rep_"+str(i+1) for i in range(n_rep)]
# imp_summary.index = X_train.columns
# imp_summary = imp_summary.reset_index()
# imp_summary.to_csv(os.path.join(f"{path}","RSF_maculopathy_risk_varimp.csv"),index=False)

# # 変数重要度
# imp_summary_tops = copy.deepcopy(imp_summary).set_index('index')
# imp_summary_tops["median"] = imp_summary_tops.median(axis="columns")
# imp_summary_tops = imp_summary_tops.sort_values("median",ascending=False)
# imp_summary_tops.to_csv(os.path.join(f"{path}","RSF_models","RSF_maculopathy_risk_varimp_descending.csv"))
