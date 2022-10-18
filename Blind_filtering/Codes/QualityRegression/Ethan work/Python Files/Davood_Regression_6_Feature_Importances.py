# Importing processing libraries
import os
import sys
import random
import mord
import numpy as np
import pandas as pd
from deslib.des import DESP
from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import make_pipeline
from joblib import parallel_backend
from sklearn import preprocessing
# Importing Pool Classifiers
from sklearn.ensemble import AdaBoostRegressor, GradientBoostingRegressor, BaggingRegressor, RandomForestRegressor
from sklearn.linear_model import Lars, LassoLars, LassoLarsCV, LarsCV, LinearRegression, \
    ElasticNet
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso, ElasticNetCV
from sklearn.metrics import accuracy_score, explained_variance_score, mean_squared_error
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVC, SVR
from sklearn.tree import DecisionTreeRegressor

from feature_selector import FeatureSelector
from ordinal_classifier import OrdinalClassifier
import pymrmr


# Importing DES and DCS techniques


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def save_results(single_prediction, actual_label, single_prediction_train, actual_train_label, save_result_out, Name):
    [single_prediction, acc, var_out, error_out, cm_out] = calculate_results(single_prediction, actual_label)
    [single_prediction_train, acc_train, var_train, error_train, cm_train] = calculate_results(single_prediction_train,
                                                                                               actual_train_label)

    save_result_out.append(
        [Name, acc, acc_train, var_out, var_train, error_out, error_train, cm_out, cm_train])

    return [save_result_out, single_prediction, single_prediction_train, acc, acc_train, var_out, var_train, error_out,
            error_train, cm_out, cm_train]


def calculate_results(single_prediction, actual_label):
    single_prediction = np.array(single_prediction)
    single_prediction[single_prediction < 1] = 1
    single_prediction[single_prediction > 5] = 5

    acc = accuracy_score(actual_label, np.round(single_prediction))
    var_out = explained_variance_score(actual_label, single_prediction)
    error_out = mean_squared_error(actual_label, single_prediction, squared=False)
    cm_out = confusion_matrix(actual_label, np.round(single_prediction))

    return [single_prediction, acc, var_out, error_out, cm_out]


def up_sample_method(neighbors1, neighbors2, neighbors3, neighbors4, neighbors5, regression_method):
    if min([neighbors1, neighbors2, neighbors3, neighbors4, neighbors5]) > 1:
        overall_pipeline = make_pipeline(
            SMOTE(random_state=42, k_neighbors=min([neighbors1, neighbors2, neighbors3, neighbors4, neighbors5])),
            regression_method)
    else:
        overall_pipeline = make_pipeline(RandomOverSampler(random_state=42), regression_method)
    return [overall_pipeline]


def fix_grouping(split_index, all_train_subjects, all_train_label):
    TrainSize = split_index[0].shape[0]
    TestSize = split_index[1].shape[0]
    for train_sub in set(all_train_subjects):
        if all_train_subjects[split_index[0]].__contains__(train_sub) & all_train_subjects[split_index[1]].__contains__(
                train_sub):
            train_loc = all_train_subjects[split_index[0]] == train_sub
            test_loc = all_train_subjects[split_index[1]] == train_sub

            # remove from training and add to test
            train_alt_1 = split_index[0][~train_loc]
            test_alt_1 = np.append(split_index[1], split_index[0][train_loc])

            t = all_train_label[train_alt_1]
            train_alt_1_dist = np.array([sum(t == 1), sum(t == 2), sum(t == 3), sum(t == 4), sum(t == 5)]) / t.shape[0]
            t = all_train_label[test_alt_1]
            test_alt_1_dist = np.array([sum(t == 1), sum(t == 2), sum(t == 3), sum(t == 4), sum(t == 5)]) / t.shape[0]

            # remove from test and add to training
            train_alt_2 = np.append(split_index[0], split_index[1][test_loc])
            test_alt_2 = split_index[1][~test_loc]

            t = all_train_label[train_alt_2]
            train_alt_2_dist = np.array([sum(t == 1), sum(t == 2), sum(t == 3), sum(t == 4), sum(t == 5)]) / t.shape[0]
            t = all_train_label[test_alt_2]
            test_alt_2_dist = np.array([sum(t == 1), sum(t == 2), sum(t == 3), sum(t == 4), sum(t == 5)]) / t.shape[
                0]

            opt1 = np.linalg.norm(np.append(train_alt_1_dist - test_alt_1_dist,
                                            [(train_alt_1.shape[0] - TrainSize) / TrainSize,
                                             (test_alt_1.shape[0] - TestSize) / TestSize]))
            opt2 = np.linalg.norm(np.append(train_alt_2_dist - test_alt_2_dist,
                                            [(train_alt_2.shape[0] - TrainSize) / TrainSize,
                                             (test_alt_2.shape[0] - TestSize) / TestSize]))
            if opt1 > opt2:
                split_index[0] = train_alt_2
                split_index[1] = test_alt_2
            else:
                split_index[0] = train_alt_1
                split_index[1] = test_alt_1
    return split_index


# file_path = 'Signal Quality.xlsx'
file_path = 'Regression Extracted Features.xlsx'
# file_path = '/Users/ethangrooby/PycharmProjects/DES_1/Complete Set.xlsx'
# annotation file
# Sheets include:
# Lung Output, Lung Features, Lung Subject Grouping
# Heart Output, Heart Features, Heart Subject Grouping
target = 'Heart'
sheet = pd.read_excel(io=file_path, sheet_name=None)
#remove_files = pd.read_excel(io=file_path, sheet_name=target + ' Sync')
tgtFea = pd.read_excel(io=file_path, sheet_name=target + ' Features')
#tgtFea = tgtFea.loc[~remove_files['Location']]
Sub = pd.read_excel(io=file_path, sheet_name=target + ' Subjects')
#Sub = Sub.loc[~remove_files['Location']]
unique_Sub = set(Sub.Subject)
tgtLab = pd.read_excel(io=file_path, sheet_name=target + ' Quality')
#tgtLab = tgtLab.loc[~remove_files['Location']]
# binary label processing
label = np.array(tgtLab['Quality'])
label = label[~np.isnan(label)]



'''
slow_features= [
            "'trunc_autocorr_se1'",
            "'trunc_autocorr_se2'",
            "'trunc_autocorr_se3'",
            "'trunc_autocorr_se4'",
            "'trunc_autocorr_se5'",
            "'strunc_autocorr_se1'",
            "'strunc_autocorr_se2'",
            "'strunc_autocorr_se3'",
            "'strunc_autocorr_se4'",
            "'strunc_autocorr_se5'",
            "'audio_rate_average_energy_1000'",
            "'audio_rate_average_energy_2000'",
            "'trunc_autocorr_svd_heart'",
            "'strunc_autocorr_svd_heart'",
            "'trunc_autocorr_svd_lung'",
            "'strunc_autocorr_svd_lung'",
            "'percentage_bad_cycles_sys_schmidt'",
            "'percentage_bad_cycles_dia_schmidt'",
            "'percentage_bad_cycles_s1_schmidt'",
            "'percentage_bad_cycles_s2_schmidt'",
            "'percentage_bad_cycles_overall_schmidt'",
            "'seg_qual_schmidt'",
            "'seg_qual_avg_schmidt'",
            "'seg_qual_inv_schmidt'",
            "'seg_qual_inv_avg_schmidt'",
            "'seg_agreement1'",
            "'seg_agreement3'",
            "'seg_hsmm_qf_s1_psd_schmidt'",
            "'seg_hsmm_qf_s2_psd_schmidt'",
            "'seg_hsmm_qf_s1_psd_ori_schmidt'",
            "'seg_hsmm_qf_s2_psd_ori_schmidt'",
            "'seg_hsmm_qf_s1_wav_schmidt'",
            "'seg_hsmm_qf_s2_wav_schmidt'",
            "'seg_hsmm_qf_s1_homo_schmidt'",
            "'seg_hsmm_qf_s2_homo_schmidt'",
            "'seg_hsmm_qf_s1_hil_schmidt'",
            "'seg_hsmm_qf_s2_hil_schmidt'",
            "'seg_hsmm_qf_s1_stft_schmidt'",
            "'seg_hsmm_qf_s2_stft_schmidt'",
            "'seg_hsmm_qf_s1_shannon_schmidt'",
            "'seg_hsmm_qf_s2_stft_schmidt'",
            "'seg_per_abnormal_rmssd_schmidt'",
            "'seg_per_abnormal_zero_schmidt'",
            "'seg_per_abnormal_sd1_schmidt'",
            "'wav_schmidt_S1S2_per_peak_window_original'",
            "'wav_schmidt_S1S2_per_peak_window_modified'",
            "'wav_schmidt_S1S2_per_peak_window_modified2'",
            "'stft_envelope_var'",
            "'stft_envelope_se'",
            "'stft_envelope_cycle_corr_avg'",
            "'stft_envelope_cycle_corr_std'",
            "'stft_envelope_hr_std'",
            "'seg_hsmm_qf_s1_stft'",
            "'seg_hsmm_qf_s2_stft'"
]

tgtFea = tgtFea.drop(slow_features, axis=1)
'''

save_result = []
arr = np.array([])
rng = np.random.RandomState(1)

# Leave-One-Out
LOO_label = []
LOO_label_train_pool = []
test_prediction = []
train_prediction = []
cv_results = []
best_num_features = []
overall_best_model_list = []
with parallel_backend('threading', n_jobs=8):
    # getting test and train set
    train_label = label
    train_data = tgtFea
    train_subjects = np.array(Sub.Subject)

    # normalization
    scaler = preprocessing.StandardScaler().fit(train_data)
    train_data = pd.DataFrame(scaler.transform(train_data), columns=train_data.columns)

    [train_data_upsampled_all, train_label_upsampled_all] = RandomOverSampler(random_state=42).fit_sample(
        train_data, train_label)

    train_data_upsampled_all.insert(0, 'quality', train_label_upsampled_all)

    important_feature_list = pymrmr.mRMR(train_data_upsampled_all, 'MID', 20)
    feature_loc = []
    for i in range(0, 20):
        feature_loc = np.append(feature_loc, train_data.columns.get_loc(important_feature_list[i]))

    k = 5
    quality_values = train_label
    subject_groupings = train_subjects
    subject_grouping_constant = train_subjects
    n = set(quality_values)
    n_len = n.__len__()
    splits_count = np.zeros([k, n_len])
    split_allocation = []
    for i in range(k):
        split_allocation.append(np.array([]))
    for idx, val in enumerate(n):
        shuffled_sub = list(set(subject_groupings[quality_values == val]))
        random.shuffle(shuffled_sub)
        for sub in shuffled_sub:
            fold_loc = np.argmin(splits_count[:, idx])
            quality_values_assign = quality_values[subject_groupings == sub]
            for idx2, val2 in enumerate(n):
                splits_count[fold_loc, idx2] += sum(quality_values_assign == val2)
            split_allocation[fold_loc] = np.append(split_allocation[fold_loc],
                                                   np.where(subject_grouping_constant == sub)[0])
            quality_values = np.delete(quality_values, subject_groupings == sub)
            subject_groupings = np.delete(subject_groupings, subject_groupings == sub)
    for i in range(k):
        random.shuffle(split_allocation[i])
        split_allocation[i] = split_allocation[i].astype(int)

    s1 = (split_allocation[0],
          np.concatenate((split_allocation[1], split_allocation[2], split_allocation[3], split_allocation[4])))
    s2 = (split_allocation[1],
          np.concatenate((split_allocation[0], split_allocation[2], split_allocation[3], split_allocation[4])))
    s3 = (split_allocation[2],
          np.concatenate((split_allocation[1], split_allocation[0], split_allocation[3], split_allocation[4])))
    s4 = (split_allocation[3],
          np.concatenate((split_allocation[1], split_allocation[2], split_allocation[0], split_allocation[4])))
    s5 = (split_allocation[4],
          np.concatenate((split_allocation[1], split_allocation[2], split_allocation[3], split_allocation[0])))

    # splits_count
    # np.sum(splits_count, axis=1)
    # train_label[split_allocation[0].astype(int)]
    # sum(train_label[split_allocation[0].astype(int)]==1)
    # set(train_subjects[split_allocation[0].astype(int)])

    # sorting out cross validation
    train_data = np.array(train_data)
    train_label = np.array(train_label)

    # Up-sample training set
    [train_data_s1, train_label_s1] = RandomOverSampler(random_state=42).fit_sample(train_data[s1[0]],
                                                                                    train_label[s1[0]])
    [train_data_s2, train_label_s2] = RandomOverSampler(random_state=42).fit_sample(train_data[s2[0]],
                                                                                    train_label[s2[0]])
    [train_data_s3, train_label_s3] = RandomOverSampler(random_state=42).fit_sample(train_data[s3[0]],
                                                                                    train_label[s3[0]])
    [train_data_s4, train_label_s4] = RandomOverSampler(random_state=42).fit_sample(train_data[s4[0]],
                                                                                    train_label[s4[0]])
    [train_data_s5, train_label_s5] = RandomOverSampler(random_state=42).fit_sample(train_data[s5[0]],
                                                                                    train_label[s5[0]])

    [test_data_s1, test_label_s1] = RandomOverSampler(random_state=42).fit_sample(train_data[s1[1]],
                                                                                  train_label[s1[1]])
    [test_data_s2, test_label_s2] = RandomOverSampler(random_state=42).fit_sample(train_data[s2[1]],
                                                                                  train_label[s2[1]])
    [test_data_s3, test_label_s3] = RandomOverSampler(random_state=42).fit_sample(train_data[s3[1]],
                                                                                  train_label[s3[1]])
    [test_data_s4, test_label_s4] = RandomOverSampler(random_state=42).fit_sample(train_data[s4[1]],
                                                                                  train_label[s4[1]])
    [test_data_s5, test_label_s5] = RandomOverSampler(random_state=42).fit_sample(train_data[s5[1]],
                                                                                  train_label[s5[1]])

    [train_data_sall, train_label_sall] = RandomOverSampler(random_state=42).fit_sample(train_data,
                                                                                        train_label)

    # merge training and then associated test set
    train_data_upsampled2 = np.append(train_data_s1, test_data_s1, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, train_data_s2, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, test_data_s2, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, train_data_s3, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, test_data_s3, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, train_data_s4, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, test_data_s4, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, train_data_s5, axis=0)
    train_data_upsampled2 = np.append(train_data_upsampled2, test_data_s5, axis=0)

    # same for labels
    train_label_upsampled2 = np.append(train_label_s1, test_label_s1)
    train_label_upsampled2 = np.append(train_label_upsampled2, train_label_s2)
    train_label_upsampled2 = np.append(train_label_upsampled2, test_label_s2)
    train_label_upsampled2 = np.append(train_label_upsampled2, train_label_s3)
    train_label_upsampled2 = np.append(train_label_upsampled2, test_label_s3)
    train_label_upsampled2 = np.append(train_label_upsampled2, train_label_s4)
    train_label_upsampled2 = np.append(train_label_upsampled2, test_label_s4)
    train_label_upsampled2 = np.append(train_label_upsampled2, train_label_s5)
    train_label_upsampled2 = np.append(train_label_upsampled2, test_label_s5)

    # create split indices
    s1tr = len(train_label_s1)
    s2tr = len(train_label_s2)
    s3tr = len(train_label_s3)
    s4tr = len(train_label_s4)
    s5tr = len(train_label_s5)

    s1te = len(test_label_s1)
    s2te = len(test_label_s2)
    s3te = len(test_label_s3)
    s4te = len(test_label_s4)
    s5te = len(test_label_s5)

    s1a2 = (np.arange(s1tr), np.arange(s1tr, s1tr + s1te))
    total = s1tr + s1te
    s2a2 = (np.arange(total, s2tr + total), np.arange(s2tr + total, s2tr + s2te + total))
    total = total + s2tr + s2te
    s3a2 = (np.arange(total, s3tr + total), np.arange(s3tr + total, s3tr + s3te + total))
    total = total + s3tr + s3te
    s4a2 = (np.arange(total, s4tr + total), np.arange(s4tr + total, s4tr + s4te + total))
    total = total + s4tr + s4te
    s5a2 = (np.arange(total, s5tr + total), np.arange(s5tr + total, s5tr + s5te + total))

    s1a = s1a2
    s2a = s2a2
    s3a = s3a2
    s4a = s4a2
    s5a = s5a2
    train_data_upsampled = train_data_upsampled2
    train_label_upsampled = train_label_upsampled2

    train_data_upsampled_full = train_data_upsampled
    train_data_full = train_data
    train_data_s1_full = train_data_s1
    train_data_s2_full = train_data_s2
    train_data_s3_full = train_data_s3
    train_data_s4_full = train_data_s4
    train_data_s5_full = train_data_s5
    test_data_s1_full = test_data_s1
    test_data_s2_full = test_data_s2
    test_data_s3_full = test_data_s3
    test_data_s4_full = test_data_s4
    test_data_s5_full = test_data_s5
    train_data_sall_full = train_data_sall
    models = []
    model_results = []
    best_models = []
    best_model_results = []
    for num_features in range(5, 16):
        print(num_features)

        train_data_upsampled = train_data_upsampled_full[:, feature_loc.astype(int)[:num_features]]
        train_data = train_data_full[:, feature_loc.astype(int)[:num_features]]
        train_data_s1 = train_data_s1_full[:, feature_loc.astype(int)[:num_features]]
        train_data_s2 = train_data_s2_full[:, feature_loc.astype(int)[:num_features]]
        train_data_s3 = train_data_s3_full[:, feature_loc.astype(int)[:num_features]]
        train_data_s4 = train_data_s4_full[:, feature_loc.astype(int)[:num_features]]
        train_data_s5 = train_data_s5_full[:, feature_loc.astype(int)[:num_features]]
        test_data_s1 = test_data_s1_full[:, feature_loc.astype(int)[:num_features]]
        test_data_s2 = test_data_s2_full[:, feature_loc.astype(int)[:num_features]]
        test_data_s3 = test_data_s3_full[:, feature_loc.astype(int)[:num_features]]
        test_data_s4 = test_data_s4_full[:, feature_loc.astype(int)[:num_features]]
        test_data_s5 = test_data_s5_full[:, feature_loc.astype(int)[:num_features]]
        train_data_sall = train_data_sall_full[:, feature_loc.astype(int)[:num_features]]

        # Alpha parameter allows dealing with large number of features through a penality with size of coefficients
        reg_ridge = RidgeCV(normalize=False, alphas=[0.1, 0.5, 1.0, 5, 10.0, 50, 100, 500, 1000],
                            cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        reg_ridge.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(-reg_ridge.best_score_)
        reg_ridge = Ridge(normalize=False, alpha=reg_ridge.alpha_)
        reg_ridge.fit(train_data_sall, train_label_sall)
        models.append(reg_ridge)

        # Good with large feature sets
        # reg_lars = LarsCV(normalize=False, max_iter=50, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_lars = LarsCV(normalize=False, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_lars.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(np.min(np.mean(reg_lars.mse_path_, axis=1)))
        reg_lars = Lars(normalize=False)
        reg_lars.fit(train_data_sall, train_label_sall)
        models.append(reg_lars)

        # The Lasso is a linear model that estimates sparse coefficients as controlled by alpha. Therefore good
        # with dealing with sets with large number of features. reg_lasso = LassoCV(normalize=False,
        # max_iter=20000, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_lasso = LassoCV(normalize=False, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_lasso.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(np.min(np.mean(reg_lasso.mse_path_, axis=1)))
        reg_lasso = Lasso(normalize=False, alpha=reg_lasso.alpha_)
        reg_lasso.fit(train_data_sall, train_label_sall)
        models.append(reg_lasso)

        # Similar to lasso but with different algorithm to fit parameters
        # reg_lassolars = LassoLarsCV(normalize=False, max_iter=50, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_lassolars = LassoLarsCV(normalize=False, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_lassolars.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(np.min(np.mean(reg_lassolars.mse_path_, axis=1)))
        reg_lassolars = LassoLars(normalize=False, alpha=reg_lassolars.alpha_)
        reg_lassolars.fit(train_data_sall, train_label_sall)
        models.append(reg_lassolars)

        # Combines the sparsity of lasso with the regularisation of ridge
        # reg_elastic = ElasticNetCV(normalize=False, max_iter=20000, cv=[s1a, s2a, s3a, s4a, s5a])
        reg_elastic = ElasticNetCV(normalize=False,
                                   l1_ratio=[0.001, 0.005, 0.01, 0.05, .1, .5, .7, .9, .95, .99, 0.995, 0.999, 1],
                                   cv=[s1a, s2a, s3a, s4a, s5a])
        reg_elastic.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(np.min(np.mean(reg_elastic.mse_path_, axis=2)))
        reg_elastic = ElasticNet(normalize=False, alpha=reg_elastic.alpha_, l1_ratio=reg_elastic.l1_ratio_)
        reg_elastic.fit(train_data_sall, train_label_sall)
        models.append(reg_elastic)

        # Requires feature selection as it relies on the independence of features
        scores = cross_val_score(LinearRegression(normalize=False), train_data_upsampled, train_label_upsampled,
                                 cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        reg_linear = LinearRegression(normalize=False).fit(train_data_sall, train_label_sall)
        # model_results.append(np.mean(reg_linear_error))
        model_results.append(-np.mean(scores))
        models.append(reg_linear)

        # Imagine would require feature selection
        param_grid = {'n_neighbors': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 'weights': ['uniform', 'distance'],
                      'algorithm': ['ball_tree', 'kd_tree', 'brute'], 'p': [1, 2]}
        new_params = {'kneighborsregressor__' + key: param_grid[key] for key in param_grid}
        calibrated_knn = KNeighborsRegressor()
        imba_pipeline = make_pipeline(calibrated_knn)
        search_knn = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                  scoring='neg_mean_squared_error',
                                  return_train_score=True)
        search_knn.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(-search_knn.best_score_)
        search_knn = search_knn.best_estimator_.named_steps['kneighborsregressor']
        search_knn.fit(train_data_sall, train_label_sall)
        models.append(search_knn)

        # Imagine would require feature selection
        param_grid = {'criterion': ['mse', 'friedman_mse', 'mae', 'poisson'],
                      'max_depth': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, None],
                      'max_features': ['auto', 'sqrt', 'log2']}
        new_params = {'decisiontreeregressor__' + key: param_grid[key] for key in param_grid}
        calibrated_tree = DecisionTreeRegressor(random_state=rng)
        imba_pipeline = make_pipeline(calibrated_tree)
        search_tree = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                   scoring='neg_mean_squared_error',
                                   return_train_score=True)
        search_tree.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(-search_tree.best_score_)
        search_tree = search_tree.best_estimator_.named_steps['decisiontreeregressor']
        search_tree.fit(train_data_sall, train_label_sall)
        models.append(search_tree)

        # Has regularisation
        new_params = [
            {'svr__C': [0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12], 'svr__kernel': ['linear']},
            {'svr__C': [0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12],
             'svr__gamma': ['scale', 'auto', 0.1, 0.01, 0.001,
                            0.0001], 'svr__kernel': ['rbf']},
        ]
        calibrated_svr = SVR(gamma='auto')
        imba_pipeline = make_pipeline(calibrated_svr)
        search_svr = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                  scoring='neg_mean_squared_error',
                                  return_train_score=True)
        search_svr.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(-search_svr.best_score_)
        search_svr = search_svr.best_estimator_.named_steps['svr']
        search_svr.fit(train_data_sall, train_label_sall)
        models.append(search_svr)

        # Has regularisation
        new_params = {'logisticat__alpha': [0, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000]}
        calibrated_logisticat = mord.LogisticAT(alpha=1.0, verbose=0, max_iter=10000)
        imba_pipeline = make_pipeline(calibrated_logisticat)
        search_logisticat = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                         scoring='neg_mean_squared_error',
                                         return_train_score=True)
        search_logisticat.fit(train_data_upsampled, train_label_upsampled.astype(int))
        model_results.append(-search_logisticat.best_score_)
        search_logisticat = search_logisticat.best_estimator_
        search_logisticat.fit(train_data_sall, train_label_sall)
        models.append(search_logisticat)

        # Has regularisation
        new_params = {'logisticit__alpha': [0, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000]}
        calibrated_logisticit = mord.LogisticIT(alpha=1.0, verbose=0, max_iter=10000)
        imba_pipeline = make_pipeline(calibrated_logisticit)
        search_logisticit = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                         scoring='neg_mean_squared_error',
                                         return_train_score=True)
        search_logisticit.fit(train_data_upsampled, train_label_upsampled.astype(int))
        model_results.append(-search_logisticit.best_score_)
        search_logisticit = search_logisticit.best_estimator_
        search_logisticit.fit(train_data_sall, train_label_sall)
        models.append(search_logisticit)

        # Has regularisation
        new_params = {'ordinalridge__alpha': [0, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000]}
        calibrated_oridge = mord.OrdinalRidge()
        imba_pipeline = make_pipeline(calibrated_oridge)
        search_oridge = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                     scoring='neg_mean_squared_error',
                                     return_train_score=True)
        search_oridge.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(-search_oridge.best_score_)
        search_oridge = search_oridge.best_estimator_
        search_oridge.fit(train_data_sall, train_label_sall)
        models.append(search_oridge)

        # Probably requires feature selection
        scores = cross_val_score(mord.LAD(max_iter=5000), train_data_upsampled, train_label_upsampled,
                                 cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        ord_lad = mord.LAD(max_iter=5000).fit(train_data_sall, train_label_sall)
        # model_results.append(np.mean(ord_lad_error))
        model_results.append(-np.mean(scores))
        models.append(ord_lad)

        # Takes too long 2
        # Has regularisation
        new_params = [
            {'ordinalclassifier__C': [0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12],
             'ordinalclassifier__kernel': ['linear']},
            {'ordinalclassifier__C': [0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12],
             'ordinalclassifier__gamma': ['scale', 'auto', 0.1, 0.01, 0.001,
                                          0.0001], 'ordinalclassifier__kernel': ['rbf']},
        ]
        ord_svc = OrdinalClassifier(SVC(probability=True))
        imba_pipeline = make_pipeline(ord_svc)
        search_osvc = GridSearchCV(imba_pipeline, param_grid=new_params, cv=[s1a, s2a, s3a, s4a, s5a],
                                   scoring='neg_mean_squared_error',
                                   return_train_score=True)
        search_osvc.fit(train_data_upsampled, train_label_upsampled)
        model_results.append(-search_osvc.best_score_)
        search_osvc = search_osvc.best_estimator_.named_steps['ordinalclassifier']
        search_osvc.fit(train_data_sall, train_label_sall)
        models.append(search_osvc)

        # Ensemble
        scores = cross_val_score(AdaBoostRegressor(loss='square'), train_data_upsampled, train_label_upsampled,
                                 cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        ensemble_ada = AdaBoostRegressor(loss='square').fit(train_data_sall, train_label_sall)
        model_results.append(-np.mean(scores))
        models.append(ensemble_ada)

        scores = cross_val_score(GradientBoostingRegressor(), train_data_upsampled, train_label_upsampled,
                                 cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        ensemble_grad = GradientBoostingRegressor().fit(train_data_sall, train_label_sall)
        model_results.append(-np.mean(scores))
        models.append(ensemble_grad)

        scores = cross_val_score(BaggingRegressor(), train_data_upsampled, train_label_upsampled,
                                 cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        ensemble_bag = BaggingRegressor().fit(train_data_sall, train_label_sall)
        model_results.append(-np.mean(scores))
        models.append(ensemble_bag)

        scores = cross_val_score(RandomForestRegressor(), train_data_upsampled, train_label_upsampled,
                                 cv=[s1a, s2a, s3a, s4a, s5a], scoring='neg_mean_squared_error')
        ensemble_forest = RandomForestRegressor().fit(train_data_sall, train_label_sall)
        model_results.append(-np.mean(scores))
        models.append(ensemble_forest)

        best_models.append(models[np.argmin(model_results)])
        best_model_results.append(np.min(model_results))
        models = []
        model_results = []

    overall_best_model = best_models[np.argmin(best_model_results)]
    num_features = np.argmin(best_model_results) + 5  # remember to change this again
    train_data = train_data_full[:, feature_loc.astype(int)[:num_features]]
    train_prediction += list(overall_best_model.predict(train_data))
    LOO_label_train_pool += list(train_label)
    cv_results += list(np.array([np.min(best_model_results)]))
    best_num_features += list(np.array([num_features]))
    overall_best_model_list.append(overall_best_model)



import pickle
with open('Regression Heart.pickle', 'wb') as f:
    pickle.dump([scaler, overall_best_model, important_feature_list,num_features], f)






