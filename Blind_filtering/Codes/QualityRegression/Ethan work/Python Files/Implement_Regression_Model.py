import pickle
import numpy as np
import pandas as pd
from numpy import savetxt

with open('Regression Lung.pickle', 'rb') as f:
    scaler, overall_best_model, important_feature_list, num_features = pickle.load(f)

    file_path = 'Davood Results.xlsx'
    target = 'Lung'
    sheet = pd.read_excel(io=file_path, sheet_name=None)
    tgtFea = pd.read_excel(io=file_path, sheet_name=target + ' Features')

    tgtFea = pd.DataFrame(scaler.transform(tgtFea), columns=tgtFea.columns)
    feature_loc = []
    for i in range(0, 20):
        feature_loc = np.append(feature_loc, tgtFea.columns.get_loc(important_feature_list[i]))

    tgtFea=np.array(tgtFea)
    tgtFea = tgtFea[:, feature_loc.astype(int)[:num_features]]

    savetxt('lung_regression_results.csv',overall_best_model.predict(tgtFea), delimiter=',')
    print(important_feature_list)

with open('Regression Heart.pickle', 'rb') as f:
    scaler, overall_best_model, important_feature_list, num_features = pickle.load(f)

    file_path = 'Davood Results.xlsx'
    target = 'Heart'
    sheet = pd.read_excel(io=file_path, sheet_name=None)
    tgtFea = pd.read_excel(io=file_path, sheet_name=target + ' Features')

    tgtFea = pd.DataFrame(scaler.transform(tgtFea), columns=tgtFea.columns)
    feature_loc = []
    for i in range(0, 20):
        feature_loc = np.append(feature_loc, tgtFea.columns.get_loc(important_feature_list[i]))

    tgtFea=np.array(tgtFea)
    tgtFea = tgtFea[:, feature_loc.astype(int)[:num_features]]

    savetxt('heart_regression_results.csv',overall_best_model.predict(tgtFea), delimiter=',')
    print(important_feature_list)
