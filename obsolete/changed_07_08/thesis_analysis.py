import argparse
import os
import random
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors, LocalOutlierFactor
from sklearn.cluster import OPTICS
from sklearn.ensemble import IsolationForest
import hdbscan

default_distance_metric = 'braycurtis'
outlier_value = 0.15

def set_initial_results_dataframe(csv_file_path):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    results_df = pd.DataFrame(index=temp_df.index)
    results_df.insert(0, 'id', temp_df.index)

    return results_df

def set_anomaly_by_id(row, outlier_list):
    if row.id in outlier_list:
        return 'anomaly'
    return 'notanomaly'

def set_anomaly_by_score(row):
    if row['anomaly_score'] == -1:
        return 'anomaly'
    return 'notanomaly'


# n_neighbors = default
# default distance metric = braycurtis
# outlier_value = 0.15
def run_LOF(csv_file_path, results_df, distance=default_distance_metric):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance == 'braycurtis' else 'precomputed'

    model = LocalOutlierFactor(metric=diststr, contamination=outlier_value, n_neighbors=len(temp_df))
    temp_df['anomaly_score'] = model.fit_predict(temp_df.values)
    temp_df['anomaly_score2'] = temp_df.apply(set_anomaly_by_score, axis=1)

    column_name = 'LOF_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = temp_df['anomaly_score2']

    return results_df

# n_neighbors = default
# default distance metric = braycurtis
# outlier_value = 0.15
def run_KNN(csv_file_path, results_df, distance=default_distance_metric):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance == 'braycurtis' else 'precomputed'

    model = NearestNeighbors(metric=diststr)
    model.fit(temp_df.values)
    distances, indexes = model.kneighbors(temp_df.values)

    threshold = pd.Series(distances.mean(axis = 1)).quantile(1 - outlier_value)
    outlier_index = np.where(distances.mean(axis = 1) > threshold)
    outlier_values = temp_df.iloc[outlier_index]
    outlier_samples = outlier_values.index.values.tolist()

    column_name = 'KNN_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = results_df.apply(set_anomaly_by_id, outlier_list=outlier_samples , axis=1)

    return results_df


# default distance metric = braycurtis
def run_OPTICS(csv_file_path, results_df, distance=default_distance_metric):
    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance == 'braycurtis' else 'precomputed'

    model = OPTICS(metric=diststr)
    model.fit(temp_df.values)

    distances = model.core_distances_
    threshold = np.quantile(distances, 1 - outlier_value)
    outlier_index = np.where(distances >= threshold)
    outlier_values = temp_df.iloc[outlier_index]
    outlier_samples = outlier_values.index.values.tolist()

    column_name = 'OPTICS_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = results_df.apply(set_anomaly_by_id, outlier_list=outlier_samples, axis=1)

    return results_df


# cluster_selection_method = leaf
# default distance metric = braycurtis
def run_HDBScan(csv_file_path, results_df, distance=default_distance_metric):
    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance == 'braycurtis' else 'precomputed'

    model = hdbscan.HDBSCAN(metric=diststr, cluster_selection_method='leaf', allow_single_cluster=True)
    model.fit(temp_df.values)

    threshold = pd.Series(model.outlier_scores_).quantile(1 - outlier_value)
    outliers = np.where(model.outlier_scores_ > threshold)[0]
    outlier_values = temp_df.iloc[outliers]
    outlier_samples = outlier_values.index.values.tolist()

    column_name = 'HDBSCAN_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = results_df.apply(set_anomaly_by_id, outlier_list=outlier_samples, axis=1)

    return results_df

# number estimators = 200
# max samples = auto
def run_IsolationForest(csv_file_path, results_df):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    model = IsolationForest(n_estimators=200, contamination=outlier_value)
    temp_df['anomaly_score'] = model.fit_predict(temp_df.values)
    temp_df['anomaly_score2'] = temp_df.apply(set_anomaly_by_score, axis=1)

    column_name = 'ISOLATIONF_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3])
    results_df[column_name] = temp_df['anomaly_score2']

    return results_df

def Main():

    parser = argparse.ArgumentParser(description='Run ML script to generate outlier classification')
    parser.add_argument("-s", "--study_name", help="STUDY NAME  (required)", required=True)
    parser.add_argument("-r", "--taxa_rank", help="TAXA RANK AGGREGATION (phylum, genus)", required=True)
    parser.add_argument("-a", "--transformation", help="OTU ABUNDANCE TRANSFORMATION (relative, clr , log10)", required=True)
    parser.add_argument("-d", "--distance_measure", help="DISTANCE USED IN R SCRIPT (gunifrac, unifrac, wunifrac)", required=True)
    args = parser.parse_args()

    random.seed(39)

    study_name = args.study_name
    taxa_rank = args.taxa_rank
    transformation = args.transformation
    distance_measure = args.distance_measure

    class1_otu = class2_otu = class1_distances = class2_distances = None

    main_run_name = f'{taxa_rank}_{transformation}_{distance_measure}'
    main_temp_folder = os.path.join('results', study_name, main_run_name, 'tempfiles')
    main_results_folder = os.path.join('results', study_name, main_run_name, 'results')

    tempfiles = os.listdir(main_temp_folder)
    for file in tempfiles:
        if all(x in file for x in [study_name, taxa_rank, transformation, distance_measure, 'class1', 'otus']):
            class1_otu = os.path.join(main_temp_folder,file)
        elif all(x in file for x in [study_name, taxa_rank, transformation, distance_measure, 'class1', 'distances']):
            class1_distances = os.path.join(main_temp_folder,file)
        elif all(x in file for x in [study_name, taxa_rank, transformation, distance_measure, 'class2', 'otus']):
            class2_otu = os.path.join(main_temp_folder,file)
        elif all(x in file for x in [study_name, taxa_rank, transformation, distance_measure, 'class2', 'distances']):
            class2_distances = os.path.join(main_temp_folder,file)

    class1_df_results = set_initial_results_dataframe(class1_otu)
    class2_df_results = set_initial_results_dataframe(class2_otu)

    # LOF
    print('RUNNING LOCAL OUTLIER FACTOR MODEL ON CLASSES 1 AND 2')
    class1_df_results = run_LOF(class1_otu, class1_df_results)
    class1_df_results = run_LOF(class1_distances, class1_df_results, distance=distance_measure)
    class2_df_results = run_LOF(class2_otu, class2_df_results)
    class2_df_results = run_LOF(class2_distances, class2_df_results, distance=distance_measure)

    # KNN
    print('RUNNING K-NEAREST NEIGHBORS MODEL ON CLASSES 1 AND 2')
    class1_df_results = run_KNN(class1_otu, class1_df_results)
    class1_df_results = run_KNN(class1_distances, class1_df_results, distance=distance_measure)
    class2_df_results = run_KNN(class2_otu, class2_df_results)
    class2_df_results = run_KNN(class2_distances, class2_df_results, distance=distance_measure)

    # OPTICS
    print('RUNNING OPTICS MODEL ON CLASS 1 AND 2')
    class1_df_results = run_OPTICS(class1_otu, class1_df_results)
    class1_df_results = run_OPTICS(class1_distances, class1_df_results, distance=distance_measure)
    class2_df_results = run_OPTICS(class2_otu, class2_df_results)
    class2_df_results = run_OPTICS(class2_distances, class2_df_results, distance=distance_measure)

    # HDBSCAN
    print('RUNNING HDBSCAN MODEL ON CLASS 1 AND 2')
    class1_df_results = run_HDBScan(class1_otu, class1_df_results)
    class1_df_results = run_HDBScan(class1_distances, class1_df_results, distance=distance_measure)
    class2_df_results = run_HDBScan(class2_otu, class2_df_results)
    class2_df_results = run_HDBScan(class2_distances, class2_df_results, distance=distance_measure)

    # Isolation Forest
    print('RUNNING ISOLATION FOREST MODEL ON CLASS 1 AND 2')
    class1_df_results = run_IsolationForest(class1_otu, class1_df_results)
    class2_df_results = run_IsolationForest(class2_otu, class2_df_results)

    print('EXPORTING RESULTS')
    class1_df_results.to_csv(os.path.join(main_results_folder, 'anomaly_class1.csv'), index=False)
    class2_df_results.to_csv(os.path.join(main_results_folder, 'anomaly_class2.csv'), index=False)


if __name__ == '__main__':
    Main()