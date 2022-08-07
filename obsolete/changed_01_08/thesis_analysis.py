import argparse
import os
import random
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors, LocalOutlierFactor
from sklearn.cluster import DBSCAN
from sklearn.ensemble import IsolationForest
import hdbscan

## SET RANDOM STATE
random.seed(39)

unifrac_distances = ['unifrac', 'wunifrac', 'gunifrac', 'va-wunifrac']

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


def run_LOF(nneighbours, distance, outlier_value, csv_file_path, results_df):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance not in unifrac_distances else 'precomputed'

    model = LocalOutlierFactor(n_neighbors=nneighbours, metric=diststr, contamination=float(outlier_value))
    temp_df['anomaly_score'] = model.fit_predict(temp_df.values)
    temp_df['anomaly_score2'] = temp_df.apply(set_anomaly_by_score, axis=1)

    column_name = 'lof_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = temp_df['anomaly_score2']

    return results_df

# n_neighbors = default
# distance metric = braycurtis
# outlier_value = 0.15


def run_KNN(nneighbours, distance, outlier_value, csv_file_path, results_df):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance not in unifrac_distances else 'precomputed'

    model = NearestNeighbors(n_neighbors=nneighbours, metric=diststr)
    model.fit(temp_df.values)
    distances, indexes = model.kneighbors(temp_df.values)

    threshold = pd.Series(distances[:, 1]).quantile(1-outlier_value)
    outlier_index = np.where(distances[:, 1] > threshold)
    outlier_values = temp_df.iloc[outlier_index]
    outlier_samples = outlier_values.index.values.tolist()

    column_name = 'knn_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = results_df.apply(set_anomaly_by_id, outlier_list=outlier_samples , axis=1)

    return results_df

def run_IsolationForest(outlier_value, csv_file_path, results_df):

    random_state = np.random.RandomState(39)
    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    model = IsolationForest(n_estimators=100, max_samples='auto', contamination=outlier_value, random_state=random_state)
    temp_df['anomaly_score'] = model.fit_predict(temp_df.values)
    temp_df['anomaly_score2'] = temp_df.apply(set_anomaly_by_score, axis=1)

    column_name = 'isof_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3])
    results_df[column_name] = temp_df['anomaly_score2']

    return results_df

def run_DBScan(eps, distance, csv_file_path, results_df):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance not in unifrac_distances else 'precomputed'

    model = DBSCAN(eps=eps, min_samples=2, metric=diststr)
    model.fit(temp_df)
    temp_df['anomaly_score'] = model.labels_
    temp_df['anomaly_score2'] = temp_df.apply(set_anomaly_by_score, axis=1)

    column_name = 'dbscan_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = temp_df['anomaly_score2']

    return results_df

def run_HDBScan(distance, csv_file_path, results_df):

    temp_df = pd.read_csv(csv_file_path, index_col=[0])
    diststr = distance if distance not in unifrac_distances else 'precomputed'

    model = hdbscan.HDBSCAN(metric=diststr, alpha=0.5, min_cluster_size=2, cluster_selection_method='eom' ,allow_single_cluster=True)
    model.fit(temp_df)
    temp_df['anomaly_score'] = model.labels_
    temp_df['anomaly_score2'] = temp_df.apply(set_anomaly_by_score, axis=1)

    column_name = 'hdbscan_' + '_'.join(csv_file_path.split('/')[-1].split('_')[1:3]) + f'_{distance}'
    results_df[column_name] = temp_df['anomaly_score2']

    return results_df


def Main():

    parser = argparse.ArgumentParser(description='Run ML script to generate outlier classification')
    parser.add_argument("-s", "--study_name", help="STUDY NAME  (required)", required=True)
    parser.add_argument("-r", "--taxa_rank", help="TAXA RANK AGGREGATION (phylum, class, order, family, genus, species)", required=True)
    parser.add_argument("-a", "--transformation", help="OTU ABUNDANCE TRANSFORMATION (counts, relative, clr (Centered Log Ratio), Log10)", required=True)
    parser.add_argument("-d", "--distance_measure", help="DISTANCE TO USE (gunifrac, unifrac, wunifrac, va-wunifrac)", required=True)
    args = parser.parse_args()

    study_name = args.study_name
    taxa_rank = args.taxa_rank
    transformation = args.transformation
    distance_measure = args.distance_measure
    nneighbours = int(args.knn_nneighbours)
    outlier_value = float(args.outlier_value)

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

    ## CLASS 1

    # KNN
    print('RUNNING KNN MODEL CLASS 1')
    class1_df_results = run_KNN(nneighbours, 'minkowski', outlier_value, class1_otu, class1_df_results)
    class1_df_results = run_KNN(nneighbours, distance_measure, outlier_value, class1_distances, class1_df_results)
    # Isolation Forest
    print('RUNNING ISOLATION FOREST MODEL CLASS 1')
    class1_df_results = run_IsolationForest(outlier_value, class1_otu, class1_df_results)
    # SKLEARN DBSCAN
    print('RUNNING DBSCAN MODEL CLASS 1')
    class1_df_results = run_DBScan(0.5, 'euclidean', class1_otu, class1_df_results)
    class1_df_results = run_DBScan(0.5, distance_measure, class1_distances, class1_df_results)
    # HDBSCAN
    print('RUNNING HDBSCAN MODEL CLASS 1')
    class1_df_results = run_HDBScan('euclidean', class1_otu, class1_df_results)
    class1_df_results = run_HDBScan(distance_measure, class1_distances, class1_df_results)

    ## CLASS 2

    # KNN
    print('RUNNING KNN MODEL CLASS 2')
    class2_df_results = run_KNN(nneighbours, 'minkowski', outlier_value, class2_otu, class2_df_results)
    class2_df_results = run_KNN(nneighbours, distance_measure, outlier_value, class2_distances, class2_df_results)
    # Isolation Forest
    print('RUNNING ISOLATION FOREST MODEL CLASS 2')
    class2_df_results = run_IsolationForest(outlier_value, class2_otu, class2_df_results)
    # SKLEARN DBSCAN
    print('RUNNING DBSCAN MODEL CLASS 2')
    class2_df_results = run_DBScan(0.5, 'euclidean', class2_otu, class2_df_results)
    class2_df_results = run_DBScan(0.5, distance_measure, class2_distances, class2_df_results)
    # HDBSCAN
    print('RUNNING HDBSCAN MODEL CLASS 2')
    class2_df_results = run_HDBScan('euclidean', class2_otu, class2_df_results)
    class2_df_results = run_HDBScan(distance_measure, class2_distances, class2_df_results)

    print('EXPORTING RESULTS')
    class1_df_results.to_csv(os.path.join(main_results_folder, 'anomaly_class1.csv'), index=False)
    class2_df_results.to_csv(os.path.join(main_results_folder, 'anomaly_class2.csv'), index=False)


if __name__ == '__main__':
    Main()