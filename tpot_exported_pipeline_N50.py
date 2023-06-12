import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler


def N50_recommend(metaFeature):
    # NOTE: Make sure that the outcome column is labeled 'target' in the data file
    reader1 = pd.read_csv('data/data_N50.csv')

    x = reader1[['x1','x2','x3','x4','x5','x6','x7','x8','x9','x10']]
    y = reader1['y']

    x_train, x_test,y_train, y_test = train_test_split(x,y, random_state=68,train_size=0.8)


    exported_pipeline = make_pipeline(
        MinMaxScaler(),
        RandomForestClassifier(bootstrap=False, criterion="gini", max_features=0.7000000000000001, min_samples_leaf=3, min_samples_split=9, n_estimators=100)
    )

    exported_pipeline.fit(x_train, y_train)

    #meta_features = test.metaFeature
    results_N50 = exported_pipeline.predict(metaFeature)
    #print("Predicted:",results)
    resultDic = {0: "Nanosv,flye",1:"picky,flye",2:"Pbsv,smartdenovo",3:"picky,mecat2",4:"Pbsv,mecat2",5:"Pbsv,flye",6:"picky,smartdenovo",7:"Nanosv,smartdenovo",
                 8:"Nanosv,mecat2",9:"picky,wtdbg2",10:"Pbsv,wtdbg2",11:"Nanosv,wtdbg2",12:"Nanosv,canu",13:"Pbsv,canu",14:"picky,canu"}

    print("the highest N50 software combination:",resultDic[results_N50[0]])
