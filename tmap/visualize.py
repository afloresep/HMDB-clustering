import pickle
import pandas as pd
prep_tmap = pd.read_pickle('HMDB_prep_tmap.pkl')
prep_tmap.to_csv('HMDB_prep_tmap.csv', index=False)
print(prep_tmap.head())
