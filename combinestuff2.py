import pandas as pd

dfnetwork = pd.read_csv("NRF2DB.csv")
dfglut = pd.read_csv("glutpath.csv")

nrf2network = dfnetwork.merge(dfglut,on='target_name')

nrf2network.to_csv('TotalData4.csv')