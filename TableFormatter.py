import pandas as pd

dfsampleinfo = pd.read_csv("HCCDB13sample.csv")
dffilter = pd.read_csv("HeatmapFilter.csv")
dfexpression = pd.read_csv("HCCDB13Expression.csv")

# Filter HCC Samples
dfsampleinfo.set_index('SAMPLE_ID',inplace=True)
dfsampleinfo = dfsampleinfo.T
dfsampleinfo = dfsampleinfo[dfsampleinfo['TYPE'] == "HCC"]
dfsampleinfo['SAMPLE_ID'] = dfsampleinfo.index


cancersamples = dfsampleinfo["SAMPLE_ID"].tolist()
cancersamples = ["Symbol"] + cancersamples
print(cancersamples)
dfexpression = dfexpression.filter(cancersamples)

print(dfexpression.head())
dfcancerexpressions = dffilter.merge(dfexpression, on="Symbol")
dfcancerexpressions.to_csv('CancerLiver13Data.csv')