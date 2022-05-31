import pandas as pd

dfsampleinfo = pd.read_csv("BRCA_Samples.csv")
dffilter = pd.read_csv("HeatmapFilter.csv")
dfexpression = pd.read_csv("TCGA_BRCA_Expression.csv")

# Filter HCC Samples
dfsampleinfo.set_index('sample',inplace=True)
dfsampleinfo = dfsampleinfo[dfsampleinfo['sample_type'] == "Primary Tumor"]
dfsampleinfo['sample'] = dfsampleinfo.index


cancersamples = dfsampleinfo["sample"].tolist()
cancersamples = ["Symbol"] + cancersamples
print(cancersamples)
dfexpression = dfexpression.filter(cancersamples)

print(dfexpression.head())
dfcancerexpressions = dffilter.merge(dfexpression, on="Symbol")
dfcancerexpressions.to_csv('AllCancerBreastData.csv')
