import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("C:\\Users\\joelc\\OneDrive\\Desktop\\GSEAplot.csv")

# use the scatterplot function
sns.scatterplot(data=df, x="Normalised Enrichment Score", y="Pathways", size="SIZE", hue="NOM p-val",alpha= 0.7, sizes=(40, 100))
plt.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize=8)

# show the graph
sns.set(rc = {'figure.figsize':(8,8)})
plt.xlim(1.2,2.2)
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(left=0.37, right=0.8)
plt.title("GSEA Results", fontsize= 15)
plt.tick_params(axis='y', labelsize=5)
plt.savefig('GSEA.jpeg',dpi=200)
plt.show()
