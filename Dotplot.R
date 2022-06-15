df <- read.csv("C:\\Users\\joelc\\OneDrive\\Desktop\\Dotplottrial.csv")

library("ggplot2")
library('easyGgplot2')
ggplot2.dotplot(data=df, xName='X',yName='log2.values', 
                stackratio=0.5, dotsize=0.02)
ggplot(df, aes(x=X, y=log2.values)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.2)





