import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 


locs = pd.read_csv("results/intersect_locs.csv",header=0,names=["start","end"]).values
nlocs = locs - locs.min()
expression = pd.read_csv("results/intersect_expression.csv",header=0,names=["bedid","expr"]).values[:,1]

print(expression.shape,locs.shape)

score = np.zeros(nlocs.max()+1)
for i, expr in enumerate(expression):
    if expr != -1:
        score[nlocs[i,0]:nlocs[i,1]+1] += expr
np.savetxt("results/score.txt",score)


# fig, ax = plt.subplots()
# ax.plot(score)
# fig.savefig("results/score.pdf")