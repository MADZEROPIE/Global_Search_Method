import pandas as pd
import matplotlib.pyplot as plt



data=pd.read_csv("./build/Grishagin0.csv")
data.columns=["x","y"]
plt.scatter(data["x"], data["y"],sizes=[2.5]*len(data["x"]), marker='ш')# , kind="scatter")
        #plt.step(data["x"],data["y"])
#plt.ylim(bottom=0,top=1) 
#plt.xlim(left=0)
    
#plt.legend()
plt.show()

