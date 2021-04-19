import pandas as pd
import matplotlib.pyplot as plt

names=["Hansen","Hill","Shekel"]
str_of_r=["2.0","3.0","4.0"]

s="res"

print(names[0]);

for j in range(0,3):
    for i in str_of_r:
    #print("r =",i)
        data=pd.read_csv("./"+s+i+names[j]+".csv")
        data.columns=["x","y"]
        plt.plot(data["x"], data["y"], label = f"r = {i}")# , kind="scatter")
        #plt.step(data["x"],data["y"])
    plt.ylim(bottom=0,top=1) 
    plt.xlim(left=0)
    
    plt.legend()
    plt.show()

