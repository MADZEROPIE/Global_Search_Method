import pandas as pd
import matplotlib.pyplot as plt

names=["Hansen_r","Hill_r","Shekel_r"]
print(len(str_of_r))

print(names[0]);

data_vec = [""]*(len(str_of_r))
maxx=0
num=len(str_of_r)
for i in range(num):
    #print("r =",i)
    data_vec[i] = pd.read_csv("./"+s+names[0]+str_of_r[i]+".csv")
    data_vec[i].columns=["x","y"]
    if data_vec[i]["x"][len(data_vec[i]["x"])-1]>maxx:
        maxx=data_vec[i]["x"][len(data_vec[i]["x"])-1]
        
print(maxx)

for i in range(len(str_of_r)):
    #data1=pd.Series([[maxx, [data_vec[i]["y"][len(data_vec[i]["y"])-1]]]) # ???
    #data1.columns=["x","y"] # ???
    # NEED TO PUSH BACK. SOMEHOW.
    data_x=list(data_vec[i]["x"])
    data_y=list(data_vec[i]["y"])
    data_x.append(maxx)
    data_y.append(data_vec[i]["y"][len(data_vec[i]["y"])-1])
    
    plt.plot(data_x, data_y, label = f"r = {str_of_r[i]}")# , kind="scatter")
    
plt.ylim(bottom=0,top=1) 
plt.xlim(left=0)
    
plt.legend()
plt.show()

