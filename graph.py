str_of_r=["2","3","4","5"]

s="Results/"

import pandas as pd
import matplotlib.pyplot as plt

names=["Hansen_r","Hill_r","Shekel_r"]
print(len(str_of_r))



### Hansen
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
plt.title("Операционная характеристика Хансен")   
plt.legend()
plt.show()


### Hill
print(names[1]);

data_vec = [""]*(len(str_of_r))
maxx=0
num=len(str_of_r)
for i in range(num):
    #print("r =",i)
    data_vec[i] = pd.read_csv("./"+s+names[1]+str_of_r[i]+".csv")
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
plt.title("Операционная характеристика Хилл")     
plt.legend()
plt.show()


### Sheckel

print(names[2]);

for i in str_of_r:
    #print("r =",i)
    data=pd.read_csv("./"+s+names[2]+i+".csv")
    data.columns=["x","y"]
    plt.plot(data["x"], data["y"], label = f"r = {i}")# , kind="scatter")
    
plt.ylim(bottom=0,top=1) 
plt.xlim(left=0)
plt.title("Операционная характеристика Шекель")    
plt.legend()
plt.show()



### SheckelConstr
names=["SheckConstr_q"]
s="Results/"

print(names[0]);

str_of_q=["1","10","100"]
for i in str_of_q:
    #print("r =",i)
    data=pd.read_csv("./"+s+names[0]+i+".csv")
    data.columns=["x","y"]
    plt.plot(data["x"], data["y"], label = f"q = {i}")# , kind="scatter")
    
plt.ylim(bottom=0,top=1) 
plt.xlim(left=0)
plt.title("Влияние параметра q на число испытаний")    
plt.legend()
plt.show()

### MultiDim

str_of_r=["2","3","4","5", "6"]
names=["GKLS_r","Grishagin_r"]

### GKLS
print(names[0]);

for i in str_of_r:
    #print("r =",i)
    data=pd.read_csv("./"+s+names[0]+i+".csv")
    data.columns=["x","y"]
    plt.plot(data["x"], data["y"], label = f"r = {i}")# , kind="scatter")
    
plt.ylim(bottom=0,top=1) 
plt.xlim(left=0)
plt.title("Операционная характеристика GKLS")        
plt.legend()
plt.show()

### Grishagin
print(names[1]);

for i in str_of_r:
    #print("r =",i)
    data=pd.read_csv("./"+s+names[1]+i+".csv")
    data.columns=["x","y"]
    plt.plot(data["x"], data["y"], label = f"r = {i}")# , kind="scatter")
    
plt.ylim(bottom=0,top=1) 
plt.xlim(left=0)
plt.title("Операционная характеристика Гришагин")    
plt.legend()
plt.show()
