import pandas as pd
import joblib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

df = pd.read_csv('SimsN.csv')
L1 = np.array(df['L1'][:])
L2 = np.array(df['L2'][:])
L3 = np.array(df['L3'][:])
T1 = np.array(df['T1'][:])
T2 = np.array(df['T2'][:])
w = np.array(df['Ux'][:])
t = np.array(df['time'][:])

param = np.array([L1,L2,L3,T1,T2,w])
scale = np.empty(shape = (2,6))
paramsc = []

for i in range(len(param)):
  scale[0][i] = np.mean(param[i])
  scale[1][i] = np.std(param[i])
  paramsc.append((param[i]-scale[0][i])/scale[1][i])

paramsc = np.array(paramsc).T
param_train, param_test, t_train, t_test = train_test_split(paramsc, t, test_size=0.3)

model = RandomForestRegressor(n_estimators=100, random_state=100)
model.fit(param_train, t_train)

t_predic = model.predict(param_test)


joblib.dump(model, "RFRmodel.joblib")


e = t_test - t_predic
l = len(t_test)
D = np.array([t_test,t_predic]).T
sorted_indices = np.argsort(D[:, 0])
Dif = D[sorted_indices]

#RMS error
RMS = np.sqrt((np.sum((e)**2))/len(t_test))
print('RMS: ', int(RMS))

#Error sum
esum = np.sum(e)
print('Error sum: ', int(esum))

#RMS for %
Range = [0.05*l, 0.1*l, 0.2*l]
Label = ["5%", "10%", "20%"]
for i in range(3):
  tot = 0
  R = int(Range[i])
  for j in range(0,R):
    tot = tot + (Dif[j,0]-Dif[j,1])**2
  rms = int(np.sqrt(tot/R))
  print("RMS for ",Label[i], ": ", rms)

#Accuracy within a %
Acc = abs(100*e/t_test)
a5,a10,a20,a50 = 0,0,0,0
for i in Acc:
  if i <= 5:
    a5 += 1
    a10 += 1
    a20 += 1
    a50 += 1
  elif i <= 10:
    a10 += 1
    a20 += 1
    a50 += 1
  elif i <= 20:
    a20 += 1
    a50 += 1
  elif i <= 50:
    a50 += 1

print('5% accuracy: ', int(100*a5/l))
print('10% accuracy: ', int(100*a10/l))
print('20% accuracy: ', int(100*a20/l))
print('50% accuracy: ', int(100*a50/l))

#Batch split accuracy
correct = 0
n = 0
for k in range(10,l-10,10):
  minp = 100000
  min = 100000
  for i in range(k-10,k):
    if t_predic[i] < minp:
      minp = t_predic[i]
      p = i
    if t_test[i] < min:
      min = t_test[i]
      true = i
  n += 1
  if true == p:
    correct += 1
print('Correct %: ', 100*correct/n)

#Absolute error as a funtion of time
emean = []
tmean = []
for k in range(20,l-20,20):
  etot = 0
  ttot = 0
  for i in range(k-20,k):
    etot += abs(Dif[i,0]-Dif[i,1])
    ttot += (Dif[i,0])
  emean.append(etot/20)
  tmean.append(ttot/20)

fig, ax = plt.subplots()
tmean = np.array(t_test)
emean = np.array(abs(e)/10)
coefficients = np.polyfit(tmean, emean, 1)
polynomial = np.poly1d(coefficients)
lfit = polynomial(tmean)
ax.set_xlabel('Time', fontsize=10)
ax.set_ylabel('Error', fontsize=10)
ax.set_title('Error as a function of time for the Scikit-Learn RFR model', fontsize=12)
plt.scatter(tmean,emean, s = 10, c = "orange")
plt.plot(tmean,lfit)
