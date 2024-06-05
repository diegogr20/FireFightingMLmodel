import numpy as np

def Parameters():
    l1 = np.linspace(10,7000, 10)
    l2 = np.linspace(0,7000, 11)
    l3 = np.linspace(0,7000, 11)
    t1 = np.linspace(0,1,2)
    t2 = np.linspace(0,2,3)

    l1grid, l2grid, l3grid, t1grid, t2grid = np.meshgrid(l1,l2,l3,t1,t2)
    Xgrid = np.array([l1grid,l2grid,l3grid,t1grid,t2grid]).reshape([5,len(t1)*len(t2)*len(l1)*len(l2)*len(l3)]).T
    X = np.array([Xgrid[0]])

    for i in range(1,len(Xgrid)):
      length = Xgrid[i,0] + Xgrid[i,1] + Xgrid[i,2]
      if (Xgrid[i,3] != Xgrid[i,4] or (Xgrid[i,3] == Xgrid[i,4] and Xgrid[i,1] == 0 and Xgrid[i,2] == 0)) and length <= 7000 and length > 0 and (Xgrid[i,1] != 0 or (Xgrid[i,1] == 0 and Xgrid[i,4] == 0)):
        X = np.vstack([X,Xgrid[i]])

    c = 3.17

    x = np.linspace(0,100,1000)
    ytot = 2*x*np.exp(-x/(c**2))/(c**2)
    y = ytot/np.trapz(ytot)

    wind = np.array([])
    dhwind = 5
    Ntot = 25

    rang = int(np.max(x)/dhwind)
    h = len(y)/rang
    it = np.linspace(0,len(y)-h,rang)
    count = 5
    for i in it:
      ydummy = np.trapz(y[int(i):int(i+h)])
      n = Ntot*ydummy
      if n < 1:
        n = 1
      else:
        n = int(n)
      windrange = np.linspace(count, count+dhwind,n)
      count = count + dhwind
      wind = np.concatenate([wind,windrange])
    winds = list(set(wind))
    winds.sort()

    Xf = np.empty(6)
    for j in range(len(winds)):
      wdummy = np.ones(len(X))*winds[j]
      Xf = np.vstack([Xf,np.concatenate([X[:,0],X[:,1],X[:,2],X[:,3],X[:,4],wdummy]).reshape(6,len(X)).T])

    Xf = np.delete(Xf,0,0)
    return Xf


