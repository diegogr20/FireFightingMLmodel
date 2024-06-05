import numpy as np

def Parameters():
    l1 = np.linspace(10,130, 3)
    l2 = np.linspace(0,630, 31)
    l3 = np.linspace(50,100, 2)
    t1 = np.linspace(0,1,2)

    l1grid, l2grid, l3grid, t1grid = np.meshgrid(l1,l2,l3,t1)
    X = np.array([l1grid,l2grid,l3grid,t1grid]).reshape([4,len(t1)*len(l1)*len(l2)*len(l3)]).T

    c = 3.17

    x = np.linspace(0,100,1000)
    ytot = 2*x*np.exp(-x/(c**2))/(c**2)
    y = ytot/np.trapz(ytot)

    wind = np.array([])
    dhwind = 5
    Ntot = 20

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

    Xf = np.empty(5)
    for j in range(len(winds)):
      wdummy = np.ones(len(X))*winds[j]
      Xf = np.vstack([Xf,np.concatenate([X[:,0],X[:,1],X[:,2],X[:,3],wdummy]).reshape(5,len(X)).T])

    Xf = np.delete(Xf,0,0)
    return Xf


