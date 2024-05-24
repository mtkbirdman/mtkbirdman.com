import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import random
import math
from scipy.interpolate import interp1d

class CST:
  def __init__(self):
    self.wl = None
    self.wu = None
    self.dz = None
    self.Node = None
    self.x = None
    self.y = None
    self.coord = None
    self.x_old = None
    self.y_old = None
    self.xt = None
    self.yt = None
    self.ler = None
    self.n_wl = None
    self.n_wu = None

  def _class_shape(self, w, x, N1, N2, dz):
    C = x**N1 * (1 - x)**N2 # Class function
    # Shape function; using Bernstein Polynomials
    S = np.zeros_like(x)
    nw = w.shape[1] - 1 # Order of Bernstein polynomials
    K = np.zeros(nw+1)
    for i in range(0, nw+1):
      K[i] = math.factorial(nw)/(math.factorial(i)*(math.factorial((nw)-(i))))
    for j in range(0, nw+1):
      S += (w[:,j]*K[j]).reshape(-1,1)*x**(j) * ((1-x)**(nw-(j)))
    y = C * S + x * dz
    return y

  def create_airfoil(self, wu=[[-1, -1, -1]], wl=[[1, 1, 1]], dz=[0], Node=101):
    self.wl = np.array(wl)
    self.wu = np.array(wu)
    self.dz = np.array(dz).reshape(-1,1)
    self.Node = Node
    N = self.dz.shape[0]

    if self.x is None: # Create x coordinate
      self.x = np.ones((N,self.Node))
      zeta = np.zeros((N,self.Node))
      for i in range(self.Node):
        zeta[:, i] = 2 * np.pi / (self.Node - 1)*i
        self.x[:, i] = 0.5 * (np.cos(zeta[:, i])+ 1)
    # Ni and N2 parameters (N1 = 0.5 and N2 = 1 for airfoil shape)
    N1 = 0.5
    N2 = 1

    center_loc = np.argmin(self.x)
    xu = self.x[:,:center_loc]
    xl = self.x[:,center_loc:]
    yu = self._class_shape(self.wu, xu, N1, N2, self.dz)
    yl = self._class_shape(self.wl, xl, N1, N2, -self.dz)

    self.x = np.concatenate([xu, xl],axis=1).reshape(N,-1)
    self.y = np.concatenate([yu, yl],axis=1).reshape(N,-1)
    self.coord = np.stack([self.x, self.y],axis=1)
    return self.coord

  def write_to_file(self, file_name = 'airfoil_shape.csv', N=0):
    df = pd.DataFrame(self.coord[N].T, columns=['x', 'y'])
    df.to_csv(file_name, header=True, index=False)

  def fit_CST(self, file_name = 'airfoil_shape.csv', n_wl=4, n_wu=4):
    data = pd.read_csv(file_name)
    self.x_old = data['x'].to_numpy()
    self.y_old = data['y'].to_numpy()
    self.Node = len(data)
    self.dz = [abs(self.y_old[-1])]
    self.n_wl = n_wl
    self.n_wu = n_wu
    # normalization
    self.x_old = self.x_old/(np.max(self.x_old)-np.min(self.x_old))
    self.x_old = self.x_old-np.min(self.x_old)
    self.y_old = self.y_old/(np.max(self.x_old)-np.min(self.x_old))
    self.x = self.x_old.reshape(1,-1)
    initial_guess = [random.uniform(-1,1) for _ in range(self.n_wl+self.n_wu)]
    result = minimize(self._objfunc, initial_guess) # minimize
    
    self.x=None
    self.create_airfoil(wl=self.wl, wu=self.wu, dz=self.dz, Node=self.Node)

  def _objfunc(self, x):
    wl = [x[:self.n_wl]]
    wu = [x[self.n_wl:self.n_wl + self.n_wu]]
    self.create_airfoil(wl=wl, wu=wu, dz=[0], Node=self.Node)
    f = np.sum((self.y * 100 - self.y_old * 100)**2)
    return f

  def get_var(self):
    return self.wu, self.wl

  def get_thickness(self):
    thickness = []
    x_tmax = []
    self.xt = np.zeros_like(self.x)
    self.yt = np.zeros_like(self.y)
    for N in range(self.coord.shape[0]):
      center_loc = np.argmin(self.coord[N,0])
      coord_u = self.coord[N, :, :center_loc+1]
      coord_l = self.coord[N, :, center_loc:]
      
      fl = interp1d(coord_l[0,:], coord_l[1,:], kind= 'linear')
      self.xt[N,:center_loc+1] = coord_u[0,:]
      self.yt[N,:center_loc+1] = coord_u[1,:]-fl(coord_u[0,:])

      fu = interp1d(coord_u[0,:], coord_u[1,:], kind='linear')
      self.xt[N,center_loc:] = coord_l[0,:]
      self.yt[N,center_loc:] = fu(coord_l[0,:])-coord_l[1,:]
      thickness.append(np.max(self.yt[N]))
      x_tmax.append(self.xt[N,np.argmax(self.yt[N])])

    return np.array(thickness),np.array(x_tmax)
  
  def get_ler(self):
    if self.yt is None: 
      self.get_thickness()
    self.ler = np.zeros(self.coord.shape[0])
    for N in range(self.coord.shape[0]):
      x0 = np.min(self.xt[N])
      y0 = self.yt[N,np.argmin(self.xt[N])]/2
      r = ((self.xt[N]-x0)**2+(-y0)**2)**0.5
      center_loc = np.argmin(self.xt[N])
      for i in range(center_loc,r.shape[0]):
        d = ((self.xt[N,:]-self.xt[N,i])**2+(self.yt[N,:]/2-y0)**2)**0.5

        if r[i]>np.min(d):
          self.ler[N] = r[i-1]
          break
    return self.ler
  
  def get_TE_angle(self):
    if self.yt is None: 
      self.get_thickness()
    TE_angle = []
    for N in range(self.coord.shape[0]):
      center_loc = np.argmin(self.coord[N,0])
      f = interp1d(self.xt[N,:center_loc+1], self.yt[N,:center_loc+1], kind= 'linear')
      x96 = 0.96
      x99 = 0.99
      y96 = f(x96)
      y99 = f(x99)
      #print(y96,y99,x96,x99)
      TE_angle.append(math.degrees(math.atan2(abs(y96-y99),abs(x96-x99))))

    return np.array(TE_angle)

  def get_section_area(self):
    if self.yt is None: 
      self.get_thickness()
    
    x_diff = np.abs(np.diff(self.xt,n=1,axis=1))
    section_area = np.sum(x_diff * self.yt[:,1:],axis=1)/2

    return section_area
  
  def plot(self,N=0):
    fig, ax = plt.subplots()
    ax.plot(self.x[N], self.y[N],label='CST')
    if self.y_old is not None:
      ax.plot(self.x_old, self.y_old,label='original')
    if self.ler is not None:
      theta = np.linspace(0, 2 * np.pi, 100)
      #中心と半径から円の座標を計算
      x = self.ler[N] + self.ler[N] * np.cos(theta)
      y = self.ler[N] * np.sin(theta)
      ax.plot(x, y,label='ler')
      ax.plot(self.xt[N], self.yt[N]/2,label='Symmetry')
    ax.set_aspect("equal", "box")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel( 'x/c')
    ax.set_ylabel( 'y/c')
    ax.spines[ 'right'] .set_visible(False)
    ax.spines[ 'top'] .set_visible(False)
    ax.yaxis.set_ticks_position( 'left')
    ax.xaxis.set_ticks_position( 'bottom')
    ax.legend()
    ax.plot
    plt.show()

if __name__ == '__main__':
  airfoil = CST()
  airfoil.fit_CST(file_name='./CST/NACA4412.csv')
  thickness, x_tmax = airfoil.get_thickness()
  ler = airfoil.get_ler()
  TE_angle = airfoil.get_TE_angle()
  section_area = airfoil.get_section_area()
  airfoil.plot(N=0)
  wu = airfoil.wu
  wl = airfoil.wl
  
  print(wu,wl)
  print(thickness,x_tmax,ler,TE_angle)
  print(section_area)
  
  airfoil.write_to_file()