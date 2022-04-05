import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

N = 1000
L = 50
x = np.linspace(-L,L,N)
y = np.linspace(-L,L,N)
X,Y = np.meshgrid(x,y)

# x = x*np.exp(1j*np.pi/4)
# y = y/2

# test integration
# y = np.trapz( x - 1)
# print(len(y))
# plt.plot(x,y)
# plt.show()
# quit()
def ds2soln(x,y):
    A = x**3
    B = x**4
    C = x**5
    D = x**3
    E = x**4
    F = x**5
    psi11 = -2*x*(C+(1+x)**2)
def dromion(x,y):
    num = (y**2-x**2)
    den = (y**2*x+x**3/3)**2-x**4*y**2-1
    return num/den

def linear_ds(x, y):
    m1 = -1
    m2 = -3
    c1 = 5
    c2 = -4
    num = m2*m1*(x+y)**2+(c1*m2+c2*m1)*(x+y)+c1*c2
    den = x**2*(m1*m2*(y**2-x**2/3)+(c1*m2+c2*m1)*y+c1*c2)**2 - 1
    return num/den

def ds1(x,y):
    num = y**2 - x**2
    den = x**2*(y**2 - 2*x**2) - 1
    return num/den

def ds(x,y):
    num = (x**2 + y**2)
    den = (y**2*x + x**3/3)**2 + 1
    # log = np.log(num)-np.log(den)
    # experiment = np.exp(log)
    return num/den

def ds2(x, y):
    num = 2*x+1j*y
    den = (1j*y+.5*x)**2*x**2+1
    return num/den

ds = ds
Z = ds(X,Y)
# Z = np.real(Z)
# Z = (X**2+Y**2)/((X**2+Y**2)**2 + 1)
# e = 20
# Z = Z[e:N-e,e:N-e]
# X = X[e:N-e,e:N-e]
# Y = Y[e:N-e,e:N-e]
# print(Z[mask])
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(121)
plt.imshow(Z)
# plt.show()

# def den(x, y):
#     return (y**2*x+x**3/3)**2+1
# plt.imshow(den(X,Y))
# plt.plot_wireframe(X,Y,Z)
# fig = plt.figure()
ax = fig.add_subplot(122, projection='3d')
x0=0+.1
y0=0
ax.plot(x, ds(x,y0), zs=y0, zdir='y', c='k', linewidth=1)
ax.plot(y, ds(x0,y), zs=x0, zdir='x', c='k', linewidth=1)
ax.plot_surface(X,Y,Z,cmap=cm.coolwarm)
ax.set_zlim([-10, 10])
ax.set_xlabel('x')
ax.set_ylabel('y')

plt.show()




# Zy = ds(0,y)
Zx = ds(x,y0)
fig  = plt.figure(figsize=(10,5))
ax = fig.add_subplot(121)
plt.plot(x, ds(x,y0))
ax = fig.add_subplot(122)
plt.plot(y, ds(x0,y))

# eps = 0.01
# py = y**2/(eps**2*y**4+1)
# plt.plot(y,py)
plt.show()
