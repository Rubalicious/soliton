import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

N = 1000
L = 20
x = np.linspace(-L,L,N)
y = np.linspace(-L,L,N)
X,Y = np.meshgrid(x,y)

# dromion parameters
m1 = 1
m2 = 1j
c1 = 1
c2 = 0



def dromion(x,y):
    # x = x-5j
    # y = 1j*y
    global m1, m2, c1, c2
    num = -m2*m1*(y**2+x**2) + c2*m1*(1j*y+x) + c1*m2*(1j*y-x) + c1*c2
    A = -m1**2*(y**2*x + x**3/3) + (c1**2 + 2*c1*m1*1j*y)*x 
    B = -m2**2*(y**2*x + x**3/3) + (c2**2 + 2*c2*m2*1j*y)*x
    den = A*B+1
    return num/den


# regular and singular solutions
def ds(x,y):
    global m1, m2, c1, c2
    num = -m1*m2*(y**2+x**2)+(c1*m2+c2*m1)*1j*y +(c1*m2-c2*m1)*x
    A = -m1*m2*(y**2*x + x**3)  + (c1*m2-c2*m1)*x**2/2 + (c1*m2+c2*m1)*1j*x*y
    den = A**2+1
    # log = np.log(num)-np.log(den)
    # experiment = np.exp(log)
    return num/den



ds = ds
Z = ds(X,Y)
Z = np.real(Z)
# Z = (X**2+Y**2)/((X**2+Y**2)**2 + 1)
# e = 20
# Z = Z[e:N-e,e:N-e]
# X = X[e:N-e,e:N-e]
# Y = Y[e:N-e,e:N-e]
# print(Z[mask])
fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(221)
plt.imshow(Z)
plt.title('top view')
# plt.show()

# def den(x, y):
#     return (y**2*x+x**3/3)**2+1
# plt.imshow(den(X,Y))
# plt.plot_wireframe(X,Y,Z)
# fig = plt.figure()
ax = fig.add_subplot(222, projection='3d')
x0=0+.1
y0=0
# ax.plot(x, ds(x,y0), zs=y0, zdir='y', c='k', linewidth=1)
# ax.plot(y, ds(x0,y), zs=x0, zdir='x', c='k', linewidth=1)
ax.plot_surface(X,Y,Z,cmap=cm.coolwarm)
# ax.set_zlim([-10, 40])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('3D view')

# plt.show()




# Zy = ds(0,y)
Zx = ds(x,y0)
# fig  = plt.figure(figsize=(10,5))
ax = fig.add_subplot(223)
plt.plot(x, ds(x,y0))
plt.title('x cross section')
ax = fig.add_subplot(224)
plt.plot(y, ds(x0,y))
plt.title('y cross section')

# eps = 0.01
# py = y**2/(eps**2*y**4+1)
# plt.plot(y,py)
plt.savefig('./figures/ds/ds_m1={}m2={}c1={}c2={}.png'.format(m1, m2, c1, c2))
plt.show()
