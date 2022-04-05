import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 
from matplotlib import cm

N = 200
L_x = 10
L_y = 10
x = np.linspace(-L_x, L_x, N)
y = np.linspace(-L_y, L_y, N)


X,Y = np.meshgrid(x, y)

def ds(x, y):
    num = (x-y)*(x+y)
    den = y**4*x**2-y**2*x**4/3 +x**6/9 + 1
    p = num/den
    return p

def ds_c(x, y):
    # y = 1j*y
    # x = 1j*x
    num = -2*(x**2+y**2)
    den = (y**2*x+x**3/3)**2 + 1
    p = num/den
    return p

def ds_const(x,y):
    c = 1
    d = 1
    num = 2*np.conj(c)*np.conj(d)*y
    den = (np.abs(c*d)*x**2 +1)*y
    p = num/den
    return p


def ds_simple(x, y):
    c = 1
    # x = x
    # y = y-3j
    num = 2*c*(x+1j*y)
    den = np.abs(c)**2*(1j*y+x/2)**2*x**2-1
    p = num/den
    return p

def ds_simple_alt(x, y):
    c = 1
    # x = x
    # y = y-3j
    num = 2*np.conj(c)*(-x+1j*y)
    den = np.abs(c)**2*(1j*y-x/2)**2*x**2-1
    p = num/den
    return p


fig = plt.figure(figsize=(9,5))
gs = GridSpec(1,2, figure = fig)
ax = fig.add_subplot(gs[0, 0] , projection='3d')

t = 0
X = X-t
phase = np.exp(1j*(X + Y))
ds = ds
Z = ds(X,Y)*phase
Z = np.real(Z)

ax.plot_surface(X,Y,Z, cmap=cm.coolwarm)
ax.set_xlabel('x')
ax.set_ylabel('y')



ax = fig.add_subplot(gs[0, 1] )
ax.imshow(Z)


# ax = fig.add_subplot(gs[1,0])
# ax.plot(x, np.real(ds(x,0)))
# ax.set_xlabel('x')
# ax.set_title('when y=0')

# ax = fig.add_subplot(gs[1,1])
# ax.plot(y, np.real(ds(0,y)))
# ax.set_title('when x=0')
# ax.set_xlabel('y')

plt.show()