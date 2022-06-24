import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 
from matplotlib import cm
from matplotlib.animation import FuncAnimation

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

def ds_lump(x, y):
    c = 1
    # x = x
    # y = y-3j
    num = 2*c*(x+1j*y)
    den = np.abs(c)**2*(1j*y+x/2)-1
    p = num/den
    return p

def ds_lump2(x, y):
    c = 1
    # x = x
    # y = y-3j
    num = 2*c*(x+1j*y)
    den = np.abs(c)**2*(1j*y+x/2)**2-1
    p = num/den
    return p

def ds_simple(x, y):
    c = 1
    # x = x
    # y = y-3j
    num = 2*c*(x+1j*y)
    den = np.abs(c)**2*(1j*y+x/2)**2*x**2+1
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
ds = ds_lump
Z = ds(X,Y)*phase
Z = np.real(Z)

ax.plot_surface(X,Y,Z, cmap=cm.coolwarm)
ax.set_xlabel('x')
ax.set_ylabel('y')



ax = fig.add_subplot(gs[0, 1] )
ax.imshow(Z)
ax.set_xticks([])
ax.set_yticks([])


# ax = fig.add_subplot(gs[1,0])
# ax.plot(x, np.real(ds(x,0)))
# ax.set_xlabel('x')
# ax.set_title('when y=0')

# ax = fig.add_subplot(gs[1,1])
# ax.plot(y, np.real(ds(0,y)))
# ax.set_title('when x=0')
# ax.set_xlabel('y')

plt.show()


def make_movie(ds):
    global X,Y
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    lam1 = 1+1j
    X = X-t
    phase = np.exp(1j*(X + Y))
    ds = ds_simple
     #  projection='3d'
    def animate(i):
        # global X,Y
        # phase = (i/60)*(np.pi*1j/4) + 3*np.pi/4
        # phase = (i/60)*2*np.pi*1j
        # lam1 = np.exp(phase)
        t = -10+i
        # lam1 = 1+1j
        
        Z = ds(X,Y)*phase
        Z = np.real(Z)
        # Z = np.real(build_solution( tau , t=t, lam=lam1, eta=eta ))
        ax.clear()
        # ax.set_zlim(-5, 30)
        # plt.title(r't={}, $\lambda$={},\eta={}$'.format(t, lam1, eta))
        # surf = ax.plot_surface(X[1:-1], Y[1:-1], Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        # cont = ax.plot_wireframe(X[1:-1],Y[1:-1],Z)
        # cont = ax.contour(X[1:-1], Y[1:-1], Z, 100) # levels=list(range(-10,10))
        # cont = ax.contour(X[1:-1], Y[1:-1], Z, 100, cmap="RdGy") # levels=list(range(-10,10))
        # cont = ax.contourf(X,Y,Z, 50) # levels=[-10, 0, 10]
        cont = plt.imshow(np.real(Z), cmap='hot', interpolation='nearest')
        # ax = fig.add_subplot(1, 2, 2)
        # ax.clear()
        # ax.set_xlim(-1.2, 1.2)
        # ax.set_ylim(-1.2, 1.2)
        # ax.axes(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
        # plt.plot(np.real(lam1), np.imag(lam1), 'k.')
        return cont

    anim = FuncAnimation(fig, animate, frames=20, interval=1) # , init_func=init, blit=True, , interval=1
    # anim.save('depth1_scattering_soliton1_lam={}_eta={}.gif'.format(lam1, eta), writer='Pillow')
    plt.show()

# make_movie(ds)