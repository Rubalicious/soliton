import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button
# plt.rcParams['text.usetex'] = True

# Define domain parameters
N = 200
window = 10
Lx = 100
Ly = 100
x = np.linspace(-window, window, N)
y = np.linspace(-window, window, N)
X,Y = np.meshgrid(x,y)
print(np.shape(X))

# define eigenvalues
eps = .2 # perturbation
lam = 1 # +eps*1j #np.exp(1j*np.pi*(1/4))
eta = 1.2


# phase
def phi(x, y, t, lam):
    return lam*x + 1j*lam**2*y - 4*lam**3*t

# polynomials of depth s = 1, 2, 3
def p1(x, y, t, lam):
    return x + 2j*lam*y - 12*lam**2*t

def p2(x, y, t, lam):
    return 2j*y-24*lam*t + np.power( p1(x, y, t, lam) , 2 )

def p3(x, y, t, lam):
    return np.power( p1(x, y, t, lam) , 3 ) + 3*np.power( p1(x, y, t, lam) , 2 )*(2j*y-24*lam*t) - 72*p1(x, y, lam, t)

# depth 1 lump
def tau(x, y, t, lam):
    alpha = 1/(2*lam)
    phase = np.exp(2*np.real(phi(x, y, t, lam + eta)))
    polynomial = np.power(x-12*lam**2*t-1/lam, 2 ) + 4*lam**2*np.power(y,2)+1/(4*lam**2)
    return alpha*phase*polynomial

# depth 1 lump
def tau_til_depth1(x, y, t, lam):
    # polynomial factor of tau
    return np.power( x - 12*lam**2*t - 1/lam, 2 ) + 4*lam**2*np.power(y , 2) + 1/(4*lam**2)


# depth 2 lump
def tau_til_depth2(x, y, t, lam):
    k = y-6*lam*t
    B = 4*lam*k
    C = 4*lam*(lam*np.power(k, 2)-6*t)
    B2C = np.power(B, 2)+ 2*C

    fac = B - 1/lam

    alpha = 2*fac
    beta = (-3/lam)*fac + B2C
    gamma = (3/lam**2)*fac - (1/lam)*B2C + 2*B*C
    delta = (-3/(2*lam**3))*fac + 2*B2C/(2*lam)**2 - B*C/lam + np.power(C, 2) + 4*np.power(y, 2)
    quartic = np.power(x, 4) + alpha*np.power(x, 3) + beta*np.power(x, 2) + gamma*x + delta
    return quartic

def depth1_scattering(x, y, t, lam, eta):
    t = tau(x, y, t, lam) + tau(x, y, t, eta)
    # tau_le = np.cos((lam**2-eta**2)*y)*(np.power(x, 2) - (2/(lam+eta) - 12*(lam**2+eta**2)*t)*x + 2/(lam+eta)**2 +12*(lam**2+eta**2)*t/(lam+eta)+(12*lam*eta*t)**2+4*eta*lam*np.power(y, 2) )
    # tau_le += 2*(lam-eta)*y*np.sin((lam**2-eta**2)*y)*(x-1/(lam+eta)+12*lam*eta*t)
    # t+= tau_le
    return t

# Rank 2 depth 1 soliton solutions - hr is higher rank
def tau_til_depth1_hr(x, y, t, lam, eta):
    t = tau(x, y, t, lam)*tau(x, y, t, eta)
    # add a mixture term
    Rphi = np.real(phi(x, y, t, lam))
    alpha = 1/(lam+eta)
    Rt = np.power(x, 2) - (2/(lam+eta)+12*(lam**2+eta**2)*t)*x
    Rt+= 2/((lam+eta)**2) + 4*eta*lam*np.power(y, 2)
    Rt+= 12*(lam**2+eta**2)*t/(lam+eta) + (12*lam*eta*t)**2
    It = 2*(lam-eta)*y*(x+12*eta*lam*t-1/(lam+eta))
    t -= np.exp(2*Rphi)*(np.power(Rt, 2) + np.power(It, 2))/(alpha**2)
    return t

# should use a decorator here
# this should be a wrapper for the tau function we use
def build_solution(tau, t=0, lam=1, eta=1):

    # TAU = tau_til_depth1(X, Y, 0, lam)
    TAU = tau(X, Y, t, lam, eta) #

    logTAU = np.log(TAU)
    U = 2*np.diff(logTAU, n=2, axis=0)

    # TAUX = np.diff(TAU, n=1, axis=0)
    # TAUXX = np.diff(TAU, n=2, axis=0)
    # # resize TAU and TAUX to dimensions of TAUXX
    # TAUX = TAUX[1:]
    # TAU = TAU[1:-1]
    # # gives a different result for some reason...
    # U2 = -2*(np.power(TAUX, 2) - TAU*TAUXX)/np.power(TAU, 2)
    return U #,U2

def plot_evolution(tau):
    before = build_solution(tau, t=-5)
    now = build_solution(tau, t=0)
    after = build_solution(tau, t=5)

    fig = plt.figure()

    axes = fig.add_subplot(131)
    axes.imshow(before)

    axes = fig.add_subplot(132)
    axes.imshow(now)

    axes = fig.add_subplot(133)
    axes.imshow(after)
    plt.show()


def plot_heatmap(tau):
    fig = plt.figure()
    window = 20
    ax = plt.axes(xlim=(100-window, 100+window), ylim=(100-window, 100+window) ) # , projection='3d'
    img = plt.imshow(build_solution(tau))
    plt.show()

def time_slider_widget(tau):
    fig = plt.figure()
    window = 10
    ax = plt.axes(xlim=(0, Lx), ylim=(100-window, 100+window) ) # , projection='3d'
    img = plt.imshow(build_solution(tau))
    plt.subplots_adjust(bottom=0.25)
    axtime = plt.axes([0.25, 0.1, 0.65, 0.03])
    time_slider = Slider(
        ax=axtime,
        label='time',
        valmin=-10,
        valmax=0,
        valinit=-5,
    )
    def update(val):
        img.set_data(build_solution(tau, t= time_slider.val, lam=lam))
        fig.canvas.draw_idle()
    time_slider.on_changed(update)
    plt.show()


PHI = phi(X, Y, 0, lam)
P1 = p1(X, Y, 0, lam)
P2 = p2(X, Y, 0, lam)
P3 = p3(X, Y, 0, lam)

plot_heatmap(depth1_scattering)
time_slider_widget(depth1_scattering)

# U1 = build_solution(tau_til_depth1)
# plt.imshow(np.real(U1))
#
# plt.show()
# plt.title("d^2_x log(tau)")
# plt.figure()
# plt.imshow(U2)
# plt.show()
# U2 = build_solution(eta)
# axes = fig.add_subplot(121)

# ax.plot_surface(X[1:-1], Y[1:-1], U1,linewidth=0, antialiased=False)
# ax.plot_surface(X[1:-1], Y[1:-1], U1,linewidth=0, antialiased=False)
# ax.plot_wireframe(X[1:-1], Y[1:-1], U1,linewidth=0.2, antialiased=False)
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# axes = fig.add_subplot(122)

# def animate(i):
#     Z = build_solution( tau_til_depth1_hr , t=i-25, lam=lam, eta=eta )
#     ax.clear()
#     # ax.set_zlim(-5, 30)
#     plt.title(r't={}, $\lambda$={},\eta={}$'.format(i, lam, eta))
#     # surf = ax.plot_surface(X[1:-1], Y[1:-1], Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#     # cont = ax.plot_wireframe(X[1:-1],Y[1:-1],Z)
#     # cont = ax.contour(X[1:-1], Y[1:-1], Z, 100) # levels=list(range(-10,10))
#     # cont = ax.contour(X[1:-1], Y[1:-1], Z, 100, cmap="RdGy") # levels=list(range(-10,10))
#     # cont = ax.contourf(X,Y,Z, 50) # levels=[-10, 0, 10]
#     cont = plt.imshow(np.real(Z), cmap='hot', interpolation='nearest')
#     return cont
#
# anim = FuncAnimation(fig, animate, frames=50) # , init_func=init, blit=True, , interval=1
# # anim.save('depth2lam1jperturbed.gif', writer='Pillow')
# plt.show()
