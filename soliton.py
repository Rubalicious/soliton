import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button, TextBox, RadioButtons
from grammian_forms import *
# plt.rcParams['text.usetex'] = True

# class Soliton(object):
#     """docstring for Soliton."""
#
#     def __init__(self, arg):
#         super(Soliton, self).__init__()
#         self.arg = arg

# Define domain parameters
N = 200
window = 100
# Lx = 5
# Ly = 100
x = np.linspace(-window, window, N)
y = np.linspace(-window, window, N)
X,Y = np.meshgrid(x,y)
X = X.astype('complex128')
Y = Y.astype('complex128')
# print(np.shape(X))

# define eigenvalues
eps = .1 # perturbation
lam = 1 #+eps # +eps*1j #np.exp(1j*np.pi*(1/4))
eta = 1

def build_solution(tau, t=0, lam=1, eta=1):
    TAU = tau(X, Y, t, lam) # , eta
    logTAU = np.log(TAU)
    U = 2*np.diff(logTAU, n=2, axis=0)
    return U

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
    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=plt.figaspect(0.5))
    U = np.real(build_solution(tau, t=0, lam=lam))
    # U = U.astype('float64')
    # =============
    # First subplot
    # =============
    # set up the axes for the first plot
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X.astype('float64')[1:-1], Y.astype('float64')[1:-1], U.astype('float64'),
                            rstride=1, cstride=1, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False
                            )
    fig.colorbar(surf, shrink=0.5, aspect=10)

    # ==============
    # Second subplot
    # ==============
    # set up the axes for the second plot
    ax = fig.add_subplot(1, 2, 2)
    ax.imshow(U)
    plt.show()

def my_first_widget():
    global X, Y

    # notes:
    # add 2D and 3D option perhaps even an expansion of window to view 2D and 3D
    # add resolution variable N
    # add a save image button
    fig = plt.figure(figsize=(10,8))
    ax = plt.axes()
    print(ax.dataLim)

    tau = forms['depth 1']
    U = np.real(build_solution(tau))

    # Rphi = np.real(phi(X, Y, 0, lam))
    img = plt.imshow(U)

    # to get rid of axis labels
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    # plt.colorbar()
    # surf = ax.plot_wireframe(X[1:-1], Y[1:-1], U, linewidth=.2)
    # convert inputs of plot_surface to float64 to use cmap functionality
    # U = U.astype('float64')
    # X = X.astype('float64')
    # Y = Y.astype('float64')
    # surf = ax.plot_surface(X[1:-1].astype('float64'), Y[1:-1].astype('float64'), U.astype('float64'), linewidth=.2, cmap='plasma') # , cmap='coolwarm'

    plt.subplots_adjust(bottom=0.25)

    # magnitude of eigenvalue slider
    axmag = plt.axes([0.2, 0.2, 0.65, 0.03])
    mag_slider = Slider(
        ax=axmag,
        label='magnitude',
        valmin=0,
        valmax=2,
        valinit=1,
    )

    # phase of eigenvalue slider
    axeigen = plt.axes([0.2, 0.15, 0.65, 0.03])
    eigen_slider = Slider(
        ax=axeigen,
        label='phase',
        valmin=0,
        valmax=2*np.pi,
        valinit=1,
    )

    # time axis slider
    axtime = plt.axes([0.2, 0.1, 0.65, 0.03])
    time_slider = Slider(
        ax=axtime,
        label='time',
        valmin=-10,
        valmax=10,
        valinit=0,
    )



    # grammian forms to choose from
    axradio = plt.axes([0.05, .58, .2, .3])
    rb = RadioButtons(
        ax=axradio,
        labels=forms.keys(),
        active=1,
        activecolor='blue'
    )

    # under development
    axtoggler = plt.axes([0.8, .8, .1, .1])
    toggle = RadioButtons(
        ax=axtoggler,
        labels=['2D', '3D'],
        active=0,
        activecolor='blue'
    )

    # Domain adjustment
    axtext1 = plt.axes([0.05+0.03, 0.5, 0.05, 0.03])
    xmin = TextBox(axtext1, label='xmin', initial=-100, color='.95', hovercolor='1', label_pad=0.1, textalignment='left')
    axtext1 = plt.axes([0.05+0.03, 0.45, 0.05, 0.03])
    ymin = TextBox(axtext1, label='ymin', initial=-100, color='.95', hovercolor='1', label_pad=0.1, textalignment='left')
    axtext1 = plt.axes([0.17+0.03, 0.5, 0.05, 0.03])
    xmax = TextBox(axtext1, label='xmax', initial=100, color='.95', hovercolor='1', label_pad=0.1, textalignment='left')
    axtext1 = plt.axes([0.17+0.03, 0.45, 0.05, 0.03])
    ymax = TextBox(axtext1, label='ymax', initial=100, color='.95', hovercolor='1', label_pad=0.1, textalignment='left')

    # Resolution adjustment
    axresolution = plt.axes([0.18, 0.38, 0.05, 0.03])
    res = TextBox(axresolution, label='Number of grid points', initial=N, color='.95', hovercolor='1', label_pad=0.1, textalignment='left')

    # eigenvalue subplot
    a = plt.axes([.77, .27, .1, .1])
    a.set_title("phase")
    circle = np.linspace(0, 2*np.pi, 100)
    a.plot(np.cos(circle), np.sin(circle))
    a.plot(np.real(lam), np.imag(lam), 'x')
    a.set_xticks([])
    a.set_yticks([])

    # save button
    axsave = plt.axes([.77, .5, .08, .05])
    axsave.set_xticks([])
    axsave.set_yticks([])
    save_button = Button(axsave, 'Save')
    # Turn off tick labels


    def save_image(event):
        print('event:{}'.format(event))
        folder = './figures/'
        filename = '{}_time={}_lam={}.png'.format(rb.value_selected.replace(' ',''), time_slider.val, np.round(np.exp(1j*eigen_slider.val), 2) )
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        # superimpose phase
        fig.savefig(folder+filename, bbox_inches=extent)
        # ax.savefig(folder+filename)
        print('image saved : {}'.format(folder+filename))

    def update(val):
        # print(ymin.text)
        tau = forms[rb.value_selected]
        eigenvalue = mag_slider.val*np.exp(1j*eigen_slider.val)
        U = np.real(build_solution(tau, t=time_slider.val, lam=eigenvalue )  )
        # Rphi = np.real(phi(X, Y, time_slider.val, eigen_slider.val))
        # img.set_data(U)
        ax.clear()

        if toggle.value_selected == '3D':
            print(ax.projection)
            ax.plot_wireframe(X[1:-1], Y[1:-1], U, linewidth=.2)
            ax.plot_surface(X[1:-1].astype('float64'), Y[1:-1].astype('float64'), U.astype('float64'), linewidth=.2, cmap='plasma')
        elif toggle.value_selected == '2D':
            # ax.axes(projection='2d')
            ax.imshow(U)
            ax.set_xticks([])
            ax.set_yticks([])
        # update phase picture
        a.clear()
        a.plot(np.cos(circle), np.sin(circle))
        a.plot(np.real(eigenvalue), np.imag(eigenvalue), 'x')
        a.set_xticks([])
        a.set_yticks([])

        fig.canvas.draw_idle()

    def update_domain(text):
        try:
            global x, y, X, Y
            x = np.linspace(int(xmin.text), int(xmax.text), int(res.text))
            y = np.linspace(int(ymin.text), int(ymax.text), int(res.text))
            X,Y = np.meshgrid(x,y)
            X = X.astype('complex128')
            Y = Y.astype('complex128')
            # tau = forms[rb.value_selected]
            # U = np.real(build_solution(tau, t=time_slider.val, lam=np.exp(1j*eigen_slider.val)))
            # img.set_data(U)
            fig.canvas.draw_idle()
            print('values have changed')
        except:
            raise('make sure entries are number')

    def figures_for_paper(event):
        tau = forms[rb.value_selected]
        eigenvalue = mag_slider.val*np.exp(1j*eigen_slider.val)
        gig = plt.figure(figsize=(8, 3))
        # gig.suptitle(rb.value_selected)
        n = 3
        dT = 3
        for i in range(n):
            ax = gig.add_subplot(1,n,i+1)
            T = -dT+i*dT
            U = np.real(build_solution(tau, t=T, lam=eigenvalue )  )
            ax.imshow(U)
            ax.set_title('t = {}'.format(T))
            ax.set_xticks([])
            ax.set_yticks([])
        # ax = gig.add_subplot(132)
        # U = np.real(build_solution(tau, t=0, lam=eigenvalue )  )
        # ax.imshow(U)
        # ax.set_title('t = 0')
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax = gig.add_subplot(133)
        # U = np.real(build_solution(tau, t=5, lam=eigenvalue )  )
        # ax.imshow(U)
        # ax.set_title('t = 5')
        # ax.set_xticks([])
        # ax.set_yticks([])

        gig.savefig('./figures/{}.png'.format(rb.value_selected))


    def figures_for_paper2(event):
        tau = forms[rb.value_selected]
        eigenvalue = mag_slider.val*np.exp(1j*eigen_slider.val)
        gig = plt.figure(figsize=(10, 5))
        ax = gig.add_subplot(151)
        U = np.real(build_solution(tau, t=-5, lam=eigenvalue )  )
        ax.imshow(U)
        ax.set_title('t = -5')
        ax = gig.add_subplot(152)
        U = np.real(build_solution(tau, t=0, lam=eigenvalue )  )
        ax.imshow(U)
        ax.set_title('t = 0')
        ax = gig.add_subplot(153)
        U = np.real(build_solution(tau, t=5, lam=eigenvalue )  )
        ax.imshow(U)
        ax.set_title('t = 5')
        gig.savefig('./figures/{}.png'.format(tau))
        ax = gig.add_subplot(154)
        U = np.real(build_solution(tau, t=5, lam=eigenvalue )  )
        ax.imshow(U)
        ax.set_title('t = 5')
        gig.savefig('./figures/{}.png'.format(tau))
        ax = gig.add_subplot(155)
        U = np.real(build_solution(tau, t=5, lam=eigenvalue )  )
        ax.imshow(U)
        ax.set_title('t = 5')
        gig.savefig('./figures/{}.png'.format(tau))
        # gig.close()
    def figures_for_paper3(event):
        tau = forms[rb.value_selected]

        gig = plt.figure(figsize=(8, 3))
        # gig.suptitle(rb.value_selected)
        n = 5
        eigenvalues = .5*np.exp(1j*np.linspace(0, np.pi/2, n))
        for i in range(n):
            ax = gig.add_subplot(1,n,i+1)
            T = 5
            U = np.real( build_solution(tau, t=T, lam=eigenvalues[i])  )
            ax.imshow(U)
            ax.set_title('({})'.format(i+1))
            ax.set_xticks([])
            ax.set_yticks([])
        gig.savefig('./figures/{}_eigenvectors1.png'.format(rb.value_selected))

    time_slider.on_changed(update)
    eigen_slider.on_changed(update)
    mag_slider.on_changed(update)
    rb.on_clicked(update)
    xmin.on_submit(update_domain)
    xmax.on_submit(update_domain)
    ymin.on_submit(update_domain)
    ymax.on_submit(update_domain)
    save_button.on_clicked(figures_for_paper)
    # toggle.on_clicked(update)

    plt.show()

def make_movie(tau):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    lam1 = 1+1j
     #  projection='3d'
    def animate(i):
        # phase = (i/60)*(np.pi*1j/4) + 3*np.pi/4
        # phase = (i/60)*2*np.pi*1j
        # lam1 = np.exp(phase)
        t = -25+i
        # lam1 = 1+1j
        Z = np.real(build_solution( tau , t=t, lam=lam1, eta=eta ))
        ax.clear()
        # ax.set_zlim(-5, 30)
        plt.title(r't={}, $\lambda$={},\eta={}$'.format(t, lam1, eta))
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

    anim = FuncAnimation(fig, animate, frames=50, interval=1) # , init_func=init, blit=True, , interval=1
    # anim.save('depth1_scattering_soliton1_lam={}_eta={}.gif'.format(lam1, eta), writer='Pillow')
    plt.show()


# PHI = phi(X, Y, 0, lam)
# P1 = p1(X, Y, 0, lam)
# P2 = p2(X, Y, 0, lam)
# P3 = p3(X, Y, 0, lam)

# plot_heatmap(tau_til_depth2)
my_first_widget()
# make_movie(tau_til_depth2)

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

# Rphi = np.real(phi(X, Y, 0, lam+eta))
# # Z = np.real(build_solution( tau_til_depth1_hr , t=0, lam=np.exp(-1j*2*np.pi*0.1), eta=eta ))
# plt.imshow(Rphi)
# plt.show()
