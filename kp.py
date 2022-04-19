from grammian_forms import *
import matplotlib.pyplot as plt

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

lam = 1
s = 3
t = 10
image = np.real(integrand(X, Y, t, lam, s)) #np.real(Psi(X,Y,t, lam, s))

# scattering potential Psi
def Psi(x, y, t, lam, s):
    if s not in [0,1,2,3]:
        raise('choose s in {0,1,2,3}')
    sth_derivative = np.diff(np.exp(phi(x, y, t, lam)), s)
    # print(np.shape(np.exp(-phi(x, y, t, lam))))
    # print(np.shape(sth_derivative))
    schur_polynomial_ps =  np.exp(-phi(x, y, t, lam))@sth_derivative
    psi = np.exp(phi(x, y, t, lam))@schur_polynomial_ps
    return psi

def integrand(x,y,t,lam,s):
    psi = Psi(x,y,t,lam,s)
    int = np.abs(psi)**2
    return int


def cumsum(list):
    new_list=[]
    j=0
    for i in range(0,len(list)):
        j+=list[i]
        new_list.append(j)
    return new_list
out = np.array([cumsum(row) for row in image])
tau = lambda wavefxn: np.array([cumsum(row) for row in image])
# print(out)
q = np.exp(phi(X, Y, t, lam))@(np.exp(-phi(X, Y, t, lam))@np.diff(np.exp(phi(X, Y, t, lam)), s))
q1 = wavefxn(X,Y, t, lam, s)
plt.imshow(np.real(tau(q1)))
plt.show()
