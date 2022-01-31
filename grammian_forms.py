import numpy as np

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

# depth 0 soliton
def tau_depth0(x, y, t, lam):
    return np.exp(2*(lam*x - 4*lam**3*t))/(2*lam)

# depth 1 lump
def tau_til_depth1(x, y, t, lam):
    # polynomial factor of tau
    return np.power( x - 12*lam**2*t - 1/lam, 2 ) + 4*lam**2*np.power(y , 2) + 1/(4*lam**2)

# depth 1 lump
def tau(x, y, t, lam):
    alpha = 1/(2*lam)
    PHI = phi(x, y, t, lam)
    Rphi = np.real(PHI)
    phase = np.exp(2*Rphi)
    polynomial = tau_til_depth1(x, y, t, lam)
    return alpha*phase*polynomial

# depth 2 lump
def tau_til_depth2(x, y, t, lam):
    k = y-6*lam*t
    B = 4*lam*k
    C = 4*lam*(lam*np.power(k, 2)-6*t)
    B2C = np.power(B, 2) + 2*C

    fac = B - 1/lam

    alpha = 2*fac
    beta = (-3/lam)*fac + B2C
    gamma = (3/lam**2)*fac - (1/lam)*B2C + 2*B*C
    delta = (-3/(2*lam**3))*fac + 2*B2C/(2*lam)**2 - B*C/lam + np.power(C, 2) + 4*np.power(y, 2)
    quartic = np.power(x, 4) + alpha*np.power(x, 3) + beta*np.power(x, 2) + gamma*x + delta
    return quartic

def depth0_scattering(x, y, t, lam, eta=1):
    c = 1/np.sqrt(lam*eta)
    PHI = phi(x, y, t, lam)
    Rphi = np.real(PHI)
    tau_til = c*np.cosh(Rphi + .5*np.log(lam/eta)) + 2*np.cos((lam**2-eta**2)*y)/(lam+eta)
    return tau_til

def depth1_scattering(x, y, t, lam, eta=1):
    out = tau(x, y, t, lam) + tau(x, y, t, eta)
    alpha = 1/(lam+eta)
    tau_le = np.power(x, 2) - (2*alpha - 12*(lam**2+eta**2)*t)*x + 2*alpha**2 + 12*(lam**2+eta**2)*alpha*t + (12*lam*eta*t)**2+4*eta*lam*np.power(y, 2)
    tau_le *= np.cos((lam**2-eta**2)*y)
    tau_le += 2*(lam-eta)*y*np.sin((lam**2-eta**2)*y)*(x-alpha+12*lam*eta*t)
    out+= tau_le
    return out

def depth01_scattering(x, y, t, lam, eta=1):
    z = x - 12*eta**2*t
    Rphi = np.real(phi(x, y, t, lam+eta))
    fac = np.exp(Rphi)
    arg = np.power(z-1/eta, 2) + 1/(4*eta**2) + 4*eta**2*np.power(y, 2)
    alpha = np.sqrt(arg/(lam*eta))
    trig1 = np.cosh(np.real(phi(x, y, t, lam-eta)) - np.log(arg)/2 - np.log(lam/eta))
    term2 = ((z-1/(lam+eta))*np.cos((eta**2-lam**2)*y) - 2*eta*y*np.sin((eta**2-lam**2)*y))/(lam+eta)
    out = fac*(alpha*trig1+term2)
    return out

# Rank 2 depth 1 soliton solutions - hr is higher rank
def tau_til_depth1_hr(x, y, t, lam, eta=1):
    out = tau(x, y, t, lam)*tau(x, y, t, eta)
    # t = t.astype('complex128')
    # add a mixture term
    alpha = 1/np.complex(lam+eta)
    Rphi = np.real(phi(x, y, t, lam+eta))
    # Rphi = Rphi.astype('complex128')
    Rt = np.power(x, 2) - (2/(lam+eta) + 12*(lam**2+eta**2)*t)*x
    # Rt = Rt.astype('complex128')
    Rt += 2*alpha**2 + 4*eta*lam*np.power(y, 2)
    Rt += 12*(lam**2+eta**2)*t*alpha + (12*lam*eta*t)**2
    It = 2*(lam-eta)*y*(x+12*eta*lam*t-alpha)
    mix = np.power(Rt, 2)
    mix += np.power(It, 2)
    phase = np.exp(2*Rphi)
    mix *= phase
    mix /= alpha**2
    out -= mix
    return out

forms = {
            'depth 0': tau_depth0,
            'depth 1': tau_til_depth1,
            'depth 2': tau_til_depth2,
            'scattering 0': depth0_scattering,
            'scattering 1': depth1_scattering,
            'scattering 01': depth01_scattering
        }
