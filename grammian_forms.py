import numpy as np
import scipy as sp
import scipy.special as special

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

# scattering potential Psi
def Psi(x, y, t, lam, s):
    if s not in [0,1,2,3]:
        raise('choose s in {0,1,2,3}')
    sth_derivative = np.diff(np.exp(phi(x, y, t, lam)), s)
    schur_polynomial_ps =  np.exp(-phi(x, y, t, lam))[s:-s][s:-s]*sth_derivative
    psi = schur_polynomial_ps*np.exp(phi(x, y, t, lam))



# def integrater():
#     def integrand(x, n=0, a=1):
#         return np.power(x, n)*np.exp(a*x)
    # result = sp.integrate.quad(integrand, -np.inf, np.inf, args=(n, a))
    # result = special.expn(-n, -x)
# def get_grammian(x, y, t, lam):
#     def integrand(x, y, t, lam, pi):
#         return np.abs(Psi(x, y, t, lam, pi))
#     def tau:
#         return sp.integrate.quad(integrand, -np.inf, z, args=(x, y, t, lam, pi) )

# depth 0 soliton
def tau_depth0(x, y, t, lam):
    return np.exp(2*(lam*x - 4*lam**3*t))/(2*lam)

# depth 1 lump
def tau_til_depth1(x, y, t, lam):
    # polynomial factor of tau
    return np.power( x - 12*lam**2*t - 1/lam, 2 ) + 4*lam**2*np.power(y , 2) + 1/(4*lam**2)

# depth 1 lump
def tau_depth1(x, y, t, lam):
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

# depth 3 lump
def tau_til_depth3(x, y, t, lam):
    # coefficients of modulus of p_3
    alp = -36*lam*(lam+2)*t
    bet = 4*lam*(lam+3)*(np.power(y, 2) + 72*lam**2*t**2) + 144*lam**3*(lam+6)*t**2 - 20*lam**2*y**2
    gam = -12*lam**2*t*(4*lam*(lam+3)*np.power(y, 2)+144*lam**3*(lam+6)*t**2) + 48*lam**3*(2*lam+9)*t*np.power(y, 2)
    det = 12*lam*y
    eps = -72*lam**2*(3*lam+5)*t*y
    zet = 8*lam**2*(lam+3)*np.power(y, 2) + 144*lam**3*(lam+6)*t**2-72+288*lam**4*(2*lam+9)*t**2*y

    # coefficients of tau
    c6 = 1
    c5 = -3/lam + 2*alp
    c4 = 30/(2*lam)**2 - 2*alp*5/(2*lam)+ alp**2 + 2*bet + det**2
    c3 = -120/(2*lam)**3+2*alp*20/(2*lam)**2 - (alp**2 + 2*bet + det**2)*4/(2*lam) + 2*(gam+alp*bet+det*eps)
    c2 = 180/(2*lam)**4 - 2*alp*60/(2*lam)**3 + (alp**2 + 2*bet + det**2)*12/(2*lam)**2
    c2 += - 2*(gam+alp*bet+det*eps)*3/(2*lam) + (bet + 2*alp*gam + eps**2 + 2*det*zet)
    c1 = -720/(2*lam)**5 +2*alp*120/(2*lam)**4 - (alp**2 + 2*bet + det**2)*24/(2*lam)**3
    c1 += 2*(gam+alp*bet+det*eps)*6/(2*lam)**2 - (bet + 2*alp*gam + eps**2 + 2*det*zet)/lam + 2*(bet*gam+eps*zet)
    c0 = 720/(2*lam)**6 - 2*alp*120/(2*lam)**5 + (alp**2 + 2*bet + det**2)*24/(2*lam)**4
    c0 += - 2*(gam+alp*bet+det*eps)*6/(2*lam)**3 + (bet + 2*alp*gam + eps**2 + 2*det*zet)*2/(2*lam)**2 - 2*(bet*gam+eps*zet)/(2*lam) + gam**2 + zet**2

    tau = c6*np.power(x, 6) + c5*np.power(x, 5) + c4*np.power(x, 4) + c3*np.power(x, 3)
    tau += c2*np.power(x, 2) + c1*x + c0
    return tau


def depth0_scattering(x, y, t, lam, eta=1):
    c = 1/np.sqrt(lam*eta)
    PHI = phi(x, y, t, lam)
    Rphi = np.real(PHI)
    tau_til = c*np.cosh(Rphi + .5*np.log(lam/eta)) + 2*np.cos((lam**2-eta**2)*y)/(lam+eta)
    return tau_til

def depth1_scattering(x, y, t, lam, eta=1):
    out = tau_til_depth1(x, y, t, lam) + tau_til_depth1(x, y, t, eta)
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

def tau_til_depth0_hr(x, y, t, lam, eta=1):
    return ((lam-eta)/(lam+eta))**2/(4*lam*eta)*np.exp(2*np.real(phi(x, y, t, lam+eta)))

# Rank 2 depth 1 soliton solutions - hr is higher rank
def tau_til_depth1_hr(x, y, t, lam, eta=1):
    out = tau_til_depth1(x, y, t, lam)*tau_til_depth1(x, y, t, eta)
    # t = t.astype('complex128')
    # add a mixture term
    alpha = 1/(lam+eta)
    Rphi = np.real(phi(x, y, t, lam+eta))
    # Rphi = Rphi.astype('complex128')
    Rt = np.power(x, 2) - (2*alpha + 12*(lam**2+eta**2)*t)*x
    # Rt = Rt.astype('complex128')
    Rt += 2*alpha**2 + 4*eta*lam*np.power(y, 2)
    Rt += 12*(lam**2 + eta**2)*t*alpha + (12*lam*eta*t)**2
    It = 2*(lam - eta)*y*(x + 12*eta*lam*t - alpha)
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
            'depth 1 alt': tau_depth1,
            'depth 2': tau_til_depth2,
            'depth 3': tau_til_depth3,
            'scattering 0': depth0_scattering,
            'scattering 1': depth1_scattering,
            'scattering 01': depth01_scattering,
            'rank 2 depth 0': tau_til_depth0_hr,
            'rank 2 depth 1': tau_til_depth1_hr
        }
