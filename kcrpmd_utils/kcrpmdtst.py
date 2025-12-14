import numpy as np

class KcrpmdTst:

    def __init__(self, beta, a, b, c, ms, ws, s0, s1, eps, Kq, Vq):

        """
        effective solvent coordinate factored 2-state spin-boson model as defined in KC-RPMD paper

             | 0.5*ms*ws^2*(s-s0)^2          Kq(q)             |
        H =  |                                                 | + I*[Vb(s,x) + Vq(q)]
             |        Kq(q)         0.5*ms*ws^2*(s-s1)^2 + eps |
             

        Vb(s,x) = sumj{0.5*Mj*wj^2*(xj-cj*s/(Mj*wj^2))^2}

        Args:
            beta ( double ): boltzmann factor [ units: a.u. ]
            a ( double ): KC-RPMD kinetic constraint parameter a [ units: a.u. ]
            b ( double ): KC-RPMD heavy-side parameter b [ units: a.u. ]
            c ( double ): KC-RPMD kinetic constraint switching parameter c [ units: a.u. ]
            ms ( double ): s coordinate mass [ units: a.u. ]
            ws ( double ): s coordinate angular frequency [ units: a.u. ]
            s0 ( double ): V0 parabola center [ units: a.u. ]
            s1 ( double ): V1 parabola center [ units: a.u. ]
            eps ( double ): V1 parabola vertical shift [ units: a.u. ]
            Kq ( vectorized function ): q coordinate coupling [ units: a.u. ]
            Vq ( vectorized function ): q coordinate potential [ units: a.u. ]
        """
        # grid parameters for s (log points centered around sdagger)
        self.lin_pts = 10000
        self.lin_width = 8
        self.log_base = 100
        self.log_pts = 10001
        self.log_width = 0.02
        # grid parameters for q (careful with q_low)
        self.q_pts = 2500
        self.q_low = -1.2
        self.q_high = 3.0
        # grid parameters for y
        self.y_pts = 2500
        self.y_low = -1.6
        self.y_high = 1.6

        self.beta = beta
        self.a = a
        self.b = b
        self.c = c
        self.ms = ms 
        self.ws = ws 
        self.s0 = s0 
        self.s1 = s1 
        self.eps = eps
        self.Kq = Kq
        self.Vq = Vq

        self.ydag = 0.0
        self.sdag = 0.5 * (self.s0 + self.s1) - self.eps / (self.ms * self.ws**2 * (self.s0 - self.s1))
        self.lam = 0.5 * self.ms * self.ws**2 * (self.s0 - self.s1)**2
        self.V0s = lambda s: 0.5 * self.ms * self.ws**2 * (s - self.s0)**2
        self.V1s = lambda s: 0.5 * self.ms * self.ws**2 * (s - self.s1)**2 + self.eps

    def log_array(self, base, pts):
        log_ar = np.zeros(2 * pts - 1)
        a = np.log(np.linspace(base**(-1), 1., pts)) / np.log(base)
        b = -np.flip(a)
        log_ar[:pts] = a
        log_ar[pts:] = b[1:]
        return log_ar

    def y_array(self):
        y_ar = np.linspace(self.y_low, self.y_high, self.y_pts)
        if not np.any(y_ar == 0.):
            y_ar = np.sort(np.append(y_ar, 0.))
        return y_ar

    def s_array(self):
        s_ar = np.zeros(self.lin_pts + 2 * self.log_pts - 1)
        s_ar[:self.lin_pts] = self.sdag + 0.5 * self.lin_width * np.linspace(-1., 1., self.lin_pts)
        s_ar[self.lin_pts:] = self.sdag + 0.5 * self.log_width * self.log_array(self.log_base, self.log_pts)
        s_ar = np.sort(s_ar)
        return s_ar

    def q_array(self):
        q_ar = np.linspace(self.q_low, self.q_high, self.q_pts)
        return q_ar

    def Pq(self, q):
        q_ar = self.q_array()
        exp_arg = -self.beta * (self.Vq(q_ar))
        exp_shift = np.max(exp_arg) - 500.
        return np.exp(-self.beta * self.Vq(q) - exp_shift) / np.trapz(np.exp(exp_arg - exp_shift), q_ar)

    ############################################################
    ####### general potentials, free energies, and rates #######
    ############################################################
    def V0(self, s, q):
        S, Q = np.meshgrid(s, q)
        return self.V0s(S) + self.Vq(Q)

    def V1(self, s, q):
        S, Q = np.meshgrid(s, q)
        return self.V1s(S) + self.Vq(Q)

    def K(self, s, q):
        S, Q = np.meshgrid(s, q)
        return abs(self.Kq(Q))
    
    def Vg(self, s, q):
        return 0.5 * (self.V0(s, q) + self.V1(s, q)) - 0.5 * np.sqrt((self.V0(s, q) - self.V1(s, q))**2
                                                                     + 4 * self.K(s, q)**2)
    def Ve(self, s, q):
        return 0.5 * (self.V0(s, q) + self.V1(s, q)) + 0.5 * np.sqrt((self.V0(s, q) - self.V1(s, q))**2
                                                                      + 4 * self.K(s, q)**2)
    def Fgs(self, s):
        q_ar = self.q_array()
        exp_arg = -self.beta * self.Vg(s, q_ar)
        exp_shift = np.max(exp_arg, axis=0, keepdims=True) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), q_ar, axis = 0)) / self.beta - exp_shift[0,:] / self.beta

    def Fg(self):
        s_ar = self.s_array()
        exp_arg = -self.beta * self.Fgs(s_ar)
        exp_shift = np.max(exp_arg) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), s_ar)) / self.beta - exp_shift / self.beta

    def kGR(self):
        q_ar = self.q_array()
        Kq_ar = self.Kq(q_ar)
        Pq_ar = self.Pq(q_ar)
        kMT_ar = 2 * np.pi * self.Kq(q_ar)**2 * np.sqrt(self.beta / (4 * np.pi * self.lam)) * np.exp(-self.beta * (self.lam + self.eps)**2 / (4 * self.lam)) * Pq_ar
        return np.trapz(kMT_ar, q_ar)

    def kBO(self):
        s_ar = self.s_array(); s_ar = s_ar[:1+np.argwhere(s_ar == self.sdag)[0,0]]
        exp_arg = -self.beta * self.Fgs(s_ar)
        exp_shift = np.max(exp_arg) - 500.
        return 1 / np.sqrt(2 * np.pi * self.beta * self.ms) * np.exp(exp_arg[-1] - exp_shift) / np.trapz(np.exp(exp_arg - exp_shift), s_ar)

    def kIF(self, kBO, kBO_0):
        q_ar = self.q_array()
        Kq_ar = self.Kq(q_ar)
        Pq_ar = self.Pq(q_ar)
        kMT_ar = 2 * np.pi * self.Kq(q_ar)**2 * np.sqrt(self.beta / (4 * np.pi * self.lam)) * np.exp(-self.beta * (self.lam + self.eps)**2 / (4 * self.lam)) * Pq_ar
        return np.trapz(kMT_ar * kBO / (kMT_ar + kBO_0), q_ar)

    ############################################################
    ####### KC-RPMD potentials, free energies, and rates #######
    ############################################################
    def set_eta_my_gammay(self):
        q_ar = self.q_array()
        exp_arg = -self.beta * self.Vq(q_ar)
        exp_shift = np.max(exp_arg) - 500.
        Pq_ar = self.Pq(q_ar)
        Kq_ar = abs(self.Kq(q_ar))
        self.eta = 2 * np.pi * (np.trapz(Kq_ar * Pq_ar, q_ar) * np.trapz(Kq_ar**2 * Pq_ar, q_ar)
                                / np.trapz(Kq_ar**3 * Pq_ar, q_ar))
        self.my = self.beta**3 * self.eta**2 / ((2*np.pi)**3) * (np.trapz(Kq_ar**3 * Pq_ar, q_ar)
                                                                 / np.trapz(Kq_ar**2 * Pq_ar, q_ar))**2
        self.gammay = 0.5 * np.sqrt((1 + np.abs(-2 * np.log(np.sqrt(self.a / np.pi) * self.eta * self.beta**2)
                                                - 4 * np.trapz(np.log(Kq_ar) * Pq_ar, q_ar))) / (self.beta * self.my))
        return None

    def VKP(self, s, q):
        V0 = self.V0(s, q)
        V1 = self.V1(s, q)
        K = self.K(s, q)
        Vg = self.Vg(s, q)
        Ve = self.Ve(s, q)

        VKP = np.zeros_like(K)
        mask_ad = (self.beta * K > 1e-3)
        mask_nad = (self.beta * K <= 1e-3) & (self.beta * np.abs(V0 - V1) > 1e-7)
        mask_hole = (self.beta * K <= 1e-3) & (self.beta * np.abs(V0 - V1) <= 1e-7)

        VKP[mask_ad] = Vg[mask_ad] - np.log(1
                                            + np.exp(-self.beta * (Ve[mask_ad] - Vg[mask_ad]))
                                            - np.exp(-self.beta * (V0[mask_ad] - Vg[mask_ad]))
                                            - np.exp(-self.beta * (V1[mask_ad] - Vg[mask_ad]))) / self.beta

        VKP[mask_nad] = 0.5 * (V0[mask_nad] + V1[mask_nad]) - np.log((self.beta * K[mask_nad])**2
                                                                     * np.sinh(0.5 * self.beta
                                                                               * (V0[mask_nad] - V1[mask_nad]))
                                                                     / (0.5 * self.beta
                                                                        * (V0[mask_nad] - V1[mask_nad]))) / self.beta

        VKP[mask_hole] = 0.5 * (V0[mask_hole] + V1[mask_hole]) - np.log((self.beta * K[mask_hole])**2) / self.beta

        return VKP

    def w(self, s, q):
        return (self.V0(s, q) - self.V1(s, q)) / self.K(s, q)

    def A(self, s, q):
        K = self.K(s, q)

        A = np.zeros_like(K)
        mask_ad = (self.beta * K > 1.0)
        mask_nad = (self.beta * K <= 1.0)

        A[mask_ad] = (self.a * np.exp(-2 * self.c * (self.beta * K[mask_ad] - 1.0))
                      / (1 + np.exp(-2 * self.c * (self.beta * K[mask_ad] - 1.0))))

        A[mask_nad] = self.a / (1 + np.exp(2 * self.c * (self.beta * K[mask_nad] - 1.0)));

        return A

    def C(self, s, q):
        K = self.K(s, q)

        C = np.zeros_like(K)
        mask_ad = (self.beta * K > 1.0)
        mask_nad = (self.beta * K <= 1.0)

        C[mask_ad] = (1 + (self.eta * np.sqrt(self.a / np.pi) * np.exp(-self.c * (self.beta * K[mask_ad] - 1.0))
                           / np.sqrt(1 + np.exp(-2 * self.c * (self.beta * K[mask_ad] - 1.0))) - 1)
                      * np.exp(-2 * self.c * (self.beta * K[mask_ad] - 1.0))
                      / (1 + np.exp(-2 * self.c * (self.beta * K[mask_ad] - 1.0))))

        C[mask_nad] = (1 + (self.eta * np.sqrt(self.a / np.pi)
                            / np.sqrt(1 + np.exp(2 * self.c * (self.beta * K[mask_nad] - 1.0))) - 1)
                       / (1 + np.exp(2 * self.c * (self.beta * K[mask_nad] - 1.0))))

        return C

    def VKC(self, theta, s, q):
        if theta == -1:
            return self.V0(s, q)
        elif theta == 0:
            return self.VKP(s, q) - np.log(self.C(s, q)) / self.beta + self.A(s, q) * self.w(s, q)**2 / self.beta
        elif theta == 1:
            return self.V1(s, q)
        else:
            print('Invalid Theta!')
            return np.sqrt(-1)

    def Fsq(self, s, q):
        exp0_arg = -self.beta * self.VKC(-1, s, q)
        exp1_arg = -self.beta * self.VKC(1, s, q)
        expKP_arg = -self.beta * self.VKC(0, s, q)
        exp_shift = np.maximum.reduce([exp0_arg, exp1_arg, expKP_arg]) - 500.
        return -np.log(np.exp(exp0_arg - exp_shift) + np.exp(exp1_arg - exp_shift) + np.exp(expKP_arg - exp_shift)) / self.beta - exp_shift / self.beta

    def Fs(self, s):
        q_ar = self.q_array()
        exp_arg = -self.beta * self.Fsq(s, q_ar)
        exp_shift = np.max(exp_arg, axis=0, keepdims=True) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), q_ar, axis = 0)) / self.beta - exp_shift[0,:] / self.beta

    def Fq(self, q):
        s_ar = self.s_array()
        exp_arg = -self.beta * self.Fsq(s_ar, q)
        exp_shift = np.max(exp_arg, axis=1, keepdims=True) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), s_ar, axis = 1)) / self.beta - exp_shift[:,0] / self.beta

    def Ftheta(self, theta):
        s_ar = self.s_array(); q_ar = self.q_array()
        exp_arg = -self.beta * self.VKC(theta, s_ar, q_ar)
        exp_shift = np.max(exp_arg) - 500.
        return -np.log(np.trapz(np.trapz(np.exp(exp_arg - exp_shift), s_ar, axis = 1), q_ar)) / self.beta - exp_shift / self.beta

    def Fthetas(self, theta, s):
        q_ar = self.q_array()
        exp_arg = -self.beta * self.VKC(theta, s, q_ar)
        exp_shift = np.max(exp_arg, axis=0, keepdims=True) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), q_ar, axis = 0)) / self.beta - exp_shift[0,:] / self.beta

    def Fthetaq(self, theta, q):
        s_ar = self.s_array()
        exp_arg = -self.beta * self.VKC(theta, s_ar, q)
        exp_shift = np.max(exp_arg, axis=1, keepdims=True) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), s_ar, axis = 1)) / self.beta - exp_shift[:,0] / self.beta

    def Vr(self, theta, y):
        return np.piecewise(y, [np.abs(y - theta) < 0.5, np.abs(y - theta) >= 0.5],
                            [lambda y: -np.log(1 / (1 + np.exp(self.b * (2 * np.abs(y - theta) - 1)))) / self.beta,
                             lambda y: (self.b * (2 * np.abs(y - theta) - 1) / self.beta
                                        - np.log(1 / (1 + np.exp(-self.b * (2 * np.abs(y - theta) - 1)))) / self.beta)])

    def Fy(self, y):
        exp0_arg = -self.beta * (self.Ftheta(-1) + self.Vr(-1, y))
        exp1_arg = -self.beta * (self.Ftheta(1) + self.Vr(1, y))
        expKP_arg = -self.beta * (self.Ftheta(0) + self.Vr(0, y))
        exp_shift = np.maximum.reduce([exp0_arg, exp1_arg, expKP_arg]) - 500.
        return -np.log(np.exp(exp0_arg - exp_shift) + np.exp(exp1_arg - exp_shift) + np.exp(expKP_arg - exp_shift)) / self.beta - exp_shift / self.beta

    def Fys(self, y, s):
        exp0_arg = -self.beta * (self.Fthetas(-1, s)[:,np.newaxis] + self.Vr(-1, y)[np.newaxis,:])
        exp1_arg = -self.beta * (self.Fthetas(1, s)[:,np.newaxis] + self.Vr(1, y)[np.newaxis,:])
        expKP_arg = -self.beta * (self.Fthetas(0, s)[:,np.newaxis] + self.Vr(0, y)[np.newaxis,:])
        exp_shift = np.maximum.reduce([exp0_arg, exp1_arg, expKP_arg]) - 500.
        return -np.log(np.exp(exp0_arg - exp_shift) + np.exp(exp1_arg - exp_shift) + np.exp(expKP_arg - exp_shift)) / self.beta - exp_shift / self.beta

    def Fyq(self, y, q):
        exp0_arg = -self.beta * (self.Fthetaq(-1, q)[:,np.newaxis] + self.Vr(-1, y)[np.newaxis,:])
        exp1_arg = -self.beta * (self.Fthetaq(1, q)[:,np.newaxis] + self.Vr(1, y)[np.newaxis,:])
        expKP_arg = -self.beta * (self.Fthetaq(0, q)[:,np.newaxis] + self.Vr(0, y)[np.newaxis,:])
        exp_shift = np.maximum.reduce([exp0_arg, exp1_arg, expKP_arg]) - 500.
        return -np.log(np.exp(exp0_arg - exp_shift) + np.exp(exp1_arg - exp_shift) + np.exp(expKP_arg - exp_shift)) / self.beta - exp_shift / self.beta

    def F(self):
        y_ar = self.y_array()
        exp_arg = -self.beta * self.Fy(y_ar)
        exp_shift = np.max(exp_arg) - 500.
        return -np.log(np.trapz(np.exp(exp_arg - exp_shift), y_ar)) / self.beta - exp_shift / self.beta

    def tst_y(self):
        y_ar = self.y_array(); y_ar = y_ar[:1+np.argwhere(y_ar == 0.)[0,0]]
        exp_arg = -self.beta * self.Fy(y_ar)
        exp_shift = np.max(exp_arg) - 500.
        return 1 / np.sqrt(2 * np.pi * self.beta * self.my) * np.exp(exp_arg[-1] - exp_shift) / np.trapz(np.exp(exp_arg - exp_shift), y_ar)

    def tst_s(self):
        s_ar = self.s_array(); s_ar = s_ar[:1+np.argwhere(s_ar == self.sdag)[0,0]]
        exp_arg = -self.beta * self.Fs(s_ar)
        exp_shift = np.max(exp_arg) - 500.
        return 1 / np.sqrt(2 * np.pi * self.beta * self.ms) * np.exp(exp_arg[-1] - exp_shift) / np.trapz(np.exp(exp_arg - exp_shift), s_ar)
 
