#import sys
import matplotlib.pyplot as plt
from scipy import stats as st
import numpy as np

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 

def Example11_1(mu_r,sig_r,mu_e,sig_e,g_a,g_b):
    deno    = sqrt( sig_r**2 + g_a**2*sig_e**2 )
    # sensitivity parameters
    alpha_r = sig_r/deno
    alpha_e = -g_a*sig_e/deno
    # reliability index
    beta    = (mu_r -g_a*mu_e - g_b)/deno
    pf      = st.norm.cdf(-beta)
    # design point
    r_d     = mu_r - beta*alpha_r*sig_r
    r_e     = mu_e - beta*alpha_e*sig_e
    #
    print(f" reliability index {beta:8.2f}\n failure probability {pf:8.3f}\n sensitivity parameter r {alpha_r:8.3f}\n sensitivity parameter e {alpha_e:8.3f}\n design point r_d  {r_d:8.3f}\n design point r_e  {r_e:8.3f}")
def Example11_3(mu_r,sig_r,mu_e,sig_e, beta,delta, alpha_r = 0.7, alpha_e = -0.8 ):
    deno    = sqrt( sig_r**2 + sig_e**2 )
    # coefficients of variation
    nu_r = sig_r / mu_r
    nu_e = sig_e / mu_e
    # failure probability
    p_f     = st.norm.cdf(-beta)
    # partial safety factors
    gamma_r = (1-delta*nu_r) / (1-beta*nu_r*alpha_r) 
    gamma_e = (1-beta*nu_e*alpha_e) / (1+delta*nu_e)
    # 
    print(f" failure probability {p_f:e}\n coefficient of variation r {nu_r:8.3f}\n coefficient of variation e {nu_e:8.3f}\n ratio standard deviations {sig_e/sig_r}\n sensitivity parameter r {alpha_r:8.3f}\n sensitivity parameter e {alpha_e:8.3f}\n partial safety factor r {gamma_r:8.3f}\n partial safety factor e {gamma_e:8.3f}")

# limit state classes
class LimStateLin():
    def __init__(self, a, b):
        self.a = a
        self.b = b
#    def Val(self, r ):
#        return (r-self.b)/self.a                                            # returns e!!!
    def Val2(self, r, e):                                                   
        return r - self.a*e-self.b                                          # returns g
    def Safe(self, r, e):                                                   # w is weighting factor for importance sampling
        g = (r-self.a*e-self.b)
        if g>0: return 0                                                    # safe
        else:   return 1.                                                   # failure 
#    def Grad(self, r, e ):
#        d = np.sqrt(1+self.a**2)
#        return 1/d, -self.a/d
    def re(self, e):                                                         # returns r
        return self.a * e
    def rde(self, e):                                                        # 1st derivative of r depending on e
        return self.a
#    def rd2e(self, e):                                                       # 2nd derivative of r depending on e
#        return 0.
class LimitStateNL():
    def __init__(self, a, e_u, M_u):
        self.a = a
        self.e_u = e_u
        self.M_u = M_u
#    def Val(self, e ):
#        x = e/self.e_u
#        return a*self.e_u * (x-x**2) + self.M_u*x**2                        # returns r!!!
    def Val2(self, r, e):
        x = e/self.e_u
        return r-a*self.e_u * (x-x**2) + self.M_u*x**2                      # returns g
    def Safe(self, r, e):
#        g = r - self.Val( e ) 
        g = r - self.re( e ) 
        if g>0: return 0                                                    # safe
        else:   return 1.                                                   # failure 
#    def Grad(self, r, e):
#        dr = 1.
#        x  = e/self.e_u 
#        de = ( -self.a*self.e_u*(1.-2.*x) -self.M_u*2.*x )/self.e_u
#        d  = np.sqrt(dr**2 + de**2)
#        return dr/d, de/d
    def re(self, e):
        x = e/self.e_u
        return a*self.e_u * (x-x**2) + self.M_u*x**2                        # returns r!!! 
    def rde(self, e):
        x  = e/self.e_u 
        de = ( -self.a*self.e_u*(1.-2.*x) -self.M_u*2.*x )/self.e_u
        return de
# polyline for limit state curve
def LimitStatePoly( G, xLow, xHigh, N ):
    l, e = np.linspace( xLow, xHigh, N), []    
    for i in l: 
        e += [G.Val(i)]
    return l, e
def LimitStatePolyE( G, eLow, eHigh, N ):
    e, r = np.linspace( eLow, eHigh, N), []    
    for i in e: 
        r += [G.re(i)]
    return r, e
# design point iteratively - obsolete
#def DesignPoint( muR, sigR, muE, sigE, G, r):
#    resi =1.
#    e   = G.Val( r )
#    while resi>1.0e-6:
#        re_ = np.array([ (r-muR)/sigR, (e-muE)/sigE ])                          # normalized limit point
#        rg, eg  = G.Grad( r, e)                                                 # gradient
#        gg_ = np.array([ rg*sigR, eg*sigE])                                     # normalized gradient
#        ggx = gg_/np.linalg.norm(gg_)                                           # unit length gradient
#        lamx= np.dot(re_,ggx)                                                   # distance
#        ddx = lamx*ggx                                                           # normalized design point
#        r   = muR + sigR*ddx[0]
#        e   = G.Val( r)
 #       resi= G.Val2( r, e)
 #       break
#    return muR + sigR*ddx[0], muE + sigE*ddx[1]                             # physical design point 
def DesignPoint_(muR, sigR, muE, sigE, G, d ):                              # determine design point by brute bisection to find zero for distance derivative -- works only for parabola like distance with single extremum
    def DisData( x ):                                                       # takes normalized action !!!
        RL  = (G.re( muE + sigE*x )-muR)/sigR                               # normalized r value
        RdL = G.rde( muE + sigE*x )*sigE/sigR                               # normalized derivative 
        DL  = np.sqrt(RL**2 + x**2)                                         # distance in normalized space
        return (RL*RdL+x)/DL                                                # derivative of distance in normalized space
    eL_, eH_, N = -d, d, 100
    for i in range(N):
        dDL = DisData( eL_ )
        dDH = DisData( eH_ )
        if dDL*dDH>=0.: raise NameError("something wrong with design point by bisection")
        e0_ = 0.5*(eL_+eH_)
        if DisData( e0_ )*dDH > 0.: eH_ = e0_
        else:                       eL_ = e0_
        if np.abs(eH_-eL_)<1.0e-6: break
    if i==N-1: raise NameError("found no design point by bisection") 
    ed = muE + sigE*e0_
    rd = G.re( ed )
    return rd, ed

def Example11_2( muR, sigR, muE, sigE, G, nlG ):
    # normal distributed random sample
    def RandomSampleNormal(mu, sig, N):
        r = np.random.rand(N)
        s = mu + sig* st.norm.ppf(r)
        return s
    # probability density for discrete points
    def ProDenseNormal( xLow, xHigh, mu, sig, N):
        l = np.linspace( xLow, xHigh, N)
        p = st.norm.pdf(l, loc=mu, scale =sig)
        return l, p
    def JointProbNormal( r, e, muR, sigR, muE, sigE ):
        return st.norm.pdf(r, loc=muR, scale =sigR) * st.norm.pdf(e, loc=muE, scale =sigE)
    # determine failure rate from samples
    def FailRate( rs, es, ww, G):
        ff, nS, nF, pp = [], 0, 0, 0
        for i, r in enumerate(rs):
            f = G.Safe( r, es[i] )
            if f==0: 
                nS += 1
            else:    
                nF += 1
                pp += ww[i]
            ff += [ f ]
        return nS, nF, pp/i, ff
    def EvalSampling( label, rs, es, muR, sigR, muE, sigE, ww, G, d, pp):
#        if Flag: rg, eg = LimitStatePoly( G,   rLow, rHigh, 5 )
#        else:    eg, rg = LimitStatePoly( G,   eLow, eHigh, 5 )
        rg_,eg_ =  LimitStatePolyE( G, eLow, eHigh, 5 )
        rd, ed = DesignPoint_( muR, sigR, muE, sigE, G, d )
        if np.abs((G.re(ed)-rd)/rd) > 1.0e-9: raise NameError("design point mismatch",np.abs((G.re(ed)-rd)/rd))
        pp.plot(rg_,eg_)
        pp.scatter([rd],[ed], s=100) #, c='red')
        pp.scatter([muR],[muE], s=70)
        pp.set_xlim( rLow, rHigh)
        pp.set_ylim( eLow, eHigh)                                                       
        nS, nF, pp, _ = FailRate( rs, es, ww, G)
        print(f"{label:s}\n design point {rd:f},{ed:f}\n {nF:d} failures and {nS:d} safe from {N:d} and {pp:f} failure probability ")
        return rd, ed
    
    N, d = 50, 3.5
    rLow, rHigh, eLow, eHigh = muR-d*sigR, muR+d*sigR, muE-d*sigE, muE+d*sigE 
    # monte carlo sampling
    label = "MC sampling"
    rs = RandomSampleNormal(muR, sigR, N)
    es = RandomSampleNormal(muE, sigE, N)
    P0, p0 = DefPlot( label )
    p0.scatter(rs,es, s=20)                       
    ww = [ 1. for i in range(N) ]
    rd, ed   = EvalSampling( label, rs, es, muR, sigR, muE, sigE, ww, G, d, p0)
    rdn, edn = EvalSampling( label, rs, es, muR, sigR, muE, sigE, ww, nlG, d, p0)
    # linear density data to plot
    N = 100
    rx, rr = ProDenseNormal( rLow, rHigh, muR, sigR, N)
    ex, ee = ProDenseNormal( eLow, eHigh, muE, sigE, N)
    P2, p2 = DefPlot('N')
    p2.plot(rx,rr) 
    p2.hist(rs, bins=10, facecolor='lightblue', alpha=0.5, histtype='stepfilled',density=True)
    p2.plot(ex,ee)
    p2.hist(es, bins=10, facecolor='lightgreen', alpha=0.5, histtype='stepfilled',density=True)

    N = 50
    # importance sampling
    label = "MC importance sampling"
    ri = RandomSampleNormal(rd, sigR, N)
    ei = RandomSampleNormal(ed, sigE, N)
    P1, p1 = DefPlot(label)                                                   # plot samples in 2D plane
    p1.scatter(ri,ei, s=20)
    ww = [ JointProbNormal( ri[i],ei[i],muR,sigR,muE,sigE )/JointProbNormal( ri[i],ei[i],rd,sigR,ed,sigE ) for i in range(N) ]
    rd, ed   = EvalSampling( label, ri, ei, muR, sigR, muE, sigE, ww, G, d, p1)
    # inmportance sampling with nl
    label = "MC importance sampling nl"
    ri = RandomSampleNormal(rdn, sigR, N)
    ei = RandomSampleNormal(edn, sigE, N)
    P2, p2 = DefPlot(label)                                                   # plot samples in 2D plane
    p2.scatter(ri,ei, s=20)
    ww = [ JointProbNormal( ri[i],ei[i],muR,sigR,muE,sigE )/JointProbNormal( ri[i],ei[i],rdn,sigR,edn,sigE ) for i in range(N) ]
    rdn, edn = EvalSampling( label, ri, ei, muR, sigR, muE, sigE, ww, nlG, d, p2)
    rg_,eg_ =  LimitStatePolyE( G, eLow, eHigh, 5 )
    p2.plot(rg_,eg_)
    # 
    plt.show()

if __name__ == '__main__':
    a = 2
    b = 0
    mu_r, sig_r = 0.27, 0.015
    mu_e, sig_e = 0.1, 0.02
    #
    case = 3
    # failure probability
    if case == 1:
        Example11_1(mu_r,sig_r,mu_e,sig_e,a,b)
    # partial safety factors
    elif case == 2:
        beta  = 3.8
        delta = 1.65
        Example11_3(mu_r,sig_r,mu_e,sig_e, beta,delta)
    # monte carlo
    elif case == 3:
        LinG = LimStateLin( a=a, b=b ) 
        nlG  = LimitStateNL(a=a, e_u=mu_e, M_u=0.23)
        Example11_2( muR=mu_r, sigR=sig_r, muE=mu_e, sigE=sig_e, G=LinG, nlG=nlG )
        
    print("finished")