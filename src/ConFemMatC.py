# This file was automatically generated by SWIG (https://www.swig.org).
# Version 4.1.1
#
# Do not make changes to this file unless you know what you are doing - modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _ConFemMatC
else:
    import _ConFemMatC

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "this":
            set(self, name, value)
        elif name == "thisown":
            self.this.own(value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)



def PrinCLT_1(vx, vy, vxy, pEps):
    return _ConFemMatC.PrinCLT_1(vx, vy, vxy, pEps)

def PrinCLT_2(vx, vy, vxy):
    return _ConFemMatC.PrinCLT_2(vx, vy, vxy)

def ElasticLTC2(ElemDim, ElemPlSt, ElemStateVar, ElemStateVarN, Eps, sig, MatM, ww, nu, Emod, Dps, Dt, rho, selfwcr_, selfbw_, selfepsct, ElemLch_, selfCrackVisc, selffct, selfepscu, DataOut):
    return _ConFemMatC.ElasticLTC2(ElemDim, ElemPlSt, ElemStateVar, ElemStateVarN, Eps, sig, MatM, ww, nu, Emod, Dps, Dt, rho, selfwcr_, selfbw_, selfepsct, ElemLch_, selfCrackVisc, selffct, selfepscu, DataOut)

def eighC(I1, J2, Eps, la, vv):
    return _ConFemMatC.eighC(I1, J2, Eps, la, vv)

def EquivStrain1C(CalcType, I1, J2, J2S, Eps, EpsD, kap_, nd, Hd, cc0, cc1, cc2, cc3, LaMax, eVec, data):
    return _ConFemMatC.EquivStrain1C(CalcType, I1, J2, J2S, Eps, EpsD, kap_, nd, Hd, cc0, cc1, cc2, cc3, LaMax, eVec, data)

def eighC2(I1, J2D, Eps, la, vv):
    return _ConFemMatC.eighC2(I1, J2D, Eps, la, vv)

def eigJacobiSym(Eps, la, vv):
    return _ConFemMatC.eigJacobiSym(Eps, la, vv)

def ViscExten3DC1(Dt, eta, Dps, ElemStateVar, ElemStateVarN, sI, Veps):
    return _ConFemMatC.ViscExten3DC1(Dt, eta, Dps, ElemStateVar, ElemStateVarN, sI, Veps)

def eigJacobiSymWrapper(Eps, la, vv):
    return _ConFemMatC.eigJacobiSymWrapper(Eps, la, vv)

def IsoDamC1(CalcType, ElemDim, ElemPlSt, PrinStrains, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM2, LiTy, cc0, cc1, cc2, cc3, RType, EpsR, kapStrength, ElemCrBws, gam2, kapUlt, edt, ed, gd, nu, Emod, Dps, eta, RegPar, ElemScaleType, sigR, CR, Dt, sVsTol, DataOut):
    return _ConFemMatC.IsoDamC1(CalcType, ElemDim, ElemPlSt, PrinStrains, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM2, LiTy, cc0, cc1, cc2, cc3, RType, EpsR, kapStrength, ElemCrBws, gam2, kapUlt, edt, ed, gd, nu, Emod, Dps, eta, RegPar, ElemScaleType, sigR, CR, Dt, sVsTol, DataOut)

def IsoDamUniaxC1(ElemDim, ElemPlSt, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM, LiTy, alpha, cc1, cc2, cc3, RType, EpsR, kapStrength, ElemCrBws, gam2, kapUlt, edt, ed, gd, nu, Emod, Dps, eta, RegPar, sigR, CR, Dt, sVsTol, DataOut):
    return _ConFemMatC.IsoDamUniaxC1(ElemDim, ElemPlSt, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM, LiTy, alpha, cc1, cc2, cc3, RType, EpsR, kapStrength, ElemCrBws, gam2, kapUlt, edt, ed, gd, nu, Emod, Dps, eta, RegPar, sigR, CR, Dt, sVsTol, DataOut)

def Rankine(PrinStrains, ElemDim, Eps_, sig0, eVec, laMax):
    return _ConFemMatC.Rankine(PrinStrains, ElemDim, Eps_, sig0, eVec, laMax)

def RankineUpd(sig0, eVec, ElemStateVarN, laMax, off):
    return _ConFemMatC.RankineUpd(sig0, eVec, ElemStateVarN, laMax, off)

def DamFunc1(kap0V, alphaV, betaV, kapOld, eta, kap, dd, Dp):
    return _ConFemMatC.DamFunc1(kap0V, alphaV, betaV, kapOld, eta, kap, dd, Dp)

def DevStiffness(m00, m01, m02, m11, m12, m22, n02, n12, n22, tbt, tbt2, obt, obt2, nn, DVD):
    return _ConFemMatC.DevStiffness(m00, m01, m02, m11, m12, m22, n02, n12, n22, tbt, tbt2, obt, obt2, nn, DVD)

def MicroPlaneC1(ElemDim, ElemPlSt, PrinStrains, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM, type, E_V, E_DM, KV0, kV1, kV2, kap0V, alphaV, betaV, RType, ElemCrBws, gam2, kapUlt, ElemScaleType, eps_ct, e0, ed, gd, nu, Emod, Dps, etaV, Dt, sVsTol, DataOut, nState, nInt, PlStressI, PlStressL, iS):
    return _ConFemMatC.MicroPlaneC1(ElemDim, ElemPlSt, PrinStrains, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM, type, E_V, E_DM, KV0, kV1, kV2, kap0V, alphaV, betaV, RType, ElemCrBws, gam2, kapUlt, ElemScaleType, eps_ct, e0, ed, gd, nu, Emod, Dps, etaV, Dt, sVsTol, DataOut, nState, nInt, PlStressI, PlStressL, iS)

def MicroPlaneC2(ElemDim, ElemPlSt, PrinStrains, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM, type, E_V, E_DM, KV0, kV1, kV2, kap0V, alphaV, betaV, RType, ElemCrBws, gam2, kapUlt, ElemScaleType, eps_ct, e0, ed, gd, nu, Emod, Dps, etaV, Dt, sVsTol, DataOut, nState, nInt, PlStressI, PlStressL, iS):
    return _ConFemMatC.MicroPlaneC2(ElemDim, ElemPlSt, PrinStrains, ElemLch, ElemStateVar, ElemStateVarN, Eps, sig, MatM, type, E_V, E_DM, KV0, kV1, kV2, kap0V, alphaV, betaV, RType, ElemCrBws, gam2, kapUlt, ElemScaleType, eps_ct, e0, ed, gd, nu, Emod, Dps, etaV, Dt, sVsTol, DataOut, nState, nInt, PlStressI, PlStressL, iS)

