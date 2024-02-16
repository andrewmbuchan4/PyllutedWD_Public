#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum

class Element(Enum):
    Placeholder = 0  # Can be used as a 'Hx' indicator if necessary!
    H  = 1
    He = 2
    Li = 3
    Be = 4
    B  = 5
    C  = 6
    N  = 7
    O  = 8
    F  = 9
    Ne = 10
    Na = 11
    Mg = 12
    Al = 13
    Si = 14
    P  = 15
    S  = 16
    Cl = 17
    Ar = 18
    K  = 19
    Ca = 20
    Sc = 21
    Ti = 22
    V  = 23
    Cr = 24
    Mn = 25
    Fe = 26
    Co = 27
    Ni = 28
    Cu = 29
    Zn = 30
    Ga = 31
    Ge = 32
    As = 33
    Se = 34
    Br = 35
    Kr = 36
    Rb = 37
    Sr = 38
    Y  = 39
    Zr = 40
    Nb = 41
    Mo = 42
    Tc = 43
    Ru = 44
    Rh = 45
    Pd = 46
    Ag = 47
    Cd = 48
    In = 49
    Sn = 50
    Sb = 51
    Te = 52
    I  = 53
    Xe = 54
    Cs = 55
    Ba = 56
    La = 57
    Ce = 58
    Pr = 59
    Nd = 60
    Pm = 61
    Sm = 62
    Eu = 63
    Gd = 64
    Tb = 65
    Dy = 66
    Ho = 67
    Er = 68
    Tm = 69
    Yb = 70
    Lu = 71
    Hf = 72
    Ta = 73
    W  = 74
    Re = 75
    Os = 76
    Ir = 77
    Pt = 78
    Au = 79
    Hg = 80
    Tl = 81
    Pb = 82
    Bi = 83
    Po = 84
    At = 85
    Rn = 86
    Fr = 87
    Ra = 88
    Ac = 89
    Th = 90
    Pa = 91
    U  = 92
    Np = 93
    Pu = 94
    Am = 95
    Cm = 96
    Bk = 97
    Cf = 98
    Es = 99
    Fm = 100
    Md = 101
    No = 102
    Lr = 103
    Rf = 104
    Db = 105
    Sg = 106
    Bh = 107
    Hs = 108
    Mt = 109
    Ds = 110
    Rg = 111
    Cn = 112
    Nh = 113
    Fl = 114
    Mc = 115
    Lv = 116
    Ts = 117
    Og = 118

    def __str__(self):
        return self.name

    def __hash__(self):
        # For performance purposes:
        # Normally python will try to hash an Enum by hashing its name (which guaranteess uniqueness)
        # In this case I know that the values are also unique, and should satisfy properties of a hash
        # So I do this, which is much quicker:
        return self.value

    def __len__(self):
        return 1

    def __lt__(self, other):
        return (self.value < other.value)

    def __le__(self, other):
        return(self.value <= other.value)

    def __gt__(self, other):
        return(self.value > other.value)

    def __ge__(self, other):
        return(self.value >= other.value)

    def __eq__(self, other):
        if not isinstance(other, Element):
            return False
        return (self.value == other.value)

    def __ne__(self, other):
        return not(self.__eq__(other))


all_elements = [
    Element.H,
    Element.He,
    Element.Li,
    Element.Be,
    Element.B,
    Element.C,
    Element.N,
    Element.O,
    Element.F,
    Element.Ne,
    Element.Na,
    Element.Mg,
    Element.Al,
    Element.Si,
    Element.P,
    Element.S,
    Element.Cl,
    Element.Ar,
    Element.K,
    Element.Ca,
    Element.Sc,
    Element.Ti,
    Element.V,
    Element.Cr,
    Element.Mn,
    Element.Fe,
    Element.Co,
    Element.Ni,
    Element.Cu,
    Element.Zn,
    Element.Ga,
    Element.Ge,
    Element.As,
    Element.Se,
    Element.Br,
    Element.Kr,
    Element.Rb,
    Element.Sr,
    Element.Y,
    Element.Zr,
    Element.Nb,
    Element.Mo,
    Element.Tc,
    Element.Ru,
    Element.Rh,
    Element.Pd,
    Element.Ag,
    Element.Cd,
    Element.In,
    Element.Sn,
    Element.Sb,
    Element.Te,
    Element.I,
    Element.Xe,
    Element.Cs,
    Element.Ba,
    Element.La,
    Element.Ce,
    Element.Pr,
    Element.Nd,
    Element.Pm,
    Element.Sm,
    Element.Eu,
    Element.Gd,
    Element.Tb,
    Element.Dy,
    Element.Ho,
    Element.Er,
    Element.Tm,
    Element.Yb,
    Element.Lu,
    Element.Hf,
    Element.Ta,
    Element.W,
    Element.Re,
    Element.Os,
    Element.Ir,
    Element.Pt,
    Element.Au,
    Element.Hg,
    Element.Tl,
    Element.Pb,
    Element.Bi,
    Element.Po,
    Element.At,
    Element.Rn,
    Element.Fr,
    Element.Ra,
    Element.Ac,
    Element.Th,
    Element.Pa,
    Element.U,
    Element.Np,
    Element.Pu,
    Element.Am,
    Element.Cm,
    Element.Bk,
    Element.Cf,
    Element.Es,
    Element.Fm,
    Element.Md,
    Element.No,
    Element.Lr,
    Element.Rf,
    Element.Db,
    Element.Sg,
    Element.Bh,
    Element.Hs,
    Element.Mt,
    Element.Ds,
    Element.Rg,
    Element.Cn,
    Element.Nh,
    Element.Fl,
    Element.Mc,
    Element.Lv,
    Element.Ts,
    Element.Og
]

element_masses = {
    Element.H: 1.008,
    Element.He: 4.002602,
    Element.Li: 6.94,
    Element.Be: 9.0121831,
    Element.B: 10.81,
    Element.C: 12.011,
    Element.N: 14.007,
    Element.O: 15.999,
    Element.F: 18.998403163,
    Element.Ne: 20.1797,
    Element.Na: 22.98976928,
    Element.Mg: 24.305,
    Element.Al: 26.9815384,
    Element.Si: 28.085,
    Element.P: 30.973761998,
    Element.S: 32.06,
    Element.Cl: 35.45,
    Element.Ar: 39.95,
    Element.K: 39.0983,
    Element.Ca: 40.078,
    Element.Sc: 44.955908,
    Element.Ti: 47.867,
    Element.V: 50.9415,
    Element.Cr: 51.9961,
    Element.Mn: 54.938043,
    Element.Fe: 55.845,
    Element.Co: 58.933194,
    Element.Ni: 58.6934,
    Element.Cu: 63.546,
    Element.Zn: 65.38,
    Element.Ga: 69.723,
    Element.Ge: 72.630,
    Element.As: 74.921595,
    Element.Se: 78.971,
    Element.Br: 79.904,
    Element.Kr: 83.798,
    Element.Rb: 85.4678,
    Element.Sr: 87.62,
    Element.Y: 88.90584,
    Element.Zr: 91.224,
    Element.Nb: 92.90637,
    Element.Mo: 95.95,
    Element.Tc: 98,
    Element.Ru: 101.07,
    Element.Rh: 102.90549,
    Element.Pd: 106.42,
    Element.Ag: 107.8682,
    Element.Cd: 112.414,
    Element.In: 114.818,
    Element.Sn: 118.710,
    Element.Sb: 121.760,
    Element.Te: 127.60,
    Element.I: 126.90447,
    Element.Xe: 131.293,
    Element.Cs: 132.90545196,
    Element.Ba: 137.327,
    Element.La: 138.90547,
    Element.Ce: 140.116,
    Element.Pr: 140.90766,
    Element.Nd: 144.242,
    Element.Pm: 145,
    Element.Sm: 150.36,
    Element.Eu: 151.964,
    Element.Gd: 157.25,
    Element.Tb: 158.925354,
    Element.Dy: 162.500,
    Element.Ho: 164.930328,
    Element.Er: 167.259,
    Element.Tm: 168.934218,
    Element.Yb: 173.045,
    Element.Lu: 174.9668,
    Element.Hf: 178.49,
    Element.Ta: 180.94788,
    Element.W: 183.84,
    Element.Re: 186.207,
    Element.Os: 190.23,
    Element.Ir: 192.217,
    Element.Pt: 195.084,
    Element.Au: 196.966570,
    Element.Hg: 200.592,
    Element.Tl: 204.38,
    Element.Pb: 207.2,
    Element.Bi: 208.98040,
    Element.Po: 209,
    Element.At: 210,
    Element.Rn: 222,
    Element.Fr: 223,
    Element.Ra: 226,
    Element.Ac: 227,
    Element.Th: 232.0377,
    Element.Pa: 231.03588,
    Element.U: 238.02891,
    Element.Np: 237,
    Element.Pu: 244,
    Element.Am: 243,
    Element.Cm: 247,
    Element.Bk: 247,
    Element.Cf: 251,
    Element.Es: 252,
    Element.Fm: 257,
    Element.Md: 258,
    Element.No: 259,
    Element.Lr: 266,
    Element.Rf: 267,
    Element.Db: 268,
    Element.Sg: 269,
    Element.Bh: 270,
    Element.Hs: 270,
    Element.Mt: 278,
    Element.Ds: 281,
    Element.Rg: 282,
    Element.Cn: 285,
    Element.Nh: 286,
    Element.Fl: 289,
    Element.Mc: 290,
    Element.Lv: 293,
    Element.Ts: 294,
    Element.Og: 294
}

def get_element_mass(element):
    return element_masses[element]

usual_elements = [
    Element.Al,
    Element.Ti,
    Element.Ca,
    Element.Ni,
    Element.Fe,
    Element.Cr,
    Element.Mg,
    Element.Si,
    Element.Na,
    Element.O,
    Element.C,
    Element.N
]

atmospheric_types = [
    Element.H,
    Element.He
]
