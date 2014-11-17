from sage.schemes.projective.projective_point import SchemeMorphism_point_abelian_variety_field
from sage.rings.integer import Integer
class EdwardsCurvePoint(SchemeMorphism_point_abelian_variety_field):
    def __init__(self, curve, v, check = True):
        point_homset = curve.point_homset()
        if v == 0:
            v = (0, 1, 1) 
        SchemeMorphism_point_abelian_variety_field.__init__(self, point_homset, v, check=True)
    def _add_(self, right):
        if self.is_zero():
            return right
        if right.is_zero():
            return self
        E = self.curve()
        d = E.get_d()
        x1, y1 = self[0], self[1]
        if self == right:
            x3 = (2*x1*y1)/(x1**2 + y1**2)
            y3 = (y1**2 - x1**2)/(2 - x1**2 - y1**2)
        else:
            x2, y2 = right[0], right[1]
            temp = d*x1*y1*x2*y2
            x3 = (x1*y2 + y1*x2)/(1+temp)
            y3 = (y1*y2 - x1*x2)/(1-temp) 
        return E([x3, y3, 1], check=False)
    def _lmul_(self, c):
        inverse = False
        if c < 0:
            c *= -1
            inverse = True
        c = c % self.order()        
        bin_c = Integer(c).binary() #taking care of negative c
        Q = self.curve()(0)
        for i in range(len(bin_c)):
            Q = Q+Q
            if bin_c[i] == '1':
                Q += self
        if inverse: 
            Q = self.curve()( (-1*Q[0],Q[1],Q[2]), check = False )
        return Q
    def __nonzero__(self):
        return not(( not bool(self[0]) ) and (self[1] == self[2]) )
    def curve(self):
        return self.scheme()
    def additive_order(self):
        try:
            return self._order
        except AttributeError: pass
             
        Q = self.curve().zero()
        i = 1
        while(1):
            Q += self
            if Q.is_zero():
                self._order = i
                break
            i += 1
        return self._order
    def weierstrass_point(self):
        try:
            return self.__weierstrass_point
        except AttributeError:
            pass
        E = self.curve()
        u = self[0]
        v = self[1]
        d = E.get_d()
        w = (d*u*u - 1)*v
        x = ( -2 * (w-1) ) / (u*u)
        y = ( 4*(w-1) + 2*(d+1)*u*u ) / (u**3)
        self.__weierstrass_point = E.weierstrass_curve()([x,y])
        return self.__weierstrass_point
    #===========================================================================
    # def __init__(self,curve,v,check = False):
    #     point_homset = curve.point_homset()
    #     if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
    #         v = list(v)
    #     elif v == 0:
    #         # some of the code assumes that E(0) has integral entries
    #         # irregardless of the base ring...
    #         #R = self.base_ring()
    #         #v = (R.zero(),R.one(),R.zero())
    #         v = (0, 1, 0)
    #     SchemeMorphism_point_abelian_variety_field.__init__(self, point_homset, v, check=False)
    #===========================================================================