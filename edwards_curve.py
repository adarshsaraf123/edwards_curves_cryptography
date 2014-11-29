#from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
#from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field
import sage.schemes.projective.projective_space as projective_space
from sage.all import QQ
import sage.schemes.plane_curves.projective_curve as plane_curve
from edwards_curve_point import EdwardsCurvePoint
import sage.rings.all as rings
from sage.structure.sequence import Sequence
from sage.rings.arith import is_square
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.functions.other import sqrt

#Q = RationalField()

class EdwardsCurve(plane_curve.ProjectiveCurve_generic):
    Element = EdwardsCurvePoint
    def __init__(self, d, extra = None):
        if extra != None:
            K, self._d = d, extra
        else: 
            if isinstance(d, (rings.Rational, rings.Integer, int, long)):
                K = QQ
            else:
                K = d.parent()
            self._d = d
        self._d = K(self._d)
        if self._d.is_zero() or self._d.is_one():
            raise ValueError("An Edwards Curve cannot have d as 1 or 0 of the base field")
        if is_square(self._d):
            raise ValueError("An Edwards Curve cannot have d which is a square in the base field")
        if K.characteristic() == 2:
            raise ValueError("the base field cannot have characteristic two")
        self.__base_ring = K
        PP = projective_space.ProjectiveSpace(2, K, names='xyz');
        x, y, z = PP.coordinate_ring().gens()
        self.f = x**2 * z**2 + y**2 * z**2 - z**4 - self._d * x**2 * y**2
        plane_curve.ProjectiveCurve_generic.__init__(self, PP, self.f)
        # super(MyEdwardsCurve,self).__init__([5,6])
        #print self.__A
    def _repr_(self):
        s = "Edwards curve defined by x^2 + y^2 = 1 "
        if self.get_d() != 1:
            s += "+ %s" %self.get_d()
        s += "x^2y^2"
        s += " over %s" %self.base_ring() 
        s = s.replace("+ -","- ")
        return s
        #return "Edwards curve defined by the equation %s over %s" %(s, self.base_ring())
    
    def get_d(self):
        return self._d
    
    def base_ring(self):
        return self.__base_ring
    
    def __call__(self, *args, **kwgs):
        return EdwardsCurvePoint(self, args[0])
    
    def points(self, limit = False):
        try:
            return self.__points
        except AttributeError: pass
        
        if not self.base_ring().is_finite():
            raise TypeError("points() is defined only for Edwards curve over finite fields")
        v = []
        for x in self.base_ring():
            test = (x**2 - 1)/(self._d*(x**2) - 1)
            if is_square(test):
                y = sqrt(test)
                v.append(self([x,y,1]))
                v.append(self([x,-y,1]))
        v.sort()
        self.__points = Sequence(v,immutable = True)
        return self.__points
    
    def random_point(self):
        while(1):
            x = self.base_ring().random_element()
            test = (x**2 - 1)/(self._d*(x**2) - 1)
            if is_square(test):
                y = sqrt(test)
                return self([x,y,1])
    
    def n_random_points(self,n):
        v = []
        while(len(v) < n):
            x = self.base_ring().random_element()
            test = (x**2 - 1)/(self._d*(x**2) - 1)
            if is_square(test):
                y = sqrt(test)
                v.append(self([x,y,1]))
        v.sort()
        return v        
        
    def torsion_points(self, n):
        torsion_points_list = []
        for p in self.points():
            if n % p.order() == 0:
                torsion_points_list.append(p)
        return torsion_points_list
    
    def weierstrass_curve(self):
        try:
            return self.__weierstrass_curve
        except AttributeError:
            pass
        d = self.get_d()
        self.__weierstrass_curve = EllipticCurve(self.base_ring(), [0, -(d+1), 0, -4*d, 4*d*(d+1)])
        return self.__weierstrass_curve    
    
    def get_s(self):
        s = []
        d = self.get_d()
        if is_square(-d):
            c = ( 2*(d-1) + 4*sqrt(-d) ) / ( 2*(d+1) )
            s_2 = 2/c
            if is_square(s_2):
                s.append(sqrt(s_2))
            
            c = ( 2*(d-1) - 4*sqrt(-d) ) / ( 2*(d+1) )
            s_2 = 2/c
            if is_square(s_2):
                s.append(sqrt(s_2))
        return s    
        
    def zero(self):
        return self(0)
    
        
    
        
        
