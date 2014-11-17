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
#Q = RationalField()

class EdwardsCurve(plane_curve.ProjectiveCurve_generic):
    Element = EdwardsCurvePoint
    def __init__(self, d, extra = None):
        if extra != None:
            K, self.__d = d, extra
        else: 
            if isinstance(d, (rings.Rational, rings.Integer, int, long)):
                K = QQ
            else:
                K = d.parent()
            self.__d = d
        self.__d = K(self.__d)
        if self.__d.is_zero() or self.__d.is_one():
            raise ValueError("An Edwards Curve cannot have d as 1 or 0 of the base field")
        if is_square(self.__d):
            raise ValueError("An Edwards Curve cannot have d which is a square in the base field")
        if K.characteristic() == 2:
            raise ValueError("the base field cannot have characteristic two")
        self.__base_ring = K
        PP = projective_space.ProjectiveSpace(2, K, names='xyz');
        x, y, z = PP.coordinate_ring().gens()
        self.f = x**2 * z**2 + y**2 * z**2 - z**4 - self.__d * x**2 * y**2
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
        return self.__d
    def base_ring(self):
        return self.__base_ring
    def __call__(self, *args, **kwgs):
        return EdwardsCurvePoint(self, args[0])
    def points(self):
        try:
            return self.__points
        except AttributeError: pass
        
        if not self.base_ring().is_finite():
            raise TypeError("points() is defined only for Edwards curve over finite fields")
        v = []
        for x in self.base_ring():
            for y in self.base_ring():
                try:
                    v.append(self([x,y,1]))
                except (TypeError,ValueError):
                    pass
        v.sort()
        self.__points = Sequence(v,immutable = True)
        return self.__points
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
        d = self.__d
        self.__weierstrass_curve = EllipticCurve(self.base_ring(), [0, -(d+1), 0, -4*d, 4*d*(d+1)])
        return self.__weierstrass_curve    
        
    def zero(self):
        return self(0)
    
        
    
        
        
