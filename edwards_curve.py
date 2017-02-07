#from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
#from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field
import sage.schemes.projective.projective_space as projective_space
from sage.all import QQ
import sage.schemes.plane_curves.projective_curve as plane_curve
from sage.schemes.elliptic_curves.edwards_curve_point import EdwardsCurvePoint
import sage.rings.all as rings
from sage.structure.sequence import Sequence
from sage.rings.arith import is_square
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.functions.other import sqrt
from sage.sets.set import Set
from sage.categories.morphism import Morphism
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories import homset
#Q = RationalField()
def is_EdwardsCurve(x):
    """
    Utility function to test if ``x`` is an instance of an Elliptic Curve class.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        sage: E = EllipticCurve([1,2,3/4,7,19])
        sage: is_EllipticCurve(E)
        True
        sage: is_EllipticCurve(0)
        False
    """
    return isinstance(x, EdwardsCurve)

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
        #relax the follwing condition for the sake of isogeny codomains
        #if is_square(self._d):
        #    raise ValueError("An Edwards Curve cannot have d which is a square in the base field")
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
        
    def zero(self):
        return self(0)


        
class EdwardsCurveIsogeny(Morphism):
    def __init__(self, E, kernel):
        if not is_EdwardsCurve(E):
            raise TypeError("E parameter must be an EdwardsCurve.")

        if not isinstance(kernel, type([1,1])) and kernel in E :
            # a single point was given, we put it in a list
            # the first condition assures that [1,1] is treated as x+1
            kernel = [E(kernel)]
        
        if isinstance(kernel, list):
            for P in kernel:
                if P not in E:
                    raise ValueError("The generators of the kernel of the isogeny must first lie on the EdwardsCurve")    
        
        self.__E1 = E #domain curve
        self.__E2 = None #codomain curve; not yet found
        
        self.__degree = None
        self.__base_field = E.base_ring()
        
        self.__kernel_gens = kernel
        self.__kernel_list = None
        self.__degree = None
        
        self.__poly_ring = PolynomialRing(self.__base_field, ['x','y'])

        self.__x_var = self.__poly_ring('x')
        self.__y_var = self.__poly_ring('y')
        
        # to determine the codomain, and the x and y rational maps
        self.__initialize_kernel_list()
        self.__compute_B()
        self.__compute_E2()
        self.__initialize_rational_maps()
        self._domain = self.__E1
        self._codomain = self.__E2

        # sets up the parent
        parent = homset.Hom(self.__E1(0).parent(), self.__E2(0).parent())
        Morphism.__init__(self, parent)
               
    def __compute_E2(self):
        """An internal function that sets the codomain of the isogeny
        """
        self.__d_cap = self.__B**8 * ( self.__E1.get_d()**self.__degree )
        print self.__d_cap
        self.__E2 = EdwardsCurve(self.__base_field, self.__d_cap)
    
    def codomain(self):
        return self.__E2
    
    def __initialize_rational_maps(self):
        pass
    
    def __initialize_kernel_list(self):
        """An internal funtion to initialize the kernel with the given generators; it also sets the degree of the isogeny
        """
        for P in self.__kernel_gens:
            if not P.has_finite_order():
                raise ValueError("The points in the kernel must be of finite order")
        
        kernel_list = Set([self.__E1(0)])
        for P in self.__kernel_gens:
            P = self.__E1(P)
            points_to_add = []
            for j in range(P.order()):
                for Q in kernel_list:
                    points_to_add.append(j*P+Q)
            kernel_list += Set(points_to_add)

        self.__kernel_list = kernel_list.list()
        self.__degree = len(self.__kernel_list)
        
    def __compute_B(self):
        """An internal function that computes the value B as it appears in the analogue of Velu's formula
        """
        self.__reduced_kernel_list = self.__kernel_list
        for p in self.__reduced_kernel_list:
            test_value = (self.__E1.base_ring().characteristic() - 1) / 2
            if p[0] > test_value:
                self.__reduced_kernel_list.remove(p)
        
        self.__reduced_kernel_list.remove(self.__E1(0))
        beta = 1
        for p in self.__reduced_kernel_list:
            beta *= p[1]
        self.__B = beta


def edwards_curve_from_weierstrass(E):
    if not isinstance(E, EllipticCurve):
        raise TypeError("This module can only be used to obtain the isomorphic Edwards model for a given Weierstrass model")
    for p in E.points():
        if p.order() == 4:
            break
    if not ( 2*p == (0,0,1) ):
        pass
    
    
    
        
        
