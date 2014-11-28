from sage.crypto.cryptosystem import PublicKeyCryptosystem
import safe_curves
from exceptions import TypeError, AttributeError
from sage.rings.integer import Integer

class EllipticCurveCryptosystem(PublicKeyCryptosystem):
    def __init__(self, curve, base_point):
        self._curve = curve
        self._base_point = self._curve(base_point)
    
    def curve(self):
        return self._curve
    
    def base_point(self):
        return self._base_point

class ECDH(EllipticCurveCryptosystem):
    def __init__(self, curve = safe_curves.CurveForTesting(), base_point):
        EllipticCurveCryptosystem.__init__(self, curve, base_point)
   
    def decrypt(self, ciphertext_pair, B_private_key):
        first, second = ciphertext_pair
        return (second - (first * B_private_key) )
    def encrypt(self, plaintext, A_private_key, B_public_key):
        r"""
        For convenience of understanding we state the encryption process to be:
            Alice wants to encrypt the 'plaintext' using
                * Alice's private key: 'A_private_key'; it is an integer
                * Bob's public key: 'B_public_key'; it is a point on the elliptic curve on which the system is based
        Returns a tuple (A_private_key*curve.base_point(), self.encode(plaintext) + A_private_key*B_public_key)
        """
        #P_m = self.encode(plaintext) # plain message point
        P_m = plaintext
        try:
            P = self.curve()((B_public_key[0],B_public_key[1]))
        except TypeError:
            raise TypeError("Bob's public key not on the curve of the cryptosystem")
        try:
            k = Integer(A_private_key)
        except TypeError:
            raise TypeError("Alice's private key not a valid integer")
        try:
            G = self.curve().base_point()
        except AttributeError:
            raise AttributeError("No base point defined for the curve")
        
        return (k*G, P_m + k*P) 
    
    def encode(self, message):
        return self.curve()(message)