import py_ecc.bn128 as b
from utils import *
from dataclasses import dataclass
from curve import *
from transcript import Transcript
from poly import Polynomial, Basis


@dataclass
class VerificationKey:
    """Verification key"""

    # we set this to some power of 2 (so that we can FFT over it), that is at least the number of constraints we have (so we can Lagrange interpolate them)
    group_order: int
    # [q_M(x)]₁ (commitment to multiplication selector polynomial)
    Qm: G1Point
    # [q_L(x)]₁ (commitment to left selector polynomial)
    Ql: G1Point
    # [q_R(x)]₁ (commitment to right selector polynomial)
    Qr: G1Point
    # [q_O(x)]₁ (commitment to output selector polynomial)
    Qo: G1Point
    # [q_C(x)]₁ (commitment to constants selector polynomial)
    Qc: G1Point
    # [S_σ1(x)]₁ (commitment to the first permutation polynomial S_σ1(X))
    S1: G1Point
    # [S_σ2(x)]₁ (commitment to the second permutation polynomial S_σ2(X))
    S2: G1Point
    # [S_σ3(x)]₁ (commitment to the third permutation polynomial S_σ3(X))
    S3: G1Point
    # [x]₂ = xH, where H is a generator of G_2
    X_2: G2Point
    # nth root of unity (i.e. ω^1), where n is the program's group order.
    w: Scalar

    """ 
    NOTE: 
    These are commitment, and CommonPreprocessedInput contains the polynomials that are committed to.
    """

    def __init__(self, pk, setup):
        # Assign the attributes of pk and setup to the VerificationKey object
        self.group_order = pk.group_order
        self.Qm = setup.commit(pk.QM)
        self.Ql = setup.commit(pk.QL)
        self.Qr = setup.commit(pk.QR)
        self.Qo = setup.commit(pk.QO)
        self.Qc = setup.commit(pk.QC)
        self.S1 = setup.commit(pk.S1)
        self.S2 = setup.commit(pk.S2)
        self.S3 = setup.commit(pk.S3)
        self.X_2 = setup.X2
        self.w = Scalar.root_of_unity(pk.group_order)
        """ 
        NOTE:
        # What are these attributes?
        1. group_order: This is the number of points in the group G1 or G2. 
        It is a power of 2 that is at least the number of constraints in the circuit. 
        It is used to perform Fourier transforms over the Lagrange basis

        2. Qm, Ql, Qr, Qo, Qc: These are commitments to the selector polynomials that define the circuit constraints. 
        They are elements of G1 that are multiplied by the corresponding wire polynomials to form the arithmetic gate equation
        
        3. S1, S2, S3: These are commitments to the permutation polynomials that enforce copy constraints between different wires. 
        They are elements of G1 that are used to construct the permutation argument
        
        4. X_2: This is a commitment to W in the group G2. 
        It is used to verify that the prover knows a polynomial that satisfies the circuit constraints 

        5. w: This is W as a scalar. It is the nth root of unity that is used to perform Fourier transforms over the Lagrange basis
        
        # Get the nth root of unity, where n is the program's group order.
        ## Method 1:
        To use the root_of_unity method, you need to know the group order of your circuit. 
        The group order is the number of points in the group G1 or G2. 
        It is usually a power of 2 that is at least the number of constraints in your circuit. 
        For example, if your circuit has 1024 constraints, you can use a group order of 2048. 
        Then you can call the root_of_unity method with the group order as an argument, 
        and it will return W as a scalar 

        you can convert it to a scalar by using the discrete logarithm function. 
        The root_of_unity method is a way to compute W from the group order, 
        but it is not the only way. You can also get W from the setup file, 
        if you have access to it. 
        The advantage of using the setup file is that you don’t need to know the group order in advance.
        
        ## Method 2:
        To get w, you need to extract the second element of the powers_of_x list, 
        which is [x]₁. Then you need to convert it to a scalar by using the discrete logarithm function. 
        For example, if you are using Python, you can use the sympy library:
            from sympy.ntheory import discrete_log
            w = discrete_log(b.field_modulus, powers_of_x[1], b.G1[0])
        """

    # More optimized version that tries hard to minimize pairings and
    # elliptic curve multiplications, but at the cost of being harder
    # to understand and mixing together a lot of the computations to
    # efficiently batch them
    def verify_proof(self, group_order: int, pf, public=[]) -> bool:
        # 4. Compute challenges
        self.compute_challenges(group_order, pf)
        # 5. Compute zero polynomial evaluation Z_H(ζ) = ζ^n - 1

        # 6. Compute Lagrange polynomial evaluation L_0(ζ)

        # 7. Compute public input polynomial evaluation PI(ζ).

        # Compute the constant term of R. This is not literally the degree-0
        # term of the R polynomial; rather, it's the portion of R that can
        # be computed directly, without resorting to elliptic cutve commitments

        # Compute D = (R - r0) + u * Z, and E and F

        # Run one pairing check to verify the last two checks.
        # What's going on here is a clever re-arrangement of terms to check
        # the same equations that are being checked in the basic version,
        # but in a way that minimizes the number of EC muls and even
        # compressed the two pairings into one. The 2 pairings -> 1 pairing
        # trick is basically to replace checking
        #
        # Y1 = A * (X - a) and Y2 = B * (X - b)
        #
        # with
        #
        # Y1 + A * a = A * X
        # Y2 + B * b = B * X
        #
        # so at this point we can take a random linear combination of the two
        # checks, and verify it with only one pairing.

        return False

    # Basic, easier-to-understand version of what's going on
    def verify_proof_unoptimized(self, group_order: int, pf, public=[]) -> bool:
        # 4. Compute challenges

        # 5. Compute zero polynomial evaluation Z_H(ζ) = ζ^n - 1

        # 6. Compute Lagrange polynomial evaluation L_0(ζ)

        # 7. Compute public input polynomial evaluation PI(ζ).

        # Recover the commitment to the linearization polynomial R,
        # exactly the same as what was created by the prover

        # Verify that R(z) = 0 and the prover-provided evaluations
        # A(z), B(z), C(z), S1(z), S2(z) are all correct

        # Verify that the provided value of Z(zeta*w) is correct

        return False

    # Compute challenges (should be same as those computed by prover)
    def compute_challenges(
        self, proof
    ) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        beta, gamma = transcript.round_1(proof.msg_1)
        alpha, _fft_cofactor = transcript.round_2(proof.msg_2)
        zeta = transcript.round_3(proof.msg_3)
        v = transcript.round_4(proof.msg_4)
        u = transcript.round_5(proof.msg_5)

        return beta, gamma, alpha, zeta, v, u
