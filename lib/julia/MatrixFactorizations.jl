module MatrixFactorizations

using QuasiHomogeneous

using PolynomialRings

export StrangeDuality, S1, T2, Dk_1, E6_1, E7_3, one_certain_unit_mf, E14_Q10

StrangeDuality() = begin
    A = @ring! ℤ[x,y,z]

    d1 = [z   y^2 x^3 0 ; y -x*z 0 x^3 ; x 0 -x*z -y^2 ; 0 x -y z  ]
    d0 = [x*z y^2 x^3 0 ; y -z   0 x^3 ; x 0 -z   -y^2 ; 0 x -y x*z]
    zz = zero(d1)

    [ zz d1; d0 zz ]
end

S1() = begin
    A = @ring! ℤ[x,y]

    q1=[x x*y; -y -x^3]
    q0=[x^3 x*y; -y -x]
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

T2() = begin
    A = @ring! ℤ[x,y]

    q1 = [x^2 y ; -y -x^2]
    q0 = [x^3 x*y; -x*y -x^3]
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

# The following are the factorizations from [22]

Dk_1(k::Integer, l::Integer) = begin
    A = @ring! ℤ[x,y,z]

    if k == 1
        q1 = [z x^2 + y^(j-2); y -z]
        q0 = q1
    elseif k % 2 == 0 && 2 <= k && k <= l - 2
        q1 = [-z 0 x*y y^(k>>1);
              0 -z y^(k-1-(k>>1)) -x;
              x y^(k>>1) z 0;
              y^(k-1-(k>>1)) -x*y 0 z
            ]
        q0 = q1
    elseif k % 2 == 1 && 3 <= k && k <= l - 2
        q1 = [-z y^(k>>1) x*y 0;
              y^(l - (k+1)>>1) z 0 -x;
              x 0 z y^((k-3)>>1);
              0 -x*y y^(l-(k>>1)) -z
            ]
        q0 = q1
    elseif k <= l
        throw(ArgumentError("not implemented yet; need to test complex numbers"))
    else
        throw(ArgumentError("k needs to be at most l"))
    end
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

E6_1() = begin
    A = @ring! ℤ[x,y,z]

    q1 = [-z 0 x^2 y^3; 0 -z y -x; x y^3 z 0; y -x^2 0 z]
    q0 = q1
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

E7_3() = begin
    A = @ring! ℚ[x,y,z]

    q1 = [-z 0 x*y -y^2 0 0 x^2 0;
          0 -z 0 y^2 0 0 0 x;
          y^2 y^2 z 0 0 -x 0 0 ;
          0 x*y 0 z -x^2 0 0 0 ;
          0 0 0 -x -z 0 0 y ;
          0 0 -x^2 0 0 -z x*y^2 y^2 ;
          x 0 0 0 -y^2 y z 0 ;       # NOTE: typo in paper [22]: it lists -y^2 y^2 z as the nonzero coefficients
          0 x^2 0 0 x*y^2 0 0 z ]
    q0 = q1
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

one_certain_unit_mf() = begin
    A = @ring! ℚ[x,y,z,u,v,w]

    q1 = [ x^3 + x^2*u + x*u^2 + u^3 + z^2        0     y^2 + y*v + v^2      z*u + u*w;
           0      x^3 + x^2*u + x*u^2 + u^3 + z^2 -z + w                           y - v;
        -y + v          z*u + u*w                           x - u                        0;
      -z + w         -y^2 - y*v - v^2                          0                     x - u]

    q0 = [  x - u         0     -y^2 - y*v - v^2     -z*u - u*w;
          0          x - u        z - w                -y + v ;
       y - v      -z*u - u*w  x^3 + x^2*u + x*u^2 + u^3 + z^2   0 ;
       z - w        y^2 + y*v + v^2   0   x^3 + x^2*u + x*u^2 + u^3 + z^2 ]

    zz = zero(q1)

    [ zz q1; q0 zz ]
end

E14_Q10() = begin
    A = @ring! ℚ[x,y,z,u,v,w]
    return (u^4*w + v^3 + w^2) - (x^4 + y^3 + x*z^2)
end

J3_0_Z13() = begin
    A = @ring! ℚ[x,y,z,u,v,w]
    return x^6*y + y^3 + z^2 - u^6 - u*v^3 - w^2
end


unit_matrix_factorization(f, left_vars, right_vars) = begin



end

using PolynomialRings
using QuasiHomogeneous: generic_matrix_factorizations
function E14_Q10_possibilities(highest_free_generator_grading_source=7, highest_free_generator_grading_target=7)
    R, vars = polynomial_ring(:x,:y,:z,:u,:v,:w)

    generic_matrix_factorizations(4,highest_free_generator_grading_source, highest_free_generator_grading_target,24,(6,8,9,3,8,12),R,:c)
end

function E14_Q10_grading()
    d = [
         12  15  16  19;
          9  12  13  16;
          8  11  12  15;
          5   8   9  12
    ]
    z = fill(-1, 4,4)
    [z d;d z]
end

function E14_Q10_matrix_and_equations()
    R = @ring! ℚ[a,b,c]
    S = @ring! R[x,y,z,u,v,w]

    κ1 = a^3//2 + a^2*b + a*b^2 - a^2*c//2 - a*b*c
    κ2 = 1 + 3a^4//4 + 3a^3*b + 4a^2*b^2 + 2a*b^3  - a^3*c - 3*a^2*b*c - 2a*b^2*c

    x15 = κ1*u^3+a*u*x+z
    x16 = v^2+v*y+y^2
    x17 = κ2*u^4//2+w-a*(-a-2b)*u^2*x//2+x^2+b*u*z
    x25 = y-v
    x26 = (-b-b^2*κ1+(c-a)*κ2//2)*u^5+(-a-2b+c)*u*w+c*u*x^2+b*(-a-b+c)*u^2*z-x*z
    x35 = (-1+(-a-2b+c)*κ1+κ2//2)*u^4-w+a*(-a-2b+2c)*u^2*x//2+x^2+(-a-b+c)*u*z

    x73 = x62 = x48 = x15
    x74 = x16; x38 = x52 = -x16
    x28 = x17; x53 = x64 = -x17
    x83 = x25; x61 = x47 = -x25
    x51 = x84 = x37 = x26
    x46 = x35; x82 = x71 = -x35

    d1 = [x15 x16 x17   0;
          x25 x26   0 x28;
          x35   0 x37 x38;
            0 x46 x47 x48]
    d2 = [x51 x52 x53   0;
          x61 x62   0 x64;
          x71   0 x73 x74;
            0 x82 x83 x84]

     z = zero(d1)

     Q = [z d1; d2 z]

     f1 = -4+3a^4+8a^3*b+8a^2*b^2-4a^3*c-8*a^2*b*c
     f2 = 4+3a^4+8a^3*b+8a^2*b^2-4a^3*c-8a^2*b*c
     g  = a^2*(a^4-8a^2*b^2-16a*b^3-8b^4+8a^2*b*c+24a*b^2*c+16b^3*c-2a^2*c^2-8a*b*c^2-8b^2*c^2)

     Q, f1, f2, g
end

function J3_0_Z13_possibilities(highest_free_generator_grading_source, highest_free_generator_grading_target)
    R, vars = polynomial_ring(:x,:y,:z,:u,:v,:w)

    generic_matrix_factorizations(4,highest_free_generator_grading_source,highest_free_generator_grading_target,18,(2,6,9,3,5,9),R,:c)
end

# ----------------------------------------------
#
# Potentials from Strange Dualities II
#
# ----------------------------------------------

Q10()           = (@ring! ℤ[x,y,z]; x^4   + y^3   + x*z^2)
Q11()           = (@ring! ℤ[x,y,z]; x^3*y + y^3   + x*z^2)
Q12(::Val{:v1}) = (@ring! ℤ[x,y,z]; x^3*z + y^3   + x*z^2)
Q12(::Val{:v2}) = (@ring! ℤ[x,y,z]; x^5   + y^3   + x*z^2)
Q12(x::Symbol)  = Q12(Val{x}())
Q12()           = Q12(:v1)
S11()           = (@ring! ℤ[x,y,z]; x^4   + y^2*z + x*z^2)
S12()           = (@ring! ℤ[x,y,z]; x^3*y + y^2*z + x*z^2)
U12(::Val{:v1}) = (@ring! ℤ[x,y,z]; x^4   + y^3   + z^3)
U12(::Val{:v2}) = (@ring! ℤ[x,y,z]; x^4   + y^3   + z^2*y)
U12(::Val{:v3}) = (@ring! ℤ[x,y,z]; x^4   + y^2*z + z^2*y)
U12(x::Symbol)  = U12(Val{x}())
U12()           = U12(:v1)
Z11()           = (@ring! ℤ[x,y,z]; x^5   + x*y^3 + z^2)
Z12()           = (@ring! ℤ[x,y,z]; y*x^4 + x*y^3 + z^2)
Z13(::Val{:v1}) = (@ring! ℤ[x,y,z]; x^3*z + x*y^3 + z^2)
Z13(::Val{:v2}) = (@ring! ℤ[x,y,z]; x^6   + y^3*x + z^2)
Z13(x::Symbol)  = Z13(Val{x}())
Z13()           = Z13(:v1)
W12(::Val{:v1}) = (@ring! ℤ[x,y,z]; x^5   + y^2*z + z^2)
W12(::Val{:v2}) = (@ring! ℤ[x,y,z]; x^5   + y^4   + z^2)
W12(x::Symbol)  = W12(Val{x}())
W12()           = W12(:v1)
W13(::Val{:v1}) = (@ring! ℤ[x,y,z]; y*x^4 + y^2*z + z^2)
W13(::Val{:v2}) = (@ring! ℤ[x,y,z]; x^4*y + y^4   + z^2)
W13(x::Symbol)  = W13(Val{x}())
W13()           = W13(:v1)
E12()           = (@ring! ℤ[x,y,z]; x^7   + y^3   + z^2)
E13()           = (@ring! ℤ[x,y,z]; y^3   + y*x^5 + z^2)
E14(::Val{:v1}) = (@ring! ℤ[x,y,z]; x^4*z + y^3   + z^2)
E14(::Val{:v2}) = (@ring! ℤ[x,y,z]; x^8   + y^3   + z^2)
E14(x::Symbol)  = E14(Val{x}())


# ----------------------------------------------
#
# Matrix factorizations from Strange Dualities II
#
# ----------------------------------------------
function StandardDuality(;substitutions...)
    @ring! ℚ[d15,d16,d17,d25,d26,d35]

    [   0    0    0    0  d15  d16  d17    0;
        0    0    0    0  d25  d26    0  d17;
        0    0    0    0  d35    0  d26 -d16;
        0    0    0    0    0  d35 -d25  d15;
      d26 -d16 -d17    0    0    0    0    0;
     -d25  d15    0 -d17    0    0    0    0;
     -d35    0  d15  d16    0    0    0    0;
        0 -d35  d25  d26    0    0    0    0](;substitutions...)
end

function find_standard_duality(W, vars...)

    gradings = QuasiHomogeneous.find_quasihomogeneous_degrees(W, vars...)
    D = quasidegree(W, gradings)

    count = 0

    for deg_d17 = 1:(D-1)
        deg_d35 = D - deg_d17
        for deg_d15 = deg_d17:(D-1)
            deg_d26 = D - deg_d15
            for deg_d16 = deg_d15:(D-1)
                deg_d25 = D - deg_d16

                coeffs = Channel() do ch
                    for c in formal_coefficients(typeof(W), :c)
                        push!(ch, c)
                    end
                end
                d17 = generic_quasihomogeneous_polynomial(deg_d17, gradings, coeffs)
                d35 = generic_quasihomogeneous_polynomial(deg_d35, gradings, coeffs)
                d15 = generic_quasihomogeneous_polynomial(deg_d15, gradings, coeffs)
                d26 = generic_quasihomogeneous_polynomial(deg_d26, gradings, coeffs)
                d16 = generic_quasihomogeneous_polynomial(deg_d16, gradings, coeffs)
                d25 = generic_quasihomogeneous_polynomial(deg_d25, gradings, coeffs)

                equation = -d17*d35 + d15*d26 - d16*d25 - W
                eqns = coefficients(equation, vars...)

                if !(1 in eqns)
                    a = maximum(max_variable_index, eqns)
                    dense_eqns = to_dense_monomials(eqns)
                    @show dense_eqns
                    if !(rem(one(eltype(dense_eqns)), groebner_basis(dense_eqns)) == 0)
                        count += 1
                    end
                end
            end
        end
    end

    info("Found $count possibilities")
end

function E14_E14()
    @ring! ℚ[a2,a3,a4,a5,
             b22,b23,b24,
             c1,
             x,y,z,u,v,w]

    # ARC: conditions on the parameters -----------------------
    b21=-1//2*(-2*b22*c1+2*b24*c1^2-2*b23*c1^3+c1^4+2*a5*c1^4)
    a1=b21+a2*c1-b22*c1-a4*c1^2+b24*c1^2+a3*c1^3-b23*c1^3+c1^4
    d11=a1*a5+a5*b21+a3*b22+a2*b23+a4*b24-a2*a5*c1-a5*b22*c1-a4*b23*c1-a3*b24*c1+a4*a5*c1^2+a3*b23*c1^2+a5*b24*c1^2-a3*a5*c1^3-a5*b23*c1^3+a5^2*c1^4
    b25=a5
    d14=a5*b25
    d155=a5*b23+a3*b25-c1*d14
    d13=a3*b23+a5*b24+a4*b25-c1*d155
    d12=a5*b22+a4*b23+a3*b24+a2*b25-c1*d13
    d10=a5+b25
    d8=a3*b21+a4*b22+a1*b23+a2*b24-c1*d11
    d7=a3+b23-c1*d10
    d5=a4*b21+a2*b22+a1*b24-a3*b21*c1-a4*b22*c1-a1*b23*c1-a2*b24*c1+a1*a5*c1^2+a5*b21*c1^2+a3*b22*c1^2+a2*b23*c1^2+a4*b24*c1^2-a2*a5*c1^3-a5*b22*c1^3-a4*b23*c1^3-a3*b24*c1^3+a4*a5*c1^4+a3*b23*c1^4+a5*b24*c1^4-a3*a5*c1^5-a5*b23*c1^5+a5^2*c1^6
    d2=a2*b21+a1*b22-c1*d5
    d4=a4+b24-a3*c1-b23*c1+2*a5*c1^2
    d1=a2+b22-c1*d4
    d9=-a3+b23-c1
    d6=-a4+b24-c1*d9
    d3=-a2+b22-c1*d6
    # ----------------------------------------------------------

    # ARC: In this case, we used directly the structure of the ansatz for a matrix factorization, which makes us have less parameters to take care of
    Q = StandardDuality(
        d15=z+w+a1*u^4+a2*x*u^3+a3*x^3*u+a4*x^2*u^2+a5*x^4,
        d16=y^2+y*v+v^2,
        d17=x^3*z+d1*u^3*w+d2*u^7+d3*z*u^3+d4*x*u^2*w+d5*x*u^6+d6*x*z*u^2+d7*x^2*u*w+d8*x^2*u^5+d9*x^2*z*u+d10*x^3*w+d11*x^3*u^4+d12*x^4*u^3+d13*x^5*u^2+d14*x^7+d155*x^6*u,
        d25=y-v,
        d26=-z+w+b21*u^4+b22*x*u^3+b23*x^3*u+b24*x^2*u^2+b25*x^4,
        d35=x+c1*u,
    )
    # ----------------------------------------------------------

    f1, f2 = (c1^4 - 2c1^2 + 2)//2, (c1^4 + 2c1^2 + 2)//2

    Q, f1, f2
end

function U12_1_U12_3()
    @ring! ℚ[a1,a2,b1,b2,c1,c2,x,y,z,u,v,w]

    Q = StandardDuality(
        d15 = z + a1*v + a2*w,
        d16 = -y*z + (-a1^2 + a1*b1 + a1*c1)*v^2 + (-2a1*a2 + b1*a2 + c1*a2 + a1*b2 + a1*c2)*v*w +
              + a2*(-a2+b2+c2)*w^2 + c1*v*z + c2*w*z,
        d17 = u^3 + u^2*x + u*x^2 + x^3,
        d25 = -y + b1*v + b2*w,
        d26 = -y*z + b1*c1*v^2 + (c1*b2 + b1*c2)*v*w + b2*c2*w^2 - (-a1+b1+c1)*v*y -
              (-a2+b2+c2)*w*y,
        d35 = x - u,
    )

    f1 = -1 + a1^2*b1 - a1*b1^2
    f2 = 2a1*b1*a2 - b1^2*a2 + a1^2*b2 - 2a1*b1*b2
    f3 = b1*a2^2 + 2a1*a2*b2 - 2b1*a2*b2 - a1*b2^2
    f4 = -1 + a2^2*b2 - a2*b2^2

    Q, f1, f2, f3, f4
end

function W12_W12()
    @ring! ℚ[a,b,c,x,y,z,u,v,w]

    f1 = a^2-a*b
    f2 = -a*b^3+4*a^3*c+4*a*b^2*c-5*a*b*c^2+2*a*c^3+a^2*(b^2-6*b*c+5*c^2)+1//4*(4+b^4-4*b^3*c+6*b^2*c^2-4*b*c^3+c^4)

    Q1 = [w+a*u*x+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2)*x^2+z u*w+c*w*x+(1//2*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)+a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3+b*x*z v^4+v^3*y+v^2*y^2+v*y^3+y^4 0
          u+(-2*a+b-c)*x -w+(-a+b)*u*x+(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2))*x^2+z 0 v^4+v^3*y+v^2*y^2+v*y^3+y^4
          v-y 0 -w+(-a+b)*u*x+(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2))*x^2+z -u*w-c*w*x+(-(1//2)*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)-a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3-b*x*z
          0 v-y -u+(2*a-b+c)*x w+a*u*x+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2)*x^2+z
         ]

    Q2 = [w+(a-b)*u*x+(-b*(-2*a+b-c)+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2))*x^2-z u*w+c*w*x+(1//2*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)+a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3+b*x*z v^4+v^3*y+v^2*y^2+v*y^3+y^4 0
          u+(-2*a+b-c)*x -w-a*u*x+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)*x^2-z 0 v^4+v^3*y+v^2*y^2+v*y^3+y^4
          v-y 0 -w-a*u*x+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)*x^2-z -u*w-c*w*x+(-(1//2)*(-a+b)*(-2*a*b+b^2+2*a*c-2*b*c+c^2)-a*(b*(-2*a+b-c)+1//2*(2*a*b-b^2-2*a*c+2*b*c-c^2)))*x^3-b*x*z
          0 v-y -u+(2*a-b+c)*x w+(a-b)*u*x+(-b*(-2*a+b-c)+1//2*(-2*a*b+b^2+2*a*c-2*b*c+c^2))*x^2-z
         ]

    z = zeros(Q1)

    Q = [z Q1; Q2 z]

    Q, f1, f2
end


end
