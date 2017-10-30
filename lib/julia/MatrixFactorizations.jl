module MatrixFactorizations

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

end
