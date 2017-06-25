module MatrixFactorizations

using PolynomialRings: polynomial_ring

export StrangeDuality, S1, T2, Dk_1, E6_1, E7_3, one_certain_unit_mf, E14_Q10

StrangeDuality() = begin
    A, (x,y,z) = polynomial_ring(Int, :x, :y, :z)

    d1 = [z   y^2 x^3 0 ; y -x*z 0 x^3 ; x 0 -x*z -y^2 ; 0 x -y z  ]
    d0 = [x*z y^2 x^3 0 ; y -z   0 x^3 ; x 0 -z   -y^2 ; 0 x -y x*z]
    zz = zero(d1)

    [ zz d1; d0 zz ]
end

S1() = begin
    A, (x,y) = polynomial_ring(Int, :x, :y)

    q1=[x x*y; -y -x^3]
    q0=[x^3 x*y; -y -x]
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

T2() = begin
    A, (x,y) = polynomial_ring(Int, :x, :y)

    q1 = [x^2 y ; -y -x^2]
    q0 = [x^3 x*y; -x*y -x^3]
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

# The following are the factorizations from [22]

Dk_1(k::Integer, l::Integer) = begin
    A, (x,y,z) = polynomial_ring(Int, :x, :y, :z)

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
    A, (x,y,z) = polynomial_ring(Int, :x, :y, :z)

    q1 = [-z 0 x^2 y^3; 0 -z y -x; x y^3 z 0; y -x^2 0 z]
    q0 = q1
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

E7_3() = begin
    A, (x,y,z) = polynomial_ring(Rational{Int}, :x, :y, :z)

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
    A, (x,y,z,u,v,w) = polynomial_ring(Rational{BigInt}, :x, :y, :z, :u, :v, :w)

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

E14_Q10() = begin A, (x,y,z,u,v,w) = polynomial_ring(Rational{Int}, :x, :y, :z, :u, :v, :w)
    return x^4 +y^3 + x*z^2 - u^4*w -v^3 - w^2
end

unit_matrix_factorization(f, left_vars, right_vars) = begin



end

using PolynomialRings
using QuasiHomogeneous: generic_matrix_factorizations
function E14_Q10_possibilities()
    R, vars = polynomial_ring(Int,:x,:y,:z,:u,:v,:w)
    F = formal_coefficients(R,:c)
    next_coeff() = take!(F)

    generic_matrix_factorizations(4,7,7,24,(6,8,9,3,8,12),R,next_coeff)
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

end
