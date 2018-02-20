hilbert(n,k) = n == 1 ? 1 : div( prod(k+i for i=1:n-1), factorial(n-1))
cumhilbert(n,k) = k < 0 ? 0 : sum(hilbert(n,j) for j=0:k)

function degrevlex_index(exponents)
   total_degree = sum(exponents)
   ret = cumhilbert(length(exponents), total_degree-1) + 1
   for (i,e) in enumerate(exponents[1:end-1])
       if e != 0
           ret += cumhilbert(length(exponents) - i, e)-1
       end
   end
   return ret
end

[degrevlex_index(t) for t in [(0,0), (0,1), (1,0), (0,2), (1,1), (2,0), (0,3), (1,2)]]

[degrevlex_index(t) for t in [(0,0,0), (0,0,1), (0,1,0), (1,0,0), (0,0,2), (0,1,1), (1,0,1), (0,2,0), (1,1,0), (2,0,0), (0,0,3)]]
