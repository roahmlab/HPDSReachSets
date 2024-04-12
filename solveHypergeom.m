function c = solveHypergeom(k, btilde, lambda, alpha, t)

% k is order
% btitle is entry of V^T*b
% lambda is Z-eigenvalue
% alpha is initial condition for c(t)
% t is time


a = (k-2)/(k-1);
z = -btilde/(lambda*alpha^(k-1));

fun =@(c) -hypergeom([1, a], a+1, -btilde/(lambda*c^(k-1)))/((k-2)*lambda*c^(k-2))+hypergeom([1, a], a+1, z)/((k-2)*lambda*alpha^(k-2))-t;

x0 = alpha;
an = fsolve(fun, x0, optimoptions('fsolve','Display','off','FunctionTolerance',1e-14,'StepTolerance',1e-14));

if sign(alpha) == -sign(an)
                error('Error: Function crosses 0');
end

c = an;
