%--------------------------------------------------------------------------
% File: class1.m
%
% Goal: Compute the coefficients a(j), b(j) of the recurrence relation
%
%             b p (x) = (x - a ) p   (x) - b   p   (x)
%              j j            j   j-1       j-1 j-2
%
%        for the various classical (normalized) orthogonal polynomials,
%        and the zero-th moment muzero = integral w(x) dx of the given
%        polynomial's weight function w(x). Since the polynomials are
%        orthonormalized, the tridiagonal matrix is guaranteed to be
%        symmetric.
%
% Use: [a,b,muzero] = class1(kind,n,alpha,beta)
%
% Input: kind - kind of classical (normalized) orthogonal polynomials
%        n - number of coefficients
%        alpha - parameter used for laguerre and jacobi polynomials
%        beta - parameter used for jacobi polynomials
%
% Output: a - coefficients
%         b - coefficients
%         muzero - the zero-th moment
%
% Recalls:gamma.m
%
% Authors:
% Date last modified: March, 2024
%--------------------------------------------------------------------------
function [a,b,muzero] = class1(kind,n,alpha,beta)
nm1 = n - 1 ;

% kind = 1: legendre polynomials p(x) on (-1, +1), w(x) = 1.
if kind==1
    muzero = 2.0;
    for   i = 1: nm1
        a(i) = 0.0;
        abi = i;
        b(i) = abi/sqrt(4*abi*abi - 1.0);
    end % for
    a(n) = 0.0;

    % kind = 2: chebyshev polynomials of the first kind t(x) on (-1, +1),
    % w(x) = 1 / sqrt(1 - x*x)
elseif kind==2
    muzero = pi;
    for   i = 1: nm1
        a(i) = 0.0;
        b(i) = 0.5;
    end % for
    b(1) = sqrt(0.5) ;
    a(n) = 0.0;

    % kind = 3: chebyshev polynomials of the second kind u(x) on (-1, +1),
    % w(x) = sqrt(1 - x*x)
elseif kind==3
    muzero = pi/2.0;
    for  i = 1 : nm1
        a(i) = 0.0d0 ;
        b(i) = 0.5d0;
    end % for i
    a(n) = 0.0d0;

    % kind = 4: hermite polynomials h(x) on (-infinity,+infinity),
    % w(x) = exp(-x**2)
elseif kind==4
    muzero = sqrt(pi);
    for i = 1 : nm1
        a(i) = 0.0;
        b(i) = sqrt(i/2);
    end % for i
    a(n) = 0.0;

    % kind = 5:  jacobi polynomials p(alpha, beta)(x) on (-1, +1),
    % w(x) = (1-x)**alpha + (1+x)**beta, alpha and beta greater than -1
elseif kind==5
    ab = alpha + beta;
    abi = 2.0 + ab;
    muzero = 2.0^(ab+1.)*gamma(alpha+1)*gamma(beta+1) / gamma(abi);
    a(1) = (beta - alpha)/abi;
    b(1) = sqrt(4.0*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi));
    a2b2 = beta*beta - alpha*alpha;
    for  i = 2: nm1
        abi = 2.0*i + ab;
        a(i) = a2b2/((abi - 2.0)*abi);
        b(i) = sqrt(4*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi));
    end % for
    abi = 2.0*n + ab;
    a(n) = a2b2/((abi - 2.0)*abi);

    % kind = 6:  laguerre polynomials l(alpha)(x) on (0, +infinity),
    % w(x) = exp(-x) * x**alpha, alpha greater than -1.
else
    muzero = gamma(alpha + 1.0);
    for   i = 1 :nm1
        a(i) = 2.0*i - 1.0 + alpha;
        b(i) = sqrt(i*(i + alpha));
    end % for
    a(n) = 2.0*n - 1 + alpha;
end % if

