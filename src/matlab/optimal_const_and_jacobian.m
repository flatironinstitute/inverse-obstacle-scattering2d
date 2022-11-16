function [z,jacwconst,jacbarwconst] = optimal_const_and_jacobian(u,f,jac,jacbar)
%OPTIMAL_CONST_AND_JACOBIAN
%
% For the non-linear least squares problem 
%
% min_{beta,z} \| u - z f(beta) \|_2^2
%
% Where the jacobian of f with respect to beta is provided, this routine
% returns the optimal z and the jacobian of 
%
% h = z(beta) f(beta) 
%
% where z(beta) is the optimal z for beta fixed.
%
% Input:
%
% u - vector-valued data 
% f - values of the function f at beta 
% jac - jacobian of f with respect to beta at the particular value of beta
% (optional)
%
% Output:
%
% z - the optimal z for the given u and f values
% jacwconst - the jacobian of h as described above (only if jac provided, 
%                   otherwise empty)
%

fdu = f'*u;
ff = f'*f;

z = fdu/ff;

if nargin < 3
    jacwconst = [];
    return
end

if nargin < 4
    jacbar = zeros(size(jac));
end

jacu = jacbar'*u;
jacbaru = jac'*u;
jacf = jacbar'*f + (f'*jac).';
jacbarf = jac'*f + (f'*jacbar).';
jacwconst = z*jac + f*(((ff*jacu - fdu*jacf)/ff^2).');
jacbarwconst = z*jacbar + f*((jacbaru/ff - fdu*jacbarf/ff^2).');
