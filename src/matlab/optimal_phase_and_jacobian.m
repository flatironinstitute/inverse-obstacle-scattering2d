function [z,jacwphase,jacbarwphase] = optimal_phase_and_jacobian(u,f,jac,jacbar)
%%OPTIMAL_PHASE_AND_JACOBIAN
%
% For the non-linear least squares problem 
%
% min_{beta,|z|=1} \| u - z f(beta) \|_2^2
%
% Where the jacobian of f with respect to beta is provided, this routine
% returns the optimal z and the jacobian of 
%
% h = z(beta) f(beta) 
%
% where z(beta) is the optimal z with |z| = 1 for beta fixed.
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
% jacwphase - the jacobian of h as described above (only if jac provided, 
%                   otherwise empty)
%

fdu = f'*u;
fdunrm = abs(fdu);

z = fdu/fdunrm;

if nargin < 3
    jacwphase = [];
    return
end

if nargin < 4
    jacbar = zeros(size(jac));
end

jacu = jacbar'*u;
jacbaru = jac'*u;
jacabs = 0.5*(conj(fdu)*(jacbar'*u)+fdu*conj(jac'*u))/fdunrm;
jacbarabs = 0.5*(conj(fdu)*(jac'*u)+fdu*conj(jacbar'*u))/fdunrm;

jacwphase = jac*z + f*( (jacu/fdunrm - fdu*jacabs/fdunrm^2).' );
jacbarwphase = jacbar*z + f*( (jacbaru/fdunrm - fdu*jacbarabs/fdunrm^2).' );
