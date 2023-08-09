function [ckcoefs,jac] = constkappa_models_convert(coefs,impedance_type)
%CONSTKAPPA_MODELS_CONVERT helper function for impedance models of the
% form 
%
%    lambda = c1 + c2*H
%
% where H is the signed curvature. The idea is that, given coefficients
% b_1,b_2, .. b_m for a given model, these should be converted to c1 and
% c2.
%
% input: 
%
% coefs - array of floats, coefficients in the chosen model
% impedance_type - string determining the model type 
%
%        impedance_type = 'constkappa' is the basic model, with 
%                     c1 = b1
%                     c2 = b2
%        impedance_type = 'antbar2' is a simpler model of Antoine-Barucq 
%          type, with 2 coefficients where 
%                     c1 = b2/sqrt(1+1i*b1)
%                     c2 = 1/(1+1i*b1)
%        impedance_type = 'antbar3' is the full model of Antoine-Barucq 
%          type, with 3 coefficients where 
%                     c1 = b2*b3/sqrt(1+1i*b1)
%                     c2 = b3/(1+1i*b1)
%
% output:
%
% ckcoefs - the coefficients c1 and c2 as given above 
% jac - the 2 x length(coefs) matrix in which row 1 is the gradient of c1
%          with respect to the bi and row 2 is the gradient of c2
%

ckcoefs = zeros(2,1);
jac = zeros(2,length(coefs));

if strcmpi(impedance_type,'constkappa')
    ckcoefs = coefs;
    jac = eye(2);
elseif strcmpi(impedance_type,'antbar2')
    aa = 1+1i*coefs(1);
    aasqrt = sqrt(aa);
    ckcoefs(1) = coefs(2)/aasqrt;
    ckcoefs(2) = 1.0/aa;
    jac(1,1) = -0.5*1i*coefs(2)/(aa*aasqrt);
    jac(1,2) = 1.0/aasqrt;
    jac(2,1) = -1.0*1i/(aa*aa);
    jac(2,2) = 0;
elseif strcmpi(impedance_type,'antbar3')
    aa = 1+1i*coefs(1);
    aasqrt = sqrt(aa);
    ckcoefs(1) = coefs(2)*coefs(3)/aasqrt;
    ckcoefs(2) = coefs(3)/aa;
    jac(1,1) = -0.5*1i*coefs(2)*coefs(3)/(aa*aasqrt);
    jac(1,2) = coefs(3)/aasqrt;
    jac(1,3) = coefs(2)/aasqrt;
    jac(2,1) = -coefs(3)*1i/(aa*aa);
    jac(2,2) = 0;
    jac(2,3) = 1/aa;
else
    error('unknown model selected');
end
