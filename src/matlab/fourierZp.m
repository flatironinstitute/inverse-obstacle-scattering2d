function zp = fourierZp(zhat,t)  % deriv func Z'
N = numel(zhat);
zp = 2*pi*fourierZ(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);
end

