function zppp = fourierZpp(zhat,t) % deriv func Z''
N = numel(zhat);
zppp = 2*pi*fourierZpp(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);
end

