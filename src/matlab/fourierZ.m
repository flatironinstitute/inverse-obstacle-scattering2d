function z = fourierZ(zhat,t)     % must work on vector of t's
t = 2*pi*t;
N =numel(zhat);  % even
z = 0*t;
for k=0:N/2
  z = z + zhat(k+1)*exp(1i*k*t);
end
for k=-N/2:-1
  z = z + zhat(k+1+N)*exp(1i*k*t);
end
end


