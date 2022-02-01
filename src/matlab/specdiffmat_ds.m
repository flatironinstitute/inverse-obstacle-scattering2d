function [D] = specdiffmat_ds(N,ds)
    tj = 2*pi/N*(0:N-1);
    D = (-1).^(1:N) .* cot(tj/2) / 2;   % note overall - sgn due to 1st row not col
    D(1) = 0;              % kill the Inf
    D = circulant(D);
    D = diag(1.0./ds')*D;
end