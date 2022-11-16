

ks = randi(20,6,1);
fun = @(thetas) exp(1i*ks*thetas(1)) + ks*thetas(2)^2;
jac = @(thetas) [1i*ks.*exp(1i*ks*thetas(1)), 2*thetas(2)*ks];

phase = exp(1i*randn());
thetatrue = randn(2,1) + 1i*randn(2,1);
u = fun(thetatrue)*phase;

[phase2,~] = optimal_phase_and_jacobian(u,fun(thetatrue),jac(thetatrue));

assert(abs(phase2-phase) < 1e-14);

thetatest = randn(2,1) + 1i*randn(2,1);
[f0,j0,j0bar] = fun2(thetatest,u,fun,jac); 

ntest = 10;

dir = randn(2,1) + 1i*randn(2,1);

for i = 1:ntest
    delta = 5^(-i)*dir;
    thetas2 = thetatest+delta;
    [f1,j1,j1bar] = fun2(thetas2,u,fun,jac);
    
    diff = (j0+j1)*delta + (j0bar+j1bar)*conj(delta);
    
    err = 1 - 0.5*diff./(f1-f0);
    fprintf('step %5.2e err %5.2e\n',norm(delta,Inf),norm(err,Inf));
end

ks = randi(20,6,1);
fun = @(thetas) exp(1i*ks*thetas(1)) + ks*thetas(2)^2;
jac = @(thetas) [1i*ks.*exp(1i*ks*thetas(1)), 2*thetas(2)*ks];

const = rand()*exp(1i*randn());
thetatrue = randn(2,1) + 1i*randn(2,1);
u = fun(thetatrue)*const;

[const2,~] = optimal_const_and_jacobian(u,fun(thetatrue),jac(thetatrue));

assert(abs(const2-const) < 1e-14);

thetatest = randn(2,1) + 1i*randn(2,1);
[f0,j0,j0bar] = fun3(thetatest,u,fun,jac); 

ntest = 10;

dir = randn(2,1) + 1i*randn(2,1);

for i = 1:ntest
    delta = 5^(-i)*dir;
    thetas2 = thetatest+delta;
    [f1,j1,j1bar] = fun3(thetas2,u,fun,jac);
    
    diff = (j0+j1)*delta + (j0bar+j1bar)*conj(delta);
    err = 1 - 0.5*diff./(f1-f0);
    fprintf('step %5.2e err %5.2e\n',norm(delta,Inf),norm(err,Inf));
end

ks = randi(20,6,1);
fun = @(thetas) exp(1i*ks*thetas(1)) + ks*thetas(2)^2;
jac = @(thetas) [1i*ks.*exp(1i*ks*thetas(1)), 2*thetas(2)*ks];

const = exp(1i*randn());
thetatrue = randn(2,1);
u = fun(thetatrue)*const;

[const2,~] = optimal_const_and_jacobian(u,fun(thetatrue),jac(thetatrue));

assert(abs(const2-const) < 1e-14);

thetatest = randn(2,1);
[f0,j0,j0bar] = fun2(thetatest,u,fun,jac); 

ntest = 10;

dir = randn(2,1);
fprintf('test (real delta)\n');

for i = 1:ntest
    delta = 5^(-i)*dir;
    thetas2 = thetatest+delta;
    [f1,j1,j1bar] = fun2(thetas2,u,fun,jac);
    
    diff = (j0+j1)*delta + (j0bar+j1bar)*conj(delta);
    err = 1 - 0.5*diff./(f1-f0);
    fprintf('step %5.2e err %5.2e\n',norm(delta,Inf),norm(err,Inf));
end


const = exp(1i*randn());
thetatrue = randn(2,1);
u = fun(thetatrue)*const;

thetas2=randn(size(thetatrue));
[f0,j0,j0bar] = fun2(thetas2,u,fun,jac);
norm(u-f0)
ntimes = 100;
reg=1;
for i = 1:ntimes
    jj = [real(j0+j0bar) imag(j0bar-j0); imag(j0+j0bar) real(j0-j0bar); reg*eye(2*length(thetas2))];
    jj2 = [real(j0+j0bar) imag(j0bar-j0); imag(j0+j0bar) real(j0-j0bar); reg/2*eye(2*length(thetas2))];
    rhs = [real(u-f0);imag(u-f0);zeros(2*length(thetas2),1)];
    delta = jj\rhs; delta = delta(1:end/2)+1i*delta(end/2+1:end);
    delta2 = jj\rhs; delta2 = delta2(1:end/2)+1i*delta2(end/2+1:end);
    thetas3 = thetas2 + delta(:);
    thetas4 = thetas2 + delta2(:);
    [f3,j3,j3bar] = fun2(thetas3,u,fun,jac);
    [f4,j4,j4bar] = fun2(thetas4,u,fun,jac);
    if (norm(u-f3) < norm(u-f0))
        f0 = f3;
        thetas2= thetas3;
        j0=j3;
        j0bar=j3bar;
        if (norm(u-f4) < norm(u-f3))
            reg = reg/2;
            f0 = f4;
            thetas2= thetas4;
            j0=j4;
            j0bar=j4bar;
        end
    else
        reg = reg*2;
    end
end
norm(u-f0)

function [f,j,jbar] = fun2(thetas,u,fun,jac)

f0 = fun(thetas);
j0 = jac(thetas);

[z,j,jbar] = optimal_phase_and_jacobian(u,f0,j0);
f = f0*z;

end

function [f,j,jbar] = fun3(thetas,u,fun,jac)

f0 = fun(thetas);
j0 = jac(thetas);

[z,j,jbar] = optimal_const_and_jacobian(u,f0,j0);
f = f0*z;

end