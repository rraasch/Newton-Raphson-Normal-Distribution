%Newton-Raphson method for N(2,1)

%In class:
%Generate a sample of size 1000 from N(2, 1).
%Use the Newton-Raphson method to estimate mu and sigma^2

%set seed
rng(1)

%random sample of n=1000 random numbers from N(2, 1)
n = 100000;
x = normrnd(2, 1, n, 1);

%set tolerance
tol = 0.000001;

%set initial values:
mu_0 = 3.5;
sigma2_0 = 0.1;
theta_0 = [mu_0; sigma2_0];

%define the function U_0 and I
ybar = mean(x);
U_0 = [n/sigma2_0*(ybar - mu_0); -n/(2*sigma2_0)+1/(2*sigma2_0^2)*sum((x-mu_0).^2)];
I_0 = [n/sigma2_0, 0; 0, n/(2*sigma2_0^2)];

%calculate first iteration
thetahat = theta_0 + inv(I_0)*U_0;

%count the iterations;
iter = 1;

while abs(thetahat(1)-theta_0(1)) > tol | abs(thetahat(2)-theta_0(2)) > tol
    theta_0 = thetahat;
    mu_0 = theta_0(1);
    sigma2_0 = theta_0(2);
    U_0 = [n/sigma2_0*(ybar - mu_0); -n/(2*sigma2_0)+1/(2*sigma2_0^2)*sum((x-mu_0).^2)];
    I_0 = [n/sigma2_0, 0; 0, n/(2*sigma2_0^2)];
    thetahat = theta_0 + inv(I_0)*U_0;
    iter = iter+1;
end

display([newline 'mu is ', num2str(thetahat(1))])
display([newline 'sigma^2 is ', num2str(thetahat(2))])
display([newline 'number of iterations is ', num2str(iter), newline])
