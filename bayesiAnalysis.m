load('coal_mine.mat')

%% data visualization
x = linspace(1851, 1963, length(coal_mine));
%histogram(coal_mine);

%% defining parameters
hyperParam = 1;
breakpoints = 1;
d = breakpoints+1;
theta = gamrnd(2,1/hyperParam);

%% gibbs sampling with MH step

N = 10000;
burn_in = 2000;
M = N + burn_in;

% time intervals
t = zeros(d+1,M);
t(1,:) = 1851;
t(d+1,:) = 1963;

% first breakpoint
t(:,1) = linspace(1851, 1963, d+1);

% lambda intensities
lambda = zeros(d,M);

% number of disasters
nDis = zeros(d,M);

% f(t)
ft = @(lambda,t,n) prod(lambda.^n .* exp(-lambda.*diff(t)) .* diff(t));
fun = @(x,y,z,n) (x.^n).*(y-z)*exp(-x.*(y-z));
%prod(lambda(:,j).^nDis(:,j).*exp(-lambda(:,j).*diff(t(:,j))).*diff(t(:,j)))

ro = 0.1;
acc = 0;

for j = 1:M-1
    for i = 1:d
        nDis(i,j) = sum(coal_mine >= t(i,j) & coal_mine <= t(i+1,j));
        %lambda(i,j) = gamrnd(nDis(i,j)+2, 1/(t(i+1,j)-t(i,j)+theta));
    end
    lambda(:,j) = gamrnd(nDis(:,j)+2, 1./(diff(t(:,j))+theta));
    
    cand(j) = randWalkProp(ro, t(:,j));
    t_star = t(:,j);
    t_star(2) = cand(j);
    
    if cand(j) > t(1,j) && cand(j) < t(3,j)
        for i = 1:d
            nDis_star(i,1) = sum(coal_mine >= t_star(i) & coal_mine <= t_star(i+1));
        end
        num = ft(lambda(:,j), t_star, nDis_star);
        den = ft(lambda(:,j), t(:,j), nDis(:,j));
        alpha = num / den;
        r = rand;
        if r <= alpha
            t(2,j+1) = cand(j);
            acc = acc + 1;
        else
            t(2,j+1) = t(2,j);
        end
    end        
    
    %{
    num = fun(lambda(2,j), t(3,j), cand(j), nDis(2,j));
    den = fun(lambda(2,j), t(3,j), t(2,j), nDis(2,j));
    alpha = num / den;
    %}        
    
    theta = gamrnd(2*(d+1), 1/(hyperParam+sum(lambda(:,j))));
end 

tau = ceil(mean(t(2, burn_in:M)));

%% statistics
histogram(cand)
str=sprintf('acceptance rate = %.3f, tau = %d', acc/M, tau);
title(str)
