clear
load('coal_mine.mat')

%% data visualization
x = linspace(1851, 1963, length(coal_mine));
%histogram(coal_mine);

%% defining parameters
hyperParam = 1;
breakpoints = 3;
d = breakpoints+1;
ro = 0.1;
acc = 0;

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

% theta parameter
theta = zeros(1,M);

% number of disasters
nDis = zeros(d,M);

% f(t)
ft = @(lambda,t,n) prod(lambda.^n .* exp(-lambda.*diff(t)) .* diff(t));

for j = 1:M-1
    
    theta(j) = gamrnd(2,1/hyperParam);
    
    % number of disaster update
    for i = 1:d
        nDis(i,j) = sum(coal_mine >= t(i,j) & coal_mine <= t(i+1,j));        
    end
    
    % lambda update
    lambda(:,j) = gamrnd(nDis(:,j)+2, 1./(diff(t(:,j))+theta(j)));
    
    % MH step
    for l = 2:d                
        cand = randWalkProp(ro, [t(l-1,j) t(l,j) t(l+1,j)]);
        t_star = t(:,j);
        t_star(l) = cand;
       
        if cand > t(l-1,j) && cand < t(l+1,j);
            for i = 1:d
                nDis_star(i,1) = sum(coal_mine >= t_star(i) & ...
                    coal_mine <= t_star(i+1));
            end
            alpha = ft(lambda(:,j), t_star, nDis_star) / ...
                    ft(lambda(:,j), t(:,j), nDis(:,j));
            if rand <= alpha
                t(l,j+1) = cand;
                acc = acc + 1;
            else
                t(l,j+1) = t(l,j);
            end            
        else
          t(l,j+1) = t(l,j);
        end
    end  
        
end 

tau = ceil(mean(t(:, burn_in:M),2));

%% displaying lambda intensities

figure
s = size(lambda,1);
for i = 1:s
    subplot(s,1,i)
    h = histfit(lambda(i,:),50,'gamma');
    h(1).FaceColor = [.6 .8 1];
    str = sprintf('lambda %d', i);
    title(str) 
end    

%% displaying theta intensities
figure
h = histfit(theta,50,'gamma');
h(1).FaceColor = [.6 .8 1];

%% displaying t distribution
figure
%t_dist = prod(t(3:6,:));
t_dist = prod(t);
histogram(t_dist);

%% statistics
% to change
if size(t,1) == 3
    histogram(cand)
    str=sprintf('acceptance rate = %.3f, tau = %d', acc/M, tau);
    title(str)
end