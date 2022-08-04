% This script reproduces Figs [3,4] in Sec. 4 Illustrative examples
% Author: Juan Felipe Florez-Ospina
% Last updated: August 04 - 2021

clearvars
close all
n = 256;

% Generate signal of interest 
signal_idx = 1;
signal_names = {'Blocks','Doppler','Random'};
switch signal_names{signal_idx}
    case 'Blocks'      
        rng(0)
        x0 = 0.75*reshape(MakeSignal('Blocks',n),n,1);
        x = x0 + 0.15*randn(n,1);
%     case 'Doppler'
%         x0 = reshape(MakeSignal('Doppler',512),n,1);
%     case 'Random'
%         x0 = 0.5*randn(n,1);
end

% Generate m undersampled noisy linear measurements
m = round(0.35*n);
rng(3);
A = randn(m,n);
y = A*x;
sigma = 0.1*(norm(y)/m);
noise = sigma*randn(m,1);
y = y + noise; 

% Define available rank-order information
p = zeros(3,n);
switch signal_names{signal_idx}
    case 'Blocks'    
        p(1,:) = 1:n;
        [~,p(2,:)] = sort(x0, 'ascend');
        [~,p(3,:)] = sort(x, 'ascend');
end


% Construct adjacency matrix of path graph induced by temporal ranks
one = ones(n,1);
W_G0 = spdiags([one one], [-1 1], n, n);
% Construct graph Laplacian
L_G0 = diag(sum(W_G0,2))-W_G0;

legends = {'Path Graph on $\{i\}_{i=1}^n$',...
           'Path Graph on $\{\hat{r}_i^{-1}\}_{i=1}^n$',...
           'Path Graph on $\{r_i^{-1}\}_{i=1}^n$'};
colors = {'b','r','k'};
fig0_legends = {'$x_i$','$x_{\hat{r}_i^{-1}}$','$x_{r^{-1}_i}$'};

for i = 1 : size(p,1)
    % Construct permutation matrix associated with above permutations
    P = eye(n,n); 
    P = P(p(i,:),:);

    figure(1)
    plot(P*x,colors{i},'LineWidth',1)
    hold on
    xlim([1 n])
    ylabel('Amplitude','Interpreter','latex','FontSize',15)
    xlabel('$i$','Interpreter','latex','FontSize',15)

    % Construct associated graph Laplacian
    L_G = P'*(L_G0*P);
    
    % Reconstruct x from y using L_G
    epsilon = sigma*sqrt((m + 1*sqrt(2*m)));
    cvx_begin
	    variables x1(n)
        minimize(x1'*(L_G*x1))
	    subject to
	    norm(A*x1-y) <= epsilon;
    cvx_end
    psnr_dB(i) = 10*log10(max(x)^2/mean((x1-x).^2));

    figure(i+1)
    axis normal    
    plot(x1,'k','LineWidth',1), hold on
    plot(x1-x,'r','LineWidth',1), 
    xlim([0 n])
    ylim([min(x1) 6])
    xlabel('i','Interpreter','latex')
    ylabel('Amplitude Estimate','Interpreter','latex')
    hold on
    leg{i} = {legends{i},strcat('Error signal,','PSNR=',num2str(psnr_dB(i)))};
end

figure(size(p,1)+2)
axis normal
plot(y,'k','LineWidth',1), hold on
xlim([0 m])
ylim([min(y) max(y)])
xlabel('i','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')
hold on
leg{size(p,1)+1} = {'Noisy Undersampled Measurements'};

figs_path = './figs\';
for p = 1:5
    ax = get(figure(p),'CurrentAxes');
    set(ax,'TickLabelInterpreter','latex')
    set(ax, 'FontSize', 20)    
    if p ~= 1
        leg_ax = legend(ax,leg{p-1});
    else
        leg_ax = legend(ax,fig0_legends);
    end
    set(leg_ax,'Interpreter','latex')
    set(leg_ax,'FontSize',15)    
    set(leg_ax,'Location','Best')   
    print(figure(p),strcat(figs_path,'Fig_IllustrativeExamples',num2str(p)),'-depsc')    
end

