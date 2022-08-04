% This script reproduces Fig 5 in Sec. 4 Illustrative examples
% Author: Juan Felipe Florez-Ospina
% Last updated: August 04 - 2021

close all

clearvars
n = 256;

% Construct DCT basis
for k = 0 : n-1
    for j = 0 : n-1    
        U(j+1,k+1) = cos((j+1/2)*pi*k/n);    
    end
end
U(:,1) = U(:,1)/sqrt(2);

% Generate signal of interest 
signal_idx = 1;
signal_names = {'Blocks','Doppler','Random'};
switch signal_names{signal_idx}
    case 'Blocks'      
        rng(0)
        x0 = 0.75*reshape(MakeSignal('Blocks',n),n,1);
        x = x0 + 0.15*randn(n,1);
end

% Define available rank-order information
p = zeros(3,n);
switch signal_names{signal_idx}
    case 'Blocks'    
        p(1,:) = 1:n;
        [~,p(2,:)] = sort(x0, 'ascend');
        [~,p(3,:)] = sort(x, 'ascend');
end


figure(1),
stem(U'*x,'LineWidth',1)
xlim([1 n/4])
ylabel('$\langle\mathbf{P}^Tu_k,x\rangle$','Interpreter','latex','FontSize',15)
xlabel('$k$','Interpreter','latex','FontSize',15)

figure(2),
P = eye(n,n); 
P = P(p(2,:),:);
stem((P'*U)'*x,'LineWidth',1)
xlim([1 n/4])
xlim([1 n/4])
ylabel('$\langle\mathbf{P}^Tu_k,x\rangle$','Interpreter','latex','FontSize',15)
xlabel('$k$','Interpreter','latex','FontSize',15)
% leg{3} = {};

figure(3),
P = eye(n,n); 
P = P(p(3,:),:);
stem((P'*U)'*x,'LineWidth',1)
xlim([1 n/4])
xlim([1 n/4])
ylabel('$\langle\mathbf{P}^Tu_k,x\rangle$','Interpreter','latex','FontSize',15)
xlabel('$k$','Interpreter','latex','FontSize',15)

figs_path = './figs\';
for p = 1:3
    ax = get(figure(p),'CurrentAxes');
    set(ax,'TickLabelInterpreter','latex')
    set(ax, 'FontSize', 20)  
    set(ax, 'YTick', []);
    print(figure(p),strcat(figs_path,'Fig_PermutedDCTSpectrum',num2str(p)),'-depsc')    
end


