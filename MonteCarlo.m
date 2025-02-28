function [f,U,n,edges] = MonteCarlo(U0,beta,gamma,N,M,U0M) %U0M mi serve per avere la matrice di partenza
set(0,'DefaultTextInterpreter','latex')
rng(1)

% eps    = 1e-2; %Parametro di riscalamento per il regime quasi-invariante
% dt     = eps; %Passo di discretizzazione temporale
% Tfin   = 1e1; %Tempo finale
% nmax   = floor(Tfin/dt); %Numero di iterazioni temporali

% hbar = waitbar(0,'','Name','Iterazioni');
for n=1:1
    % waitbar(n/nmax,hbar,sprintf('$n$ = %d / %d',n,nmax));
    % U     = U(randperm(N));
    U1    = U0; % per farle interagire le divido in due gruppi
    U1new = zeros(length(U1));
    Theta_b = binornd(1,beta,N/2,1); % per quando N/2 non è intero
    Theta_g = binornd(1,gamma,N/2,1);
    for j=1:size(U0M,2) %chiedo un ciclo su tutte le classi
        U2= U0M(:,j);
        U2new= zeros(length(U2));
        for p = 1:min(numel(U1),numel(U2))
            if (U1(p)<=-1) && (abs(U2(p))<=1) % interazione S-I
                    % !!!!!!!
                    % nutro dubbi riguardo l'implementazione delle
                    % interazioni con la matrice di contato
                    % !!!!!!!
                    UU1 = ((1-beta)*U1(p) + beta*(U1(p)+2))*M(j); % stati post interazione
                    UU2 = U2(p);
                    U1new(p) = (1-Theta_b(p)).*U1(p) + Theta_b(p).*(UU1); % aggiornamento temporale degli stati
                    U2new(p) = (1-Theta_b(p)).*U2(p) + Theta_b(p).*(UU2);
                elseif (abs(U1(p))<=1) && (abs(U2(p))<=1) % interazione I-I
                    UU1 = ((1-gamma)*U1(p) + gamma*(U1(p)+2))*M(j); % stati post interazione
                    UU2 = U2(p);
                    U1new(p) = (1-Theta_g(p)).*U1(p) + Theta_g(p).*(UU1);
                    U2new(p) = (1-Theta_g(p)).*U2(p) + Theta_g(p).*(UU2);
                else
                    U1new(p) = U1(p);
                    U2new(p) = U2(p);
            end
         end
    end
   
    U = U1new;
end
% close(hbar)

h = histogram(U,'Visible','off');
% ricavo la distribuzione prendendo i valori degli istogrammi
f = h.Values;
n = h.NumBins;
edges = h.BinEdges;
% figure(1)
% plot(h.BinLimits(1):h.BinWidth:h.BinLimits(2)-h.BinWidth,h.Values)
% legend(sprintf('Distribution of class i at time step %d',t))
xlabel('$u$')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',1.2);
set(gca,'TickLabelInterpreter','latex')
end