function [f,U,bins,edges] = MonteCarlo(U0M,beta,gamma,N_c,M) % U0M mi serve per avere la matrice di partenza
set(0,'DefaultTextInterpreter','latex')
rng(1)

a     = U0M(:,N_c);
U1    = a(~isnan(a));clear a; % per farle interagire le divido in due gruppi
U1new = zeros(numel(U1),1);

for j=1:size(U0M,2) % chiedo un ciclo su tutte le classi
    a       = U0M(:,j);
    U2      = a(~isnan(a)); clear a;
    U2new   = zeros(numel(U2),1);
    pp      = min(numel(U1),numel(U2));
    Theta_b = binornd(1,beta,pp,1);
    Theta_g = binornd(1,gamma,pp,1);
    U1      = U1(randperm(numel(U1))); % così teniamo conto del fatto che non tutte le particelle interagiscono per via del numero diverso di individui all'interno di classi diverse
    U2      = U2(randperm(numel(U2)));
    for p = 1:pp
        if (U1(p)<=-1) && (abs(U2(p))<=1) % interazione S-I
            % !!!!!!!
            % nutro dubbi riguardo l'implementazione delle
            % interazioni con la matrice di contato
            % !!!!!!!
            UU1      = ((1-beta)*U1(p) + beta*(U1(p)+2))*M(j); % stati post interazione
            UU2      = U2(p);
            U1new(p) = (1-Theta_b(p)).*U1(p) + Theta_b(p).*(UU1); % aggiornamento temporale degli stati
            U2new(p) = (1-Theta_b(p)).*U2(p) + Theta_b(p).*(UU2);
        elseif (abs(U1(p))<=1) && (abs(U2(p))<=1) % interazione I-I
            UU1      = ((1-gamma)*U1(p) + gamma*(U1(p)+2))*M(j); % stati post interazione
            UU2      = U2(p);
            U1new(p) = (1-Theta_g(p)).*U1(p) + Theta_g(p).*(UU1);
            U2new(p) = (1-Theta_g(p)).*U2(p) + Theta_g(p).*(UU2);
        else
            U1new(p) = U1(p);
            U2new(p) = U2(p);
        end
    end
end

U = U1new; % l'interazione modifica solo la classe i

h = histogram(U,'Visible','off');
% ricavo la distribuzione prendendo i valori degli istogrammi
f = h.Values;
bins = h.NumBins;
edges = h.BinEdges;
xlabel('$u$')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',1.2);
set(gca,'TickLabelInterpreter','latex')
end
