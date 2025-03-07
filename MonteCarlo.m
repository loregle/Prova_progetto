function [f,U,bins,edges] = MonteCarlo(n,U0M,beta,gamma,N_c,M) % U0M mi serve per avere la matrice di partenza
% set(0,'DefaultTextInterpreter','latex')
rng(1)

a     = U0M(:,N_c);
U1    = a(~isnan(a));clear a; % per farle interagire le divido in due gruppi
U1new = zeros(numel(U1),1);
no=max(M);

for j=1:size(U0M,2) % chiedo un ciclo su tutte le classi
    a       = U0M(:,j);
    U2      = a(~isnan(a)); clear a;
    U2new   = zeros(numel(U2),1);
    pp      = min(numel(U1),numel(U2));
    Theta_b = binornd(1,beta,pp,1);
    Theta_g = binornd(1,gamma,pp,1);
    Theta_c=binornd(1,1,pp,1);

    if N_c==7
            Theta_c(1)=1; %altrimenti l'epidemia non parte
            Theta_c(2)=1;
    end

    %debug
    disp(['N_c = ', num2str(N_c), ' - Popolazione coinvolta: ', num2str(pp)]);
    disp(['Theta_c attivi: ', num2str(sum(Theta_c)), '/', num2str(pp)]);
    disp(['Theta_b attivi: ', num2str(sum(Theta_b)), '/', num2str(pp)]);
    disp(['Theta_g attivi: ', num2str(sum(Theta_g)), '/', num2str(pp)]);


    for p = 1:pp
        
        Theta_b(p)=Theta_b(p)*(U1(p)<-1)*(-1<U2(p))*(U2(p)<1); %mi assicuro che le interazioni ci siano tra le classi giuste
        Theta_g(p)=Theta_g(p)*(-1<U1(p))*(U1(p)<1)*(-1<U2(p))*(U2(p)<1);
        

        if Theta_c(p)==0
                U1new(p)=U1(p);
                U2new(p)=U2(p);
        else
            if (U1(p)<=-1) && (abs(U2(p))<=1) % interazione S-I
                UU1      = (U1(p)+2); % stati post interazione
                UU2      = U2(p);
                U1new(p) = (1-Theta_b(p)).*U1(p) + Theta_b(p).*(UU1); % aggiornamento temporale degli stati
                U2new(p) = (1-Theta_b(p)).*U2(p) + Theta_b(p).*(UU2);
                
            elseif (abs(U1(p))<=1) && (abs(U2(p))<=1) % interazione I-I
                UU1      = (U1(p)+2); % stati post interazione
                UU2      = U2(p);
                U1new(p) = (1-Theta_g(p)).*U1(p) + Theta_g(p).*(UU1);
                U2new(p) = (1-Theta_g(p)).*U2(p) + Theta_g(p).*(UU2);
            else
                U1new(p) = U1(p);
                U2new(p) = U2(p);
            end
        end
    end
end

U = U1new; % l'interazione modifica solo la classe i
disp(['U (stati aggiornati) min/max: ', num2str(min(U)), ' / ', num2str(max(U))]);

h = histogram(U,'Visible','off');
% ricavo la distribuzione prendendo i valori degli istogrammi
f     = h.Values;
bins  = h.NumBins;
edges = h.BinEdges;
% xlabel('$u$')
% set(0,'DefaultAxesFontSize',18)
% set(0,'DefaultLineLineWidth',1.2);
% set(gca,'TickLabelInterpreter','latex')
end
