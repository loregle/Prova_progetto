
clearvars; close all;clc
format short;

N_class = [441148; 5666380; 5962570; 6570438; 8061698; 9619516; 7964818; 6141544; 4543122];

% Percentuale di infetti iniziali per fascia d'età
percent_infetti = [0.12, 0.10, 0.09, 0.09, 0.10, 0.11, 0.12, 0.08, 0.07];  % 0-9, 10-19, 20-29, ..., 80+
% Numero iniziale di infetti
num_infetti_iniziali = 888;
rim_iniziali         = 46;

L    = 100;
Tfin = 7;   % Parametri stimati per la settimana dal 28 febbraio al 5 marzo

% Inizializzazione dei suscettibili, degli infetti e dei rimossi nelle classi, rispettando i dati reali
U0 = nan(max(N_class), numel(N_class));
for i = 1:numel(N_class)
    U0(1:N_class(i),i) = -3 + (-2.01 + 3) * rand(N_class(i),1); % Impedisce valori ≥ -1 e quidi non ci sono infetti iniziali non voluti;
    U0(1:floor(num_infetti_iniziali*percent_infetti(i)),i) = -.5;
    if i == 7 || i == 8
        U0(1:23,i) = 2;
    end
    a  = U0(1:N_class(i),i);
    a  = a(randperm(N_class(i)));
    U0(1:N_class(i),i) = a;
end


M =  [19.2 4.8 3.0 7.1 3.7 3.1 2.3 1.4 1.4;
    4.8 42.4 6.4 5.4 7.5 5.0 1.8 1.7 1.7;
    3.0 6.4 20.7 9.2 7.1 6.3 2.0 0.9 0.9;
    7.1 5.4 9.2 16.9 10.1 6.8 3.4 1.5 1.5;
    3.7 7.5 7.1 10.1 13.1 7.4 2.6 2.1 2.1;
    3.1 5.0 6.3 6.8 7.4 10.4 3.5 1.8 1.8;
    2.3 1.8 2.0 3.4 2.6 3.5 7.5 3.2 3.2;
    1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2;
    1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2];


% parametri
beta       = 0.275;
gamma      = 0.01;
w          = 1/14;
hbar       = waitbar(0,'','Name','Time iterations');
S=zeros(Tfin,numel(N_class));
I=zeros(Tfin,numel(N_class));
R=zeros(Tfin,numel(N_class));

%inizio iterazioni temporali

for t = 1:Tfin
   for N_c = 1:numel(N_class)
        waitbar(N_c/numel(N_class),hbar,sprintf('Time step: %d/%d.\n $N_c$ = %d / %d',t,Tfin,N_c,numel(N_class)));


        %Montecarlo%


        a     = U0(:,N_c);
        U1    = a(~isnan(a));clear a; % per farle interagire le divido in due gruppi
        U1new = -2*ones(numel(U1),1); %-2 è stato scelto per non creare infetti artificisli durante le interazioni.
        M_tilda=M(:,N_c);
        no=norm(M_tilda);

        for j=1:size(U0,2) % chiedo un ciclo su tutte le classi
            a       = U0(:,j);
            U2      = a(~isnan(a)); clear a;
            U2new   = -2*ones(numel(U2),1);  
            pp      = min(numel(U1),numel(U2));
            Theta_b = binornd(1,beta,pp,1);
            Theta_g = binornd(1,gamma,pp,1);
            Theta_c = binornd(1,M_tilda(j)/no,pp,1);


            %debug
            disp(['N_c = ', num2str(N_c), ' - Popolazione coinvolta: ', num2str(pp)]);
            disp(['Theta_c attivi: ', num2str(sum(Theta_c)), '/', num2str(pp)]);
            disp(['Theta_b attivi: ', num2str(sum(Theta_b)), '/', num2str(pp)]);
            disp(['Theta_g attivi: ', num2str(sum(Theta_g)), '/', num2str(pp)]);


            for p = 1:pp
                if Theta_c(p)==0
                    U1new(p)=U1(p);
                    U2new(p)=U2(p);
                else
                    if (U1(p)<-1) && (abs(U2(p))<=1) % interazione S-I
                        UU1      = U2(p); % stati post interazione
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
            U1(1:pp)     = U1new(1:pp);
            U1(pp+1:end) = U1(pp+1:end);
        end

        % l'interazione modifica solo la classe i

        U = U1;

        disp(['U (stati aggiornati) min/max: ', num2str(min(U)), ' / ', num2str(max(U))]);

        %Healing Process%

        if abs(U(N_class(N_c)))<1
            U(N_class(N_c)) = U(N_class(N_c))+w; %la moltiplicazione per dt è implicita in quanto dt=1 e la velocità è stimata in giorni
        end

        U0(1:N_class(N_c),N_c) = U;

        % Ricavo S, I, R
        for i = 1:N_class(N_c)
            if U0(i, N_c) < -1
                S(t,N_c) = S(t,N_c) + 1;
            elseif abs(U0(i, N_c)) < 1
                I(t,N_c) = I(t,N_c) + 1;
            else
                R(t,N_c) = R(t,N_c) + 1;
            end
        end

    end
end



Tot       = sum(S,2)+sum(I,2)+sum(R,2);
confirmed = sum(I,2);
removed   = sum(R,2);


data_Inf=[888, 1128, 1694, 2036, 2502,3089,3858];
data_Rem = [46, 46, 83, 149, 160, 276, 414];

% Grafico Infetti
figure(1)

plot(1:Tfin, confirmed, 1:Tfin, data_Inf, 'r');
xlabel('Days');
ylabel('Confirmed');
title('Epidemic Progression: Infected Individuals');
legend("Simulated data", "Test data")

% Grafico Removed
figure(2)
plot(1:Tfin, removed, 1:Tfin, data_Rem, 'r')
xlabel('Days')
ylabel('Removed')
title('Epidemic Progression: Removed Individuals')
legend("Simulated data", "Test data")

% Grafici Infetti classi
figure(3)

plot(1:Tfin, I(:,1))
hold on
plot(1:Tfin, I(:,2))
plot(1:Tfin, I(:,3))
plot(1:Tfin, I(:,4))
plot(1:Tfin, I(:,5))
plot(1:Tfin, I(:,6))
plot(1:Tfin, I(:,7))
plot(1:Tfin, I(:,8))
plot(1:Tfin, I(:,9))
legend('Class 1','Class 2','Class 3','Class 4','Class 5','Class 6','Class 7','Class 8','Class 9',s'Location','best')
xlabel('Days')
ylabel('Confirmed')
title('Epidemic Progression: Infected Individuals by classes')


