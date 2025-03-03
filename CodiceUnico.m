%Prova di unificazione del codice che tiene in considerazione soltanto gli
%stati microscopici. 


clearvars; close all;clc
format short;

N_class = [441148; 5666380; 5962570; 6570438; 8061698; 9619516; 7964818; 6141544; 4543122];

% Percentuale di infetti iniziali per fascia d'età
percent_infetti = [0.12, 0.10, 0.09, 0.09, 0.10, 0.11, 0.12, 0.08, 0.07];  % 0-9, 10-19, 20-29, ..., 80+

% Numero iniziale di infetti
num_infetti_iniziali = 1000;

% Inizializzazione degli infetti nelle classi
U0 = nan(max(N_class), numel(N_class));
for i = 1:numel(N_class)
    num_infetti = round(num_infetti_iniziali * percent_infetti(i));
    U0(1:num_infetti, i) = -1+ (0+1) * rand(num_infetti, 1); 
    U0(num_infetti+1:N_class(i),i)=-3 + (-1.01 + 3) * rand(N_class(i) - num_infetti, 1); % Assegno valori per i suscettibili
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


L    = 100;
Tfin = 5;

% parametri
beta       = 0.3;
gamma      = 0.01;
w          = 1/14;
t_lock     = 9;         %data iniziale è il 29 febbraio, il lockdown il 9 marzo
data       = cell(Tfin,1); 
hbar       = waitbar(0,'','Name','Time iterations');
flag       = 0;
S=zeros(Tfin,1);
I=zeros(Tfin,1);
R=zeros(Tfin,1);

%inizio iterazioni temporali

for t = 1:Tfin
    if (t>t_lock)&&(~flag)
    M(1:3,1:3) = M(1:3,1:3) * 0.3;  % Giovani  %riduzioni fatte in base ad un crierio basato sul lockdown
    M(4:6,4:6) = M(4:6,4:6) * 0.5;  % Adulti
    M(7:9,7:9) = M(7:9,7:9) * 0.7;  % Anziani
    beta=0.075;                     %il lockdown fa calare la possibilità di trasmissione
    flag = 1;
    end
    for N_c = 1:numel(N_class)
        waitbar(N_c/numel(N_class),hbar,sprintf('Time step: %d/%d.\n $N_c$ = %d / %d',t,Tfin,N_c,numel(N_class)));
        
        
        %Montecarlo%


        a     = U0(:,N_c);
        U1    = a(~isnan(a));clear a; % per farle interagire le divido in due gruppi
        U1new = zeros(numel(U1),1);
        M_tilda=M(:,N_c);
        no=max(M_tilda);

        for j=1:size(U0,2) % chiedo un ciclo su tutte le classi
        a       = U0(:,j);
        U2      = a(~isnan(a)); clear a;
        U2new   = zeros(numel(U2),1);
        pp      = min(numel(U1),numel(U2));
        Theta_b = binornd(1,beta,pp,1);
        Theta_g = binornd(1,gamma,pp,1);
        Theta_c = binornd(1,max(M_tilda(j)/no,0.25),pp,1);

       

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
                if (U1(p)<=-1) && (abs(U2(p))<=1) % interazione S-I
                    UU1      = U2(p); % stati post interazione
                    UU2      = U2(p);
                    U1new(p) = (1-Theta_b(p)).*U1(p) + Theta_b(p).*(UU1); % aggiornamento temporale degli stati
                    U2new(p) = (1-Theta_b(p)).*U2(p) + Theta_b(p).*(UU2);
                
                elseif (abs(U1(p))<=1) && (abs(U2(p))<=1) % interazione I-I
                    UU1      = U1(p)+2; % stati post interazione
                    UU2      = U2(p);
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

     % l'interazione modifica solo la classe i

    disp(['U (stati aggiornati) min/max: ', num2str(min(U)), ' / ', num2str(max(U))]);

    %Healing Process%

    if U(N_class(N_c))>-1
        U(N_class(N_c))=U(N_class(N_c))+w;
    end

    U0(1:N_class(N_c),N_c) = U;
    
    %Ricavo S, I, R
     for i = 1:N_class(N_c)
    if U0(i, N_c) < -1
        S(t) = S(t) + 1;
    elseif abs(U0(i, N_c)) <= 1
        I(t) = I(t) + 1;
    else
        R(t) = R(t) + 1;
    end
    end

    end
end


% Grafico
figure;
plot(1:Tfin, I);
xlabel('Giorni');
ylabel('Infetti');
title('Andamento epidemico');