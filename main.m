clearvars; close all;clc
format short;
% numero popolazione per ogni fascia d'età
N_class = [441148; 5666380; 5962570; 6570438; 8061698; 9619516; 7964818; 6141544; 4543122];

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
Nu   = 1000;                            % Numero di punti spaziali
Tmax = 50;                            % Tempo massimo di simulazione
Tfin = 15;
CFL  = 0.9;                            % Numero di Courant-Friedrichs-Lewy
% condizione iniziale, matrice in cui ogni colonna rappresenta una fascia
% d'età
U0 = nan(max(N_class),numel(N_class));
for N_c = 1:numel(N_class)
    U0(1:N_class(N_c),N_c) = -1.5+(-1.01+1.5)*rand(N_class(N_c),1); % ho ridotto il numero di infetti iniziali (da 9.1 a 9.01)
end

%Infetti il 29 febbraio
U0(1:5,2) = 0;  % 5 infetti tra i giovani (fascia 19-50)
U0(6:20,5) = 0; % 15 infetti tra gli adulti (50-70)
U0(21:30,8) = 0; % 10 infetti tra gli anziani (70+)


% parametri
beta       = 0.25;
gamma      = 0.24;
w          = 1/14;
t_lock     = 30;
data       = cell(Tfin,1); 
hbar       = waitbar(0,'','Name','Time iterations');
flag       = 0;



for t = 1:Tfin
    if (t>t_lock)&&(~flag)
        M(1:3,1:3) = M(1:3,1:3) * 0.3;  % Giovani  %riduzioni fatte in base ad un crierio basato sul lockdown          M(4:6,4:6) = M(4:6,4:6) * 0.5;  % Adulti
        M(7:9,7:9) = M(7:9,7:9) * 0.7;  % Anziani
        beta=0.075;                     %il lockdown fa calare la possibilità di trasmissione
        flag = 1;
    end
    for N_c = 1:numel(N_class)
        waitbar(N_c/numel(N_class),hbar,sprintf('Time step: %d/%d.\n $N_c$ = %d / %d',t,Tfin,N_c,numel(N_class)));
        % PRIMO PASSO
        [f_new_tilda,U,num_bins,edges] = MonteCarlo(t,U0,beta,gamma,N_c,M(:,N_c));close all
        % SECONDO PASSO
        f_new                          = PassoUpwind(L,num_bins,Tmax,CFL,w,f_new_tilda,3);
        % AGGIORNAMENTO
        for c = 1:numel(U)
            if U(c)>-1
                U(c) = U(c)+w;
            end
        end
        U0(1:N_class(N_c),N_c) = U; 
        
        % CALCOLO S, I, R
        S=0;
        I=0;
        R=0;

           for i = 1:N_class(N_c)
                if U0(i, N_c) < -1
                    S = S + 1;
                elseif abs(U0(i, N_c)) <= 1
                    I = I + 1;
                else
                    R = R + 1;
                end
            end

        
        T           = table(t, S, I, R, 'RowNames', {'Value'},...
                       'VariableNames', {'Time Step', 'Susceptible',...
                       'Infected', 'Removed'});
        data{t,N_c} = T; % time step sulle righe, classi sulle colonne

        % plot della distribuzione della classe N_c
        % plot(edges(1:end-1),f_new)
        % legend(sprintf('Distribution plot of class %d', N_c),'Location','best')
        
    end
    % pause
end
close(hbar)
% check popolazione conservata
% check infetti            %da qui in poi ho il dubbio che ci sia un
% problema di calcolo
Tot       = zeros(numel(N_class),1);
confirmed = zeros(numel(N_class),1);
removed   = zeros(numel(N_class),1);
for t=1:numel(N_class)
    Tot(t)       = sum(data{end,t}.Susceptible+data{end,t}.Infected+data{end,t}.Removed);
    confirmed(t) = sum(data{end,t}.Infected);
    removed(t)   = sum(data{end,t}.Removed);
end
confirmed  = sum(confirmed);
removed    = sum(removed);
err        = abs(N_class-Tot)./N_class*100
err_weight = sum(err.*N_class)/sum(N_class)