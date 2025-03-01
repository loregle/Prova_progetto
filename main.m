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
Tfin = 50;
CFL  = 0.9;                            % Numero di Courant-Friedrichs-Lewyc
w    = 1;
% condizione iniziale, matrice in cui ogni colonna rappresenta una fascia
% d'età
U0 = nan(max(N_class),numel(N_class));
for N_c = 1:numel(N_class)
    U0(1:N_class(N_c),N_c) = -L+(L-.99)*rand(N_class(N_c),1); % ho ridotto il numero di infetti iniziali (da 9.1 a 9.01)
end

% parametri
beta   = 1.6e-8;
gamma  = 0.24;
lambda = 0.1;
t_lock = 14;
data   = cell(Tfin,1); 
hbar   = waitbar(0,'','Name','Time iterations');
flag   = 0;
for t = 1:Tfin
    if (t>t_lock)&&(~flag)
        M(1:6,1:6) = M(1:6,1:6)./lambda;
        flag       = 1;
    end
    for N_c = 1:numel(N_class)
        waitbar(N_c/numel(N_class),hbar,sprintf('Time step: %d.\n $N_c$ = %d / %d',t,N_c,numel(N_class)));
        % PRIMO PASSO
        [f_new_tilda,U,num_bins,edges] = MonteCarlo(U0,beta,gamma,N_c,M(:,N_c));close all
        % SECONDO PASSO
        f_new = PassoUpwind(L,num_bins,Tmax,CFL,w,f_new_tilda);
        % AGGIORNAMENTO
        for c = 1:numel(U)
            if U(c)>-1
                U(c) = U(c)+w;
            end
        end
        U0(1:N_class(N_c),N_c) = U; 
        % CALCOLO S, I, R
        S = sum(f_new(1:find(edges==-1)));
        if isempty(find(edges==1, 1))
            I = sum(f_new(find(edges==-1):end));
            R = 0;
        else
            r = find(edges==1);
            if r>num_bins
                r = num_bins;
            end
            I = sum(f_new(find(edges==-1):r));
            R = sum(f_new(r:end));
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
Tot       = zeros(numel(N_class),1);
confirmed = zeros(numel(N_class),1);
for t=1:numel(N_class)
    Tot(t) = sum(data{end,t}.Susceptible+data{end,t}.Infected+data{end,t}.Removed);
    confirmed(t) = sum(data{end,t}.Infected+data{end,t}.Removed);
end
confirmed = sum(confirmed);
err = abs(N_class-Tot)./N_class*100