function f = PassoUpwind(L,Nu,Tmax,CFL,w,f_0,order)
% f_0 deve essere un vettore che mi dà le condizoni iniziali del ciclo che
% voglo mettere in atto.

% Discretizzazione spaziale
du = 2*L / Nu;
uu = linspace(-L, L, Nu);

dt = CFL * (du /max(abs(w)));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento
f  = f_0; % Evoluzione nel tempo con schema UDS
%Passo Upwind
idx = 2:Nu-1;
if order ~=1
    idx = 3:Nu-2;
end
f_new = f;  % Array temporaneo per aggiornamento
for j = idx
    if uu(j)>-1 % questo evita di perdere massa numerica
    switch order
        case 1 % UDS
            if w>0
                g_jm12 = w * f(j-1);     % Flusso tra j-1 e j
                g_jp12 = w * f(j);   % Flusso tra j e j+1
            else
                g_jm12 = w * f(j);
                g_jp12 = w * f(j+1);
            end
        case 2 % LUDS
            if w>0
                g_jm12 = w * (1.5*f(j-1)-0.5*f(j-2));
                g_jp12 = w * (1.5*f(j)-0.5*f(j-1));
            else
                g_jm12 = w * (1.5*f(j)-0.5*f(j+1));
                g_jp12 = w * (1.5*f(j+1)-0.5*f(j));
            end
        case 3 % QUICK
            if w>0
                g_jm12 = w * (1/3*f(j)+5/6*f(j-1)-1/6*f(j-2));
                g_jp12 = w * (1/3*f(j+1)+5/6*f(j)-1/6*f(j-1));
            else
                g_jm12 = w * (1/3*f(j-1)+5/6*f(j)-1/6*f(j+1));
                g_jp12 = w * (1/3*f(j)+5/6*f(j+1)-1/6*f(j+2));
            end
    end
        
        % Aggiornamento della soluzione
        f_new(j) = f(j) - (dt/du) * (g_jp12 - g_jm12);
    end
end
% Aggiorna la soluzione
f = f_new;
