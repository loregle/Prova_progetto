function f = PassoUpwind(L,Nu,Tmax,CFL,w,f_0)
% f_0 deve essere un vettore che mi dà le condizoni iniziali del ciclo che
% voglo mettere in atto.

% Discretizzazione spaziale
du = 2*L / Nu;
uu = linspace(-L, L, Nu);

dt = CFL * (du /max(abs(w(uu))));
Nt = ceil(Tmax / dt);  % Numero di passi temporali
dt = Tmax / Nt;        % Ricalcolo per adattamento

f=f_0;% Evoluzione nel tempo con schema UDS
%Passo Upwind
f_new = f;  % Array temporaneo per aggiornamento
for j = 2:Nu-1
    if uu(j)>-1 % questo evita di perdere massa numerica
        % Flusso numerico per w > 0
        g_jm12 = w(uu(j)) * f(j-1);     % Flusso tra j-1 e j
        g_jp12 = w(uu(j)) * f(j);   % Flusso tra j e j+1
        % Aggiornamento della soluzione
        f_new(j) = f(j) - (dt/du) * (g_jp12 - g_jm12);
    end
end
% Aggiorna la soluzione
f = f_new;
