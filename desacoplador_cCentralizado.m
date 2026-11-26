clear; clc; close all
s = tf('s');

%% 1. PROCESSO
G11 = 3.04/s;
G12 = -278.28/(s*(30*s+1)*(6*s+1));
G21 = 0.052/s;
G22 = 319.47/(s*(30*s+1)*(6*s+1));
G = [G11 G12; G21 G22];

% ganhos dos integradores: G = K/s * dinamica
K11 = 3.04;    K12 = -278.28;
K21 = 0.052;   K22 = 319.47;
K = [K11 K12; K21 K22];

fprintf('Processo analisado: \n'); G

%% 2. PARAMETROS DE SINTONIA E SIMULACAO
lambda1  = 5;    lambda2  = 0.5;     % malhas principais
lambda12 = 20;   lambda21 = 20;      % acoplamentos (so no centralizado)
alpha    = 0.1;                      % filtro da derivada do PD
Tf_min   = 0.001;                    % filtro que torna D21 proprio

% centralizado de Morilla (Caso 3, processos SISO equivalentes)
ordem_mor = 2;                       % 1 = PI | 2 = PID (necessario pela dinamica (30s+1)(6s+1))
z_mor     = [0.1 0.02];              % zero z_i de lbar_i = (s + z_i)/s (z_2 = 0.1 instabiliza)
delta_mor = [1 1];                   % amortecimento desejado (1 = critico, wn = 2*z_i)
N_mor     = 10;                      % filtro da derivada: Tf = |Kd/(Kp*N)|
w_mor     = logspace(-3, 0, 300).';  % faixa de frequencia do ajuste de m_ji
t_d       = 10;                      % instante da perturbacao de carga (teste separado)
d_amp     = [1e-3 0];                % degrau de perturbacao somado a u1 e u2

% caracterizacao da planta com desacoplador ideal (desacoplador e centralizado)
K1 = 3.0853;   tau1 = 0;    t0_1 = 0;
K2 = 324.23;   tau2 = 33;   t0_2 = 3;

amp_degrau = 1;
t_r1 = 5;
t_r2 = 150;
Tsim = 300;
Ts   = 0.05;

t  = (0:Ts:Tsim)';
r1 = amp_degrau * (t >= t_r1);
r2 = amp_degrau * (t >= t_r2);
r  = [r1, r2];

%% 3. DESACOPLADOR IDEAL
% u' = Dm*u, Dm = [1 -D12; -D21 1] -> G*Dm diagonal
% D12 e D21 sem sinal porque sao os numD12/numD21 do Simulink
D12 = minreal(G12/G11);
D21 = minreal(G21/G22 * tf(1, [Tf_min 1])^2);   % filtro: G21/G22 e improprio
Dm  = [1 -D12; -D21 1];

fprintf('Desacoplador ideal:\n'); Dm
H = minreal(G*Dm);
fprintf('G*Dm (deve ser diagonal):\n'); H

fprintf('Parametros caracterizados\n');
fprintf('Malha 1: K1=%.4f | tau1=%g | t0_1=%g\n', K1, tau1, t0_1);
fprintf('Malha 2: K2=%.4f | tau2=%g | t0_2=%g\n\n', K2, tau2, t0_2);

% sintonia lambda: PD com filtro, Kc = 1/(K*(lambda + t0)), tD = tau
Kc1_prime = 1/(K1*(lambda1 + t0_1));
Kc2_prime = 1/(K2*(lambda2 + t0_2));
C_d = [pd_filtrado(Kc1_prime, tau1, alpha) 0; 0 pd_filtrado(Kc2_prime, tau2, alpha)];

fprintf('Kc1_prime = %.6f | Kc2_prime = %.6f\n', Kc1_prime, Kc2_prime);

[y_com, c_com, u_com, L_com, S_com, T_com] = malha_fechada(G, Dm, C_d, r, t);
fprintf('Polos (desacoplador ideal): '); fprintf('%.4f  ', pole(T_com)); fprintf('\n');

%% 4. DESACOPLADOR INVERTIDO
% eqs. (2.19)-(2.20): u'1 = u1 + D12*u'2 ;  u'2 = u2 + D21*u'1  ->  u' = (I - Doff)^-1 u
% D12 e D21 calculados com o modelo com erro G_hat

usar_erro_ident = true;                 % false -> modelo perfeito (desacoplamento total)
erroK = 2*[0.50 1.05; 2.08 0.13];       % erro relativo em cada Gij

if ~usar_erro_ident
    erroK = zeros(2);
end
G_hat = [G11*(1+erroK(1,1))  G12*(1+erroK(1,2));
         G21*(1+erroK(2,1))  G22*(1+erroK(2,2))];

D12_inv  = minreal(-G_hat(1,2)/G_hat(1,1));
D21_inv  = minreal(-G_hat(2,1)/G_hat(2,2) * tf(1, [Tf_min 1])^2);
Doff_inv = [0 D12_inv; D21_inv 0];
Pre_inv  = feedback(tf(eye(2)), Doff_inv, +1);   % u' = Pre_inv*u

F_hat = minreal(G_hat * Pre_inv);   % planta desacoplada pelo modelo (diagonal)
F_inv = minreal(G * Pre_inv);       % planta desacoplada real (fora da diagonal = acoplamento residual)

fprintf('\n--- Desacoplador invertido ---\n');
fprintf('D12_inv(s) = -G^12/G^11:\n'); D12_inv
fprintf('D21_inv(s) = -G^21/G^22 (com filtro):\n'); D21_inv
fprintf('F^(s) do modelo (deve ser diag(G^11, G^22)):\n'); zpk(F_hat)
fprintf('F(s) real:\n'); zpk(F_inv)

%% 4.1 Ensaio em malha aberta (planta real com e sem desacoplador)
Tma  = 600;
t_u1 = 10;   t_u2 = 300;   du = 1;

ta = (0:Ts:Tma)';
Na = length(ta);
ua = du*[ta >= t_u1, ta >= t_u2];

ya  = lsim(F_inv, ua, ta);     % com desacoplador
yas = lsim(G,     ua, ta);     % sem desacoplador

%% 4.2 Identificacao: Fit3 (Cardenas, 2024) x minimos quadrados
%   ya  -> planta + desacoplador invertido (F = G*Pre_inv) -> sintonia do PD (4.3)
%   yas -> planta sem desacoplador (G)                      -> centralizado identificado (5.1)
ruido = 0;                                  % desvio padrao do ruido de medicao (0 = sem ruido)
rng(0);
ym  = ya  + ruido*randn(Na,2);
yms = yas + ruido*randn(Na,2);

Kint = @(sys) real(1j*1e-6*squeeze(freqresp(sys, 1e-6)));   % lim s->0 de s*sys(s)
KF_real = Kint(F_inv);

% janelas das inclinacoes [antes do degrau; depois do degrau]
% comecam 200 s apos o degrau (transitorio de 30 s + 6 s ja passou)
jan_u1 = [1, t_u1 - 1;            t_u1 + 200, t_u2 - 5];   % degrau em u1
jan_u2 = [t_u1 + 200, t_u2 - 5;   t_u2 + 200, Tma - 1];    % degrau em u2

% --- Fit3: K* = diferenca das inclinacoes, tau e t0 por 28,3% e 63,2% da derivada
[Kp1_f3, tau1_f3, t01_f3] = fit3(ta, ym(:,1), t_u1, jan_u1, du, Ts);
[Kp2_f3, tau2_f3, t02_f3] = fit3(ta, ym(:,2), t_u2, jan_u2, du, Ts);

K_f3 = zeros(2);   tau_f3 = zeros(2);   t0_f3 = zeros(2);
jan = {jan_u1, jan_u2};   t_u = [t_u1 t_u2];
for i = 1:2
    for j = 1:2
        [K_f3(i,j), tau_f3(i,j), t0_f3(i,j)] = fit3(ta, yms(:,i), t_u(j), jan{j}, du, Ts);
    end
end

% --- minimos quadrados: y = K*integral(u) - K*theta*u
% integrador com atraso theta: y = K*integral(u(t - theta)) ~ K*integral(u) - K*theta*u
% vale depois do transitorio -> fora do ajuste os mesmos 200 s apos cada degrau do Fit3
PHI = [Ts*cumtrapz(ua), ua];
usa = ~(ta >= t_u1 & ta < t_u1 + 200) & ~(ta >= t_u2 & ta < t_u2 + 200);

th_a  = PHI(usa,:) \ ym(usa,:);             % [K; -K*theta] de cada saida
H_mq  = th_a(1:2,:).';
th_mq = -th_a(3:4,:).' ./ H_mq;             % atraso equivalente de cada elemento

th_s  = PHI(usa,:) \ yms(usa,:);
K_mq  = th_s(1:2,:).';
t0_mq = -th_s(3:4,:).' ./ K_mq;

Kp1_mq = H_mq(1,1);   t01_mq = th_mq(1,1);
Kp2_mq = H_mq(2,2);   t02_mq = th_mq(2,2);

cab   = @() fprintf('%-6s %12s %12s %12s %14s %12s\n', '', 'Real', 'Fit3', 'MQ', 'erro Fit3 %', 'erro MQ %');
linha = @(nome, re, f3, mq) fprintf('%-6s %12.4f %12.4f %12.4f %14.4f %12.4f\n', ...
    nome, re, f3, mq, 100*abs((f3 - re)/re), 100*abs((mq - re)/re));

fprintf('\n--- Identificacao: Fit3 x minimos quadrados (ruido = %g) ---\n', ruido);
fprintf('Planta + desacoplador invertido (F = G*Pre_inv) -> usada na sintonia do PD\n');
cab();
linha('F11', KF_real(1,1), Kp1_f3, Kp1_mq);
linha('F22', KF_real(2,2), Kp2_f3, Kp2_mq);
fprintf('Fit3: malha 1 tau = %.4f  t0 = %.4f | malha 2 tau = %.4f  t0 = %.4f\n', tau1_f3, t01_f3, tau2_f3, t02_f3);
fprintf('MQ  : malha 1 atraso eq. = %.4f      | malha 2 atraso eq. = %.4f\n', t01_mq, t02_mq);
fprintf('      (comparar com tau + t0 do Fit3: %.4f | %.4f)\n', tau1_f3 + t01_f3, tau2_f3 + t02_f3);

fprintf('\nPlanta sem desacoplador (G) -> usada no centralizado identificado (5.1)\n');
cab();
nomesK = {'K11', 'K12'; 'K21', 'K22'};
for i = 1:2
    for j = 1:2
        linha(nomesK{i,j}, K(i,j), K_f3(i,j), K_mq(i,j));
    end
end
fprintf('Fit3 tau:\n');                  disp(tau_f3);
fprintf('Fit3 t0:\n');                   disp(t0_f3);
fprintf('MQ atraso equivalente:\n');     disp(t0_mq);
fprintf('(real: G11 e G21 = 0 | G12 e G22 = 30 + 6 = 36 s de constante de tempo)\n');

%% 4.3 Sintonia PD, eq. (15): Kc = 1/(K*(lambda + t0)), tD = tau
ident_sintonia = 'MQ';     % 'MQ' ou 'Fit3'

if strcmp(ident_sintonia, 'MQ')
    % modelo integrador + atraso (sem tau) -> o PD vira P
    Ks1 = Kp1_mq;   tD1_inv = 0;         t01 = t01_mq;
    Ks2 = Kp2_mq;   tD2_inv = 0;         t02 = t02_mq;
else
    Ks1 = Kp1_f3;   tD1_inv = tau1_f3;   t01 = t01_f3;
    Ks2 = Kp2_f3;   tD2_inv = tau2_f3;   t02 = t02_f3;
end
fprintf('\nSintonia com %s: K1 = %.4f  t0 = %.4f | K2 = %.4f  t0 = %.4f\n', ident_sintonia, Ks1, t01, Ks2, t02);

Kc1_inv = 1/(Ks1*(lambda1 + t01));
Kc2_inv = 1/(Ks2*(lambda2 + t02));
C_inv = [pd_filtrado(Kc1_inv, tD1_inv, alpha) 0; 0 pd_filtrado(Kc2_inv, tD2_inv, alpha)];

fprintf('Kc1_inv = %.6f   tD1 = %.4f\n', Kc1_inv, tD1_inv);
fprintf('Kc2_inv = %.6f   tD2 = %.4f\n', Kc2_inv, tD2_inv);

%% 4.4 Malha fechada
[y_inv, c_inv, v_inv, L_inv, S_inv, T_inv] = malha_fechada(G, Pre_inv, C_inv, r, t);
fprintf('Polos (desacoplador invertido): '); fprintf('%.4f  ', pole(T_inv)); fprintf('\n');

%% 5. CONTROLADOR CENTRALIZADO
% Gc = diag(PD_1, PD_2) * M,  M = K^-1 * N,  N_ij = 1/(lambda_ij + t0_i)
detK = K11*K22 - K12*K21;
fprintf('\ndet(K) = %.4f | cond(K) = %.2f\n', detK, cond(K));

Lam  = [lambda1 lambda12; lambda21 lambda2];
t0_i = [t0_1; t0_2];
N    = 1 ./ (Lam + t0_i);

fprintf('Matriz N = [1/(lambda_ij + t0_i)]:\n'); disp(N);
fprintf('Autovalores de N (devem ter parte real > 0): ');
fprintf('%.4f  ', eig(N)); fprintf('\n');
if any(real(eig(N)) <= 0)
    warning('Lambdas escolhidos tornam a malha instável (autovalor de N <= 0).');
end

M = K \ N;
fprintf('Ganho do controlador centralizado (M = K^-1 N):\n'); disp(M);

tau_j = [tau1; tau2];
Gpd   = [pd_filtrado(1, tau_j(1), alpha) 0; 0 pd_filtrado(1, tau_j(2), alpha)];
Gc_centralizado = Gpd * M;
fprintf('Controlador centralizado Gc(s) = diag(PD_1, PD_2) * M:\n'); Gc_centralizado

[yC, c_cent, u_cent, L_cent, S_cent, T_cent] = malha_fechada(G, tf(eye(2)), Gc_centralizado, r, t);
fprintf('Polos (centralizado): '); fprintf('%.4f  ', pole(T_cent)); fprintf('\n');

%% 5.1 Centralizado projetado com o modelo identificado: MQ x Fit3
% mesmo metodo da secao 5, trocando K, t0_i e tau_i pelos identificados (diagonal de G)
% MQ nao estima tau -> PD vira P | os dois sao simulados na planta real G
t0_i_mq  = diag(t0_mq);
M_mq     = K_mq \ (1 ./ (Lam + t0_i_mq));
Gc_mq    = tf(M_mq);

t0_i_f3  = diag(t0_f3);   tau_i_f3 = diag(tau_f3);
M_f3     = K_f3 \ (1 ./ (Lam + t0_i_f3));
Gc_f3    = [pd_filtrado(1, tau_i_f3(1), alpha) 0; 0 pd_filtrado(1, tau_i_f3(2), alpha)] * M_f3;

fprintf('\n--- Centralizado com modelo identificado ---\n');
fprintf('M (caracterizacao original):\n');   disp(M);
fprintf('M (MQ):\n');                        disp(M_mq);
fprintf('M (Fit3):\n');                      disp(M_f3);
fprintf('tau_i Fit3 = [%.4f %.4f] | t0_i Fit3 = [%.4f %.4f] | t0_i MQ = [%.4f %.4f]\n', ...
    tau_i_f3, t0_i_f3, t0_i_mq);

[y_cmq, c_cmq, ~, ~, ~, T_cmq] = malha_fechada(G, tf(eye(2)), Gc_mq, r, t);
[y_cf3, c_cf3, ~, ~, ~, T_cf3] = malha_fechada(G, tf(eye(2)), Gc_f3, r, t);
fprintf('Polos (centralizado MQ):   '); fprintf('%.4f  ', pole(T_cmq)); fprintf('\n');
fprintf('Polos (centralizado Fit3): '); fprintf('%.4f  ', pole(T_cf3)); fprintf('\n');

figure('Color','w', 'Name', 'Centralizado: modelo identificado por MQ x Fit3', 'NumberTitle', 'off');
sgtitle('Centralizado projetado com MQ x Fit3 (simulado na planta real)');
for m = 1:2
    subplot(2,2,m);
    plot(t, r(:,m), 'k-', 'LineWidth', 1.2); hold on;
    plot(t, yC(:,m), 'k:', 'LineWidth', 1.2);
    plot(t, y_cmq(:,m), 'b-', 'LineWidth', 1.4);
    plot(t, y_cf3(:,m), 'r--', 'LineWidth', 1.4);
    grid on; box on;
    ylabel(sprintf('y_%d', m)); title(sprintf('Saída y_%d', m));
    legend(sprintf('setpoint r_%d', m), 'Caracterização original', 'MQ', 'Fit3', 'Location', 'best');

    subplot(2,2,m+2);
    plot(t, c_cent(:,m), 'k:', 'LineWidth', 1.2); hold on;
    plot(t, c_cmq(:,m), 'b-', 'LineWidth', 1.4);
    plot(t, c_cf3(:,m), 'r--', 'LineWidth', 1.4);
    grid on; box on;
    xlabel('Tempo [s]'); ylabel(sprintf('u_%d', m)); title(sprintf('Sinal de controle u_%d', m));
    legend('Caracterização original', 'MQ', 'Fit3', 'Location', 'best');
end

%% 6. CONTROLADOR CENTRALIZADO DE MORILLA (CASO 3, PROCESSOS SISO EQUIVALENTES)
% malha desejada: l_i = k_i*lbar_i/s, lbar_i = (s + z_i)/s
%   -> h_i = k_i(s + z_i)/(s^2 + k_i s + k_i z_i)  (eq. 29)
% k_i = 4*z_i*delta_i^2 (eq. 30)
% processo SISO equivalente: g~_ij = +-det(G)/g_(3-i,3-j)
% m_ji = g~_ij*exp(tau~_ij s)/lbar_i, ajustado por b0/(a2 s^2 + s + a0) (eq. 34)
% controlador: k_ji(s) = k_i*D(j,i)(s)*exp(-tau(k_ji) s),  D = Kp + Ki/s + Kd s/(Tf s + 1)

% atrasos, eqs. (12)-(16) e (32)
Tau_tot = totaldelay(G);
tau_n   = min(Tau_tot(1,1) + Tau_tot(2,2), Tau_tot(1,2) + Tau_tot(2,1));
tau_til = zeros(2);
for i = 1:2
    for j = 1:2
        tau_til(i,j) = tau_n - Tau_tot(3-i, 3-j);
    end
end
tau_mor = max(tau_til, [], 2);

lbar  = {(s + z_mor(1))/s, (s + z_mor(2))/s};
k_mor = 4*z_mor.*delta_mor.^2;

jw    = 1j*w_mor;
Gw    = freqresp(G, w_mor);
gw    = @(a, b) squeeze(Gw(a, b, :));
detGw = gw(1,1).*gw(2,2) - gw(1,2).*gw(2,1);

Kp_mor = zeros(2);   Ki_mor = zeros(2);   Kd_mor = zeros(2);   tau_pid = zeros(2);
m_real = cell(2);    m_aprox = cell(2);   D_mor = cell(2);
for i = 1:2
    lbar_w = squeeze(freqresp(lbar{i}, w_mor));
    for j = 1:2
        sinal = 1 - 2*(i ~= j);
        gtil  = sinal*detGw ./ gw(3-i, 3-j);
        m     = gtil .* exp(jw*tau_til(i,j)) ./ lbar_w;

        % minimos quadrados: m = b0/(a2 s^2 + s + a0)  ->  b0 - a0*m - a2*s^2*m = s*m
        if ordem_mor == 1
            A = [ones(size(m)), -m];
        else
            A = [ones(size(m)), -m, -m.*jw.^2];
        end
        x  = [real(A); imag(A)] \ [real(m.*jw); imag(m.*jw)];
        b0 = x(1);   a0 = x(2);   a2 = 0;
        if ordem_mor == 2, a2 = x(3); end

        Kp_mor(j,i) = 1/b0;
        Ki_mor(j,i) = a0/b0;
        Kd_mor(j,i) = a2/b0;
        m_real{j,i}  = m;
        m_aprox{j,i} = b0./(a2*jw.^2 + jw + a0);

        Dji = Kp_mor(j,i) + Ki_mor(j,i)/s;
        if abs(Kd_mor(j,i)) > 1e-6*abs(Kp_mor(j,i))
            Tf_d = abs(Kd_mor(j,i)/(Kp_mor(j,i)*N_mor));
            Dji  = Dji + Kd_mor(j,i)*s/(Tf_d*s + 1);
        else
            Kd_mor(j,i) = 0;                          % derivada desprezivel -> PI
        end
        D_mor{j,i}   = minreal(tf(Dji), 1e-9);
        tau_pid(j,i) = tau_mor(i) - tau_til(i,j);
    end
end

K_mor = [D_mor{1,1}*k_mor(1)*exp(-tau_pid(1,1)*s), D_mor{1,2}*k_mor(2)*exp(-tau_pid(1,2)*s);
         D_mor{2,1}*k_mor(1)*exp(-tau_pid(2,1)*s), D_mor{2,2}*k_mor(2)*exp(-tau_pid(2,2)*s)];

fprintf('\n--- Centralizado de Morilla (Caso 3) ---\n');
fprintf('k1 = %.4f   k2 = %.4f\n', k_mor);
for j = 1:2
    for i = 1:2
        fprintf('D%d%d: Kp = %.6f   Ki = %.6f   Kd = %.6f\n', j, i, Kp_mor(j,i), Ki_mor(j,i), Kd_mor(j,i));
    end
end

[y_mor, c_mor, u_mor, L_mor, S_mor, T_mor] = malha_fechada(G, tf(eye(2)), K_mor, r, t);
fprintf('Polos (Morilla): '); fprintf('%.4f  ', pole(T_mor)); fprintf('\n');

%% 7. CONTROLADOR DESCENTRALIZADO
% um PD por malha (y1-u1, y2-u2), sem desacoplador, sintonizado so com a diagonal Gii
% RGA(1,1) > 0 confirma o pareamento; perto de 1 -> pouca interacao
RGA = K .* inv(K).';
fprintf('\n--- Controlador descentralizado ---\n');
fprintf('RGA:\n'); disp(RGA);

Kc1_dec = 1/(K11*(lambda1 + t0_1));   tD1_dec = tau1;
Kc2_dec = 1/(K22*(lambda2 + t0_2));   tD2_dec = tau2;
fprintf('Kc1_dec = %.6f   tD1 = %g | Kc2_dec = %.6f   tD2 = %g\n', Kc1_dec, tD1_dec, Kc2_dec, tD2_dec);

C_dec = [pd_filtrado(Kc1_dec, tD1_dec, alpha) 0; 0 pd_filtrado(Kc2_dec, tD2_dec, alpha)];
[y_dec, c_dec, u_dec, L_dec, S_dec, T_dec] = malha_fechada(G, tf(eye(2)), C_dec, r, t);
fprintf('Polos (descentralizado): '); fprintf('%.4f  ', pole(T_dec)); fprintf('\n');

%% 8. COMPARACAO FINAL
nomes = {'Desacoplador', 'Desac. invertido', 'Centralizado', 'Morilla', 'Descentralizado'};
ne    = numel(nomes);

Y  = {y_com, y_inv, yC, y_mor, y_dec};
Cc = {c_com, c_inv, c_cent, c_mor, c_dec};     % u  = saida do controlador
U  = {u_com, v_inv, u_cent, u_mor, u_dec};     % u' = entrada da planta (sem desacoplador: u = u')
S_mimo = {S_com, S_inv, S_cent, S_mor, S_dec};
T_mimo = {T_com, T_inv, T_cent, T_mor, T_dec};

%% 8.1 Resposta do controlador
figure('Color','w', 'Name', 'Resposta do controlador', 'NumberTitle', 'off');
sgtitle('Resposta do controlador (u = saída do controlador, u'' = entrada da planta)');
for j = 1:ne
    subplot(ne,1,j);
    plot(t, Cc{j}(:,1), 'k:', t, Cc{j}(:,2), 'k-.', t, U{j}(:,1), t, U{j}(:,2), 'LineWidth', 1.2);
    grid on; ylabel('controle'); title(nomes{j});
    legend('u_1', 'u_2', 'u''_1', 'u''_2', 'Location', 'best');
end
xlabel('Tempo [s]');

%% 8.2 Desempenho no tempo
% IAE  : IAE de seguimento
% SII  : IAE da interacao
% Dmax : desvio maximo por interacao (mesmas janelas do SII)
% Mp   : sobressinal [%]
% ts   : tempo de acomodacao (2%)
ind = [];
for j = 1:ne
    ind = [ind, indices(t, Y{j}, r, t_r1, t_r2, amp_degrau)]; %#ok<AGROW>
end

campos_t  = {'IAE', 'SII', 'Dmax', 'Mp', 'ts'};
rotulos_t = {'IAE seguimento', 'SII', 'Desvio max interacao', 'Mp [%]', 'ts 2% [s]'};
imprime_tabela('Desempenho no tempo', nomes, ind, campos_t, rotulos_t, 'malha');

%% 8.3 Esforco de controle (entrada da planta u')
% pico |u'| ; ISU = integral u'^2
ef = [];
for j = 1:ne
    ef = [ef, esforco(t, U{j})]; %#ok<AGROW>
end

campos_u  = {'pico', 'ISU'};
rotulos_u = {'Pico |u''|', 'ISU'};
imprime_tabela('Esforco de controle', nomes, ef, campos_u, rotulos_u, 'entrada');

%% 8.4 Processo como um todo (MIMO)
% gamma = sigma_max/sigma_min dos ganhos dos integradores (baixa frequencia)
fprintf('\n--- Processo como um todo ---\n');
fprintf('RGA(1,1) do processo                   = %.4f\n', K11*K22/detK);
fprintf('gamma(K) processo                      = %.4f\n', cond(K));
fprintf('gamma planta vista pelo controlador:\n');
fprintf('   Desacoplador (G*Dm)                 = %.4f\n', cond(Kint(G*Dm)));
fprintf('   Desac. invertido (G*Pre_inv, real)  = %.4f\n', cond(KF_real));
fprintf('   Centralizado (G)                    = %.4f\n', cond(K));
fprintf('   Morilla (G)                         = %.4f\n', cond(K));
fprintf('   Descentralizado (G)                 = %.4f\n', cond(K));

MsS = cellfun(@getPeakGain, S_mimo);
MsT = cellfun(@getPeakGain, T_mimo);
sv  = struct('sigS', num2cell(MsS), 'sigT', num2cell(MsT));
imprime_tabela('Valores singulares maximos', nomes, sv, {'sigS', 'sigT'}, ...
    {'Ms = max sigma(S)', 'Mt = max sigma(T)'}, '');

%% 8.5 Figuras de comparacao
estilo = {'--', '-', '-.', '-', '-'};
cores  = [0.55 0.25 0.60; 0.10 0.25 0.65; 0.85 0.45 0.10; 0.80 0.15 0.20; 0.20 0.60 0.30];

% modelo de referencia: resposta desejada pela sintonia lambda, 1/((lambda + t0)s + 1)
T_ref = [1/((lambda1 + t0_1)*s + 1) 0; 0 1/((lambda2 + t0_2)*s + 1)];
S_ref = eye(2) - T_ref;
y_ref = lsim(T_ref, r, t);

w = logspace(-3, 1, 1000);
sv_max = @(sys) max(sigma(sys, w), [], 1);

figure('Color','w', 'Name', 'Sensibilidade: máximo valor singular de S', 'NumberTitle', 'off');
semilogx(w, sv_max(S_ref), 'k:', 'LineWidth', 1.2); hold on;
for j = 1:ne
    semilogx(w, sv_max(S_mimo{j}), estilo{j}, 'Color', cores(j,:), 'LineWidth', 1.2);
end
grid on; xlabel('\omega [rad/s]'); ylabel('Máximo valor singular de S(j\omega)');
title('Sensibilidade');
legend(['Modelo de referência', arrayfun(@(j) sprintf('%s (M_s = %.3f)', nomes{j}, MsS(j)), 1:ne, 'UniformOutput', false)], 'Location', 'best');

figure('Color','w', 'Name', 'Sensibilidade complementar: máximo valor singular de T', 'NumberTitle', 'off');
semilogx(w, sv_max(T_ref), 'k:', 'LineWidth', 1.2); hold on;
for j = 1:ne
    semilogx(w, sv_max(T_mimo{j}), estilo{j}, 'Color', cores(j,:), 'LineWidth', 1.2);
end
grid on; xlabel('\omega [rad/s]'); ylabel('Máximo valor singular de T(j\omega)');
title('Sensibilidade complementar');
legend(['Modelo de referência', arrayfun(@(j) sprintf('%s (M_t = %.3f)', nomes{j}, MsT(j)), 1:ne, 'UniformOutput', false)], 'Location', 'best');

figure('Color','w', 'Name', 'Resposta ao degrau em malha fechada', 'NumberTitle', 'off');
sgtitle('Resposta ao degrau em malha fechada');
for m = 1:2
    subplot(2,1,m);
    plot(t, y_ref(:,m), 'k:', 'LineWidth', 1.2); hold on;
    for j = 1:ne
        plot(t, Y{j}(:,m), estilo{j}, 'Color', cores(j,:), 'LineWidth', 1.2);
    end
    grid on; ylabel(sprintf('y_%d', m)); title(sprintf('Malha %d — saída y_%d', m, m));
end
xlabel('Tempo [s]');
subplot(2,1,1); legend(['Modelo de referência', nomes], 'Location', 'best');

%% 8.6 Um grafico por indice, todos os metodos
% janelas (as mesmas da funcao indices):
%   seguimento: malha 1 em [t_r1, t_r2) | malha 2 em [t_r2, fim)
%   interacao : malha 1 em [t_r2, fim)  | malha 2 em [t_r1, t_r2)
Wseg = [t >= t_r1 & t < t_r2,  t >= t_r2];
Wint = [t >= t_r2,  t >= t_r1 & t < t_r2];

jan_seg = [t_r1, t_r2;  t_r2, Tsim];
jan_int = [t_r2, Tsim;  t_r1, t_r2];

tit_seg = {'Malha 1 — degrau em r_1', 'Malha 2 — degrau em r_2'};
tit_int = {'Malha 1 — efeito do degrau em r_2', 'Malha 2 — efeito do degrau em r_1'};
tit_u   = {'Entrada da planta u''_1', 'Entrada da planta u''_2'};

val = @(campo, S) reshape([S.(campo)], 2, []).';   % ne x 2

IAEc = cell(1,ne);  SIIc = cell(1,ne);  Dev = cell(1,ne);  ISUc = cell(1,ne);
for j = 1:ne
    e = r - Y{j};
    IAEc{j} = integral_acum(t, abs(e), Wseg);
    SIIc{j} = integral_acum(t, abs(e), Wint);
    ISUc{j} = integral_acum(t, U{j}.^2, true(size(e)));

    Dev{j} = Y{j} - r;   Dev{j}(~Wint) = NaN;           % desvio por interacao
end

figura_indice('IAE de seguimento (\int|e| dt acumulado)', t, IAEc, jan_seg, tit_seg, ...
    {'IAE_1', 'IAE_2'}, nomes, estilo, cores, 'IAE', val('IAE', ind), []);

figura_indice('SII (\int|e| dt acumulado na interação)', t, SIIc, jan_int, tit_int, ...
    {'SII_1', 'SII_2'}, nomes, estilo, cores, 'SII', val('SII', ind), []);

figura_indice('Desvio máximo por interação (D_{max})', t, Dev, jan_int, tit_int, ...
    {'y_1 - r_1', 'y_2 - r_2'}, nomes, estilo, cores, 'D_{max}', val('Dmax', ind), []);

figura_indice('Sobressinal (M_p)', t, Y, jan_seg, tit_seg, ...
    {'y_1', 'y_2'}, nomes, estilo, cores, 'M_p [%]', val('Mp', ind), r);

ax = figura_indice('Tempo de acomodação 2% (t_s)', t, Y, jan_seg, tit_seg, ...
    {'y_1', 'y_2'}, nomes, estilo, cores, 't_s [s]', val('ts', ind), r);
for m = 1:2
    yline(ax(m), amp_degrau*[0.98 1.02], 'k:', 'HandleVisibility', 'off');
end

figura_indice('Pico do sinal de controle |u''|', t, U, [0 Tsim; 0 Tsim], tit_u, ...
    {'u''_1', 'u''_2'}, nomes, estilo, cores, 'pico', val('pico', ef), []);

figura_indice('ISU (\int u''^2 dt acumulado)', t, ISUc, [0 Tsim; 0 Tsim], tit_u, ...
    {'ISU_1', 'ISU_2'}, nomes, estilo, cores, 'ISU', val('ISU', ef), []);

figure('Color','w', 'Name', 'Saída do controlador (u): todos os métodos', 'NumberTitle', 'off');
sgtitle('Saída do controlador (u)');
for m = 1:2
    subplot(2,1,m);
    for j = 1:ne
        plot(t, Cc{j}(:,m), estilo{j}, 'Color', cores(j,:), 'LineWidth', 1.2); hold on;
    end
    xline(t_r1, ':', 'r_1');   xline(t_r2, ':', 'r_2');
    grid on; box on; ylabel(sprintf('u_%d', m)); title(sprintf('Malha %d — saída do controlador u_%d', m, m));
    legend(nomes, 'Location', 'best');
end
xlabel('Tempo [s]');

%% 9. VARIAVEIS PARA O SIMULINK
[numG11, denG11] = tfdata(G11, 'v');
[numG12, denG12] = tfdata(G12, 'v');
[numG21, denG21] = tfdata(G21, 'v');
[numG22, denG22] = tfdata(G22, 'v');

% PID paralelo: P = Kc, D = Kc*tD, N = 1/(alpha*tD)  (sem derivada: N = 100)
Nf = @(tD) (tD > 0)/(alpha*max(tD, eps)) + (tD == 0)*100;

% desacoplador ideal: blocos D_12/D_21 (sinal - no somador)
[numD12, denD12] = tfdata(D12, 'v');
[numD21, denD21] = tfdata(D21, 'v');
P1_D = Kc1_prime;   I1_D = 0;   D1_D = Kc1_prime*tau1;   N1_D = Nf(tau1);
P2_D = Kc2_prime;   I2_D = 0;   D2_D = Kc2_prime*tau2;   N2_D = Nf(tau2);

% centralizado: linha i do controlador usa o PD da malha i
P11 = M(1,1);   I11 = 0;   D11   = tau_j(1)*M(1,1);   N11 = Nf(tau_j(1));
P12 = M(1,2);   I12 = 0;   D12_c = tau_j(1)*M(1,2);   N12 = Nf(tau_j(1));
P21 = M(2,1);   I21 = 0;   D21_c = tau_j(2)*M(2,1);   N21 = Nf(tau_j(2));
P22 = M(2,2);   I22 = 0;   D22   = tau_j(2)*M(2,2);   N22 = Nf(tau_j(2));

% desacoplador invertido: D12_inv recebe u'2 e soma em u'1 ; D21_inv recebe u'1 e soma em u'2
% (D12_inv estritamente proprio -> sem laco algebrico)
[numD12_inv, denD12_inv] = tfdata(D12_inv, 'v');
[numD21_inv, denD21_inv] = tfdata(D21_inv, 'v');
P1_inv = Kc1_inv;   I1_inv = 0;   D1_inv = Kc1_inv*tD1_inv;   N1_inv = Nf(tD1_inv);
P2_inv = Kc2_inv;   I2_inv = 0;   D2_inv = Kc2_inv*tD2_inv;   N2_inv = Nf(tD2_inv);

% descentralizado
P1_dec = Kc1_dec;   I1_dec = 0;   D1_dec = Kc1_dec*tD1_dec;   N1_dec = Nf(tD1_dec);
P2_dec = Kc2_dec;   I2_dec = 0;   D2_dec = Kc2_dec*tD2_dec;   N2_dec = Nf(tD2_dec);

% Morilla, forma 1: um PID por elemento k_ji(s) = k_i*D(j,i)(s), recebe e_i e soma em u_j
P11_pi = Kp_mor(1,1)*k_mor(1);   I11_pi = Ki_mor(1,1)*k_mor(1);   D11_pi = Kd_mor(1,1)*k_mor(1);
P12_pi = Kp_mor(1,2)*k_mor(2);   I12_pi = Ki_mor(1,2)*k_mor(2);   D12_pi = Kd_mor(1,2)*k_mor(2);
P21_pi = Kp_mor(2,1)*k_mor(1);   I21_pi = Ki_mor(2,1)*k_mor(1);   D21_pi = Kd_mor(2,1)*k_mor(1);
P22_pi = Kp_mor(2,2)*k_mor(2);   I22_pi = Ki_mor(2,2)*k_mor(2);   D22_pi = Kd_mor(2,2)*k_mor(2);
Nmor = @(j,i) (Kd_mor(j,i) ~= 0)*N_mor*abs(Kp_mor(j,i)/(Kd_mor(j,i) + (Kd_mor(j,i) == 0))) ...
              + (Kd_mor(j,i) == 0)*100;   % N = 1/Tf
N11_pi = Nmor(1,1);   N12_pi = Nmor(1,2);   N21_pi = Nmor(2,1);   N22_pi = Nmor(2,2);

% Morilla, forma 2 (modelo Método de Morilaa)
%   C11, C22 (PID, so P)         : ganho k_i da malha -> e_i*k_i
%   D11..D22 (Transfer Fcn)      : D(j,i)(s), recebe e_i*k_i e soma em u_j
%   atraso_dji (Transport Delay) : tau_pid(j,i), atraso de cada D(j,i)
%   atraso_gij (Transport Delay) : Tau_tot(i,j), atraso de cada Gij da planta
C11_P = k_mor(1);
C22_P = k_mor(2);
[numD11_mor, denD11_mor] = tfdata(D_mor{1,1}, 'v');
[numD12_mor, denD12_mor] = tfdata(D_mor{1,2}, 'v');
[numD21_mor, denD21_mor] = tfdata(D_mor{2,1}, 'v');
[numD22_mor, denD22_mor] = tfdata(D_mor{2,2}, 'v');

fprintf('\n--- Variaveis Simulink ---\n');
fprintf('C1_D: P=%.6f D=%.6f N=%.4f | C2_D: P=%.6f D=%.6f N=%.4f\n', P1_D, D1_D, N1_D, P2_D, D2_D, N2_D);
fprintf('Morilla: C11_P=%.6f | C22_P=%.6f\n', C11_P, C22_P);
fprintf('D11_mor: num = [%s] | den = [%s]\n', num2str(numD11_mor), num2str(denD11_mor));
fprintf('D12_mor: num = [%s] | den = [%s]\n', num2str(numD12_mor), num2str(denD12_mor));
fprintf('D21_mor: num = [%s] | den = [%s]\n', num2str(numD21_mor), num2str(denD21_mor));
fprintf('D22_mor: num = [%s] | den = [%s]\n', num2str(numD22_mor), num2str(denD22_mor));
fprintf('Centralizado: P11=%.6f P12=%.6f P21=%.6f P22=%.6f\n', P11, P12, P21, P22);
fprintf('              D11=%.6f D12_c=%.6f D21_c=%.6f D22=%.6f\n', D11, D12_c, D21_c, D22);
fprintf('D12_inv: num = [%s] | den = [%s]\n', num2str(numD12_inv), num2str(denD12_inv));
fprintf('D21_inv: num = [%s] | den = [%s]\n', num2str(numD21_inv), num2str(denD21_inv));
fprintf('C1_inv: P=%.6f D=%.6f N=%.4f | C2_inv: P=%.6f D=%.6f N=%.4f\n', ...
    P1_inv, D1_inv, N1_inv, P2_inv, D2_inv, N2_inv);

%% FUNCOES LOCAIS (precisam ficar no fim do script)
function C = pd_filtrado(Kc, tD, alpha)
    % PD com filtro na derivada: Kc*(1 + tD*s/(alpha*tD*s + 1))
    s = tf('s');
    if tD > 0
        C = Kc*(1 + tD*s/(alpha*tD*s + 1));
    else
        C = tf(Kc);
    end
end

function [y, u, up, L, S, T] = malha_fechada(G, P, C, r, t)
    % P = desacoplador ; u = saida do controlador ; up = entrada da planta
    L  = minreal(ss(G*P*C), [], false);
    S  = feedback(ss(eye(2)), L);
    T  = feedback(L, ss(eye(2)));
    y  = lsim(T, r, t);
    u  = lsim(C*S, r, t);
    up = lsim(P*C*S, r, t);
end

function [Kf, tau, t0] = fit3(t, y, t_deg, jan, du, Ts)
    % Fit3 (Cardenas, 2024) para processo integrador
    % ganho, eq. (2): K* = (s2 - s1)/du  (diferenca das inclinacoes antes/depois do degrau)
    % tau e t0, eq. (4): instantes de 28,3% e 63,2% da variacao da derivada da saida
    ii   = @(tt) round(tt/Ts) + 1;
    incl = @(a, b) (y(ii(b)) - y(ii(a)))/(b - a);
    s1 = incl(jan(1,1), jan(1,2));
    s2 = incl(jan(2,1), jan(2,2));
    Kf = (s2 - s1)/du;

    dy = [diff(y)/Ts; 0];
    N  = numel(t);
    k  = ii(t_deg);
    while k < N && abs(dy(k) - s1) < 0.283*abs(s2 - s1), k = k + 1; end
    t1 = t(k) - t_deg;
    while k < N && abs(dy(k) - s1) < 0.632*abs(s2 - s1), k = k + 1; end
    t2 = t(k) - t_deg;
    tau = 1.5*(t2 - t1);   t0 = t2 - tau;
end

function c = integral_acum(t, x, W)
    % integral acumulada de x so dentro da janela W (fora dela = NaN, nao plota)
    c = nan(size(x));
    for m = 1:size(x,2)
        c(W(:,m),m) = cumtrapz(t(W(:,m)), x(W(:,m),m));
    end
end

function ax = figura_indice(titulo, t, X, jan, tits, ylabs, nomes, estilo, cores, txt, val, ref)
    % uma figura por indice: 2 malhas, todos os metodos, valor do indice na legenda
    % ref = setpoint (ou [])
    figure('Color','w', 'Name', titulo, 'NumberTitle', 'off');
    sgtitle(titulo);
    ax = gobjects(1,2);
    for m = 1:2
        ax(m) = subplot(2,1,m);   hold on;
        leg = {};
        if ~isempty(ref)
            plot(t, ref(:,m), 'k-', 'LineWidth', 1.5);
            leg{end+1} = sprintf('setpoint r_%d', m);
        end
        for j = 1:numel(X)
            plot(t, X{j}(:,m), estilo{j}, 'Color', cores(j,:), 'LineWidth', 1.3);
            leg{end+1} = sprintf('%s (%s = %.3g)', nomes{j}, txt, val(j,m));
        end
        grid on;   box on;   xlim(jan(m,:));
        ylabel(ylabs{m});   title(tits{m});
        legend(leg, 'Location', 'best');
    end
    xlabel('Tempo [s]');
end

function ind = indices(t, y, r, t_r1, t_r2, amp)
    e  = r - y;
    w1 = t >= t_r1 & t < t_r2;    % janela do degrau em r1
    w2 = t >= t_r2;               % janela do degrau em r2

    ind.IAE  = [trapz(t(w1), abs(e(w1,1))), trapz(t(w2), abs(e(w2,2)))];
    ind.SII  = [trapz(t(w2), abs(e(w2,1))), trapz(t(w1), abs(e(w1,2)))];
    ind.Dmax = [max(abs(e(w2,1))), max(abs(e(w1,2)))];
    ind.Mp   = [max(0, (max(y(w1,1)) - amp)/amp*100), ...
                max(0, (max(y(w2,2)) - amp)/amp*100)];
    ind.ts   = [t_acomodacao(t(w1), y(w1,1), amp) - t_r1, ...
                t_acomodacao(t(w2), y(w2,2), amp) - t_r2];
end

function ts = t_acomodacao(t, y, ref)
    % primeiro instante a partir do qual |y - ref| fica dentro de 2% de ref
    k = find(abs(y - ref) > 0.02*abs(ref), 1, 'last');
    if isempty(k)
        ts = t(1);
    elseif k == numel(t)
        ts = NaN;                 % nao acomodou dentro da janela
    else
        ts = t(k+1);
    end
end

function ef = esforco(t, u)
    ef.pico = max(abs(u));
    ef.ISU  = trapz(t, u.^2);
end

function imprime_tabela(titulo, nomes, S, campos, rotulos, sufixo)
    fprintf('\n--- %s ---\n', titulo);
    fprintf('%-28s', 'Indice'); fprintf('%18s', nomes{:}); fprintf('\n');
    for c = 1:numel(campos)
        n = numel(S(1).(campos{c}));
        for m = 1:n
            if n > 1
                fprintf('%-28s', sprintf('%s %s %d', rotulos{c}, sufixo, m));
            else
                fprintf('%-28s', rotulos{c});
            end
            for j = 1:numel(S)
                fprintf('%18.4f', S(j).(campos{c})(m));
            end
            fprintf('\n');
        end
    end
end
