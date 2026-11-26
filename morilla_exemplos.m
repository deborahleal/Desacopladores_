clear; clc; close all;
s = tf('s');

%% 1) PLANTA

% EXEMPLO 1: Wood-Berry
% G = [ 12.8/(16.7*s+1)*exp(-1*s),   -18.9/(21.0*s+1)*exp(-3*s);
%        6.6/(10.9*s+1)*exp(-7*s),   -19.4/(14.4*s+1)*exp(-3*s) ];

% EXEMPLO 2: Quatro tanques
% T1 = 1; T2 = 1.5; T3 = 0.5; T4 = 0.5;     % constantes de tempo
% alfa1 = 3.7;  alfa2 = 4.7;
% gama1 = 0.5;  gama2 = 0.3;                % gama1+gama2 < 1 -> fase não mínima
% G = [ gama1*alfa1/(T1*s+1),                 (1-gama2)*alfa1/((T1*s+1)*(T3*s+1));
%       (1-gama1)*alfa2/((T2*s+1)*(T4*s+1)),  gama2*alfa2/(T2*s+1) ];

% EXEMPLO 3: Twin Rotor MIMO System, eq. (45)
% y1 = guinada (yaw), y2 = arfagem (pitch); u1, u2 = tensões nos rotores
G = [ 34.88/(1.43*s^3 + 1.07*s^2 + 2.05*s + 1.40),   13.36*s/(0.38*s^3 + 1.02*s^2 + 0.59*s + 1.40);
      (19.88*s + 13.50)/((8.43*s^2 + 7.26*s + 1)*s),  342.90/((2.24*s^2 + 6.21*s + 1)*s) ];

%% 2) PARÂMETROS DE PROJETO
%   Caso 1 -> lbar_i = 1
%   Caso 2 -> lbar_i = -(s - z)/(s + z)
%   Caso 3 -> lbar_i = (s + z_i)/s

% EXEMPLO 1: Wood-Berry, Caso 1 nas duas malhas
% caso = [1 1];
% z_l  = [0 0];
% lbar = {tf(1), tf(1)};

%  EXEMPLO 2: Quatro tanques, Caso 2 nas duas malhas
% caso = [2 2];
% z_l  = [1.06 1.06];                % zero no SPD
% lbar = {-(s - z_l(1))/(s + z_l(1)), -(s - z_l(2))/(s + z_l(2))};

% EXEMPLO 3: Twin Rotor, Caso 1 na malha 1 e Caso 3 na malha 2
% (o zero no SPD em z ~ 6,38 é desprezado: está longe dos polos dominantes)
caso = [1 3];
z_l  = [0 1/20];                   % z_2 = 1/20 (eq. 30)
lbar = {tf(1), (s + z_l(2))/s};

% EXEMPLO 1 e EXEMPLO 2: PI
% ordem = 1;                        % 1 = PI ; 2 = PID
% EXEMPLO 3: PID (o artigo diz que PI não é indicado aqui, modos oscilatórios)
ordem = 2;
N_filtro = 10;                      % filtro da derivada

% Elementos do desacoplador a desprezar (zerar) depois do ajuste
% EXEMPLO 1 e EXEMPLO 2: nenhum
% zerar = false(2);
% EXEMPLO 3: o artigo despreza D12 na eq. (46) (desacoplamento parcial)
zerar = [false true; false false];

% Faixa de frequência para o ajuste de m_ji
w_min = 1e-3;
N_w   = 300;
% w_max = 0.3;                      % EXEMPLO 1
% w_max = 1;                        % EXEMPLO 2 (ordem de grandeza do zero z)
w_max = 0.3;                        % EXEMPLO 3 (reproduz a eq. 46)

% Especificação para os ganhos k_i (um critério por malha):
%   'Am'    -> margem de ganho Am_i > 1          (Caso 1 com atraso)
%   'phim'  -> margem de fase  phim_i < 90 graus  (Caso 1 com atraso)
%   'T'     -> constante de tempo de MF T_i       (Caso 1 sem atraso)
%   'delta' -> amortecimento delta_i              (Casos 2 e 3 sem atraso)

%EXEMPLO 1: Am1 = 5, Am2 = 3
% criterio = {'Am', 'Am'};
% espec    = [5 3];

% EXEMPLO 2: delta1 = delta2 = 0,4
% criterio = {'delta', 'delta'};
% espec    = [0.4 0.4];

% EXEMPLO 3: T1 = 20 (-> k1 = 0,05, eq. 22) e delta2 = 1 (-> k2 = 0,2, eq. 30)
criterio = {'T', 'delta'};
espec    = [20 1];

% Simulação (r1 em t = 0 e r2 em t = t_r2)

% EXEMPLO 1
% t_final = 120;  t_r2 = 60;  unid_tempo = 'min';

%  EXEMPLO 2
% t_final = 60;   t_r2 = 30;  unid_tempo = 'min';

% EXEMPLO 3 (Fig. 8)
t_final = 300;  t_r2 = 150; unid_tempo = 's';

dt = 0.05;

%% 3) ATRASOS
theta = totaldelay(G);              % atraso total de cada elemento g_ij

% atraso do numerador comum n(s), eq. (12)
tau_n = min(theta(1,1) + theta(2,2), theta(1,2) + theta(2,1));

% atraso de cada processo equivalente g~_ij, eqs. (13)-(14)
tau_til = zeros(2);
for i = 1:2
    for j = 1:2
        tau_til(i,j) = tau_n - theta(3-i, 3-j);
    end
end

% atraso da malha aberta l_i, eq. (16)
tau_i = max(tau_til, [], 2);

% atraso de cada controlador, eq. (32)
tauK = zeros(2);                    % tauK(j,i) = atraso de k_ji
for i = 1:2
    for j = 1:2
        tauK(j,i) = tau_i(i) - tau_til(i,j);
    end
end

fprintf(' ATRASOS \n');
fprintf('tau(n) = %.4g\n', tau_n);
fprintf('tau~_ij =\n'); disp(tau_til);
fprintf('tau_i = [%.4g  %.4g]\n', tau_i);
fprintf('tau(k_ji) (linha j, coluna i) =\n'); disp(tauK);

%% 4) PROCESSOS EQUIVALENTES g~_ij
w   = logspace(log10(w_min), log10(w_max), N_w).';
jw  = 1j*w;
Gw  = freqresp(G, w);                          % 2 x 2 x N_w
g11 = squeeze(Gw(1,1,:));  g12 = squeeze(Gw(1,2,:));
g21 = squeeze(Gw(2,1,:));  g22 = squeeze(Gw(2,2,:));
detG = g11.*g22 - g12.*g21;
Gel  = {g11, g12; g21, g22};

gtil = cell(2);
for i = 1:2
    for j = 1:2
        if i == j, sinal = 1; else, sinal = -1; end
        gtil{i,j} = sinal*detG ./ Gel{3-i, 3-j};
    end
end

%% 4b) FUNÇÕES DE TRANSFERÊNCIA DOS PROCESSOS SISO EQUIVALENTES g~_ij
% Gr = parte racional de G (sem atrasos), tf 2x2. Também é usada no
% Simulink: blocos g_ij -> Gr.Numerator{i,j} / Gr.Denominator{i,j}
[Nc_g, Dc_g] = tfdata(G);           % (Dc_g: o nome Dc fica reservado ao Simulink)
Gr = tf(Nc_g, Dc_g);
a_d = theta(1,1) + theta(2,2);
b_d = theta(1,2) + theta(2,1);
P1  = Gr(1,1)*Gr(2,2);
P2  = Gr(1,2)*Gr(2,1);

gtil_sys = cell(2);                 % g~_ij como sistema (com atrasos)
fprintf('\n FUNÇÕES DE TRANSFERÊNCIA SISO EQUIVALENTES g~_ij(s)\n');
for i = 1:2
    for j = 1:2
        if i == j, sinal = 1; else, sinal = -1; end
        gk  = Gr(3-i, 3-j);
        thk = theta(3-i, 3-j);
        R1  = minreal(sinal*P1/gk, 1e-6);   d1 = a_d - thk;
        R2  = minreal(sinal*P2/gk, 1e-6);   d2 = b_d - thk;

        if abs(d1 - d2) < 1e-9
            % mesmo atraso nos dois termos: g~_ij = (R1 - R2) e^{-d1 s}
            Rij = minreal(R1 - R2, 1e-6);
            gtil_sys{i,j} = Rij*exp(-d1*s);
            fprintf('\ng~_%d%d(s) = e^{-%.4g s} *', i, j, d1);
            disp(Rij);
        else
            % atrasos diferentes: g~_ij = R1 e^{-d1 s} - R2 e^{-d2 s}
            gtil_sys{i,j} = ss(R1*exp(-d1*s)) - ss(R2*exp(-d2*s));
            fprintf('\ng~_%d%d(s) = e^{-%.4g s}*[ A(s) - B(s) e^{-%.4g s} ], com\n', ...
                    i, j, min(d1,d2), abs(d1 - d2));
            if d1 <= d2
                fprintf('A(s) ='); disp(R1); fprintf('B(s) ='); disp(R2);
            else
                fprintf('A(s) = -1 *'); disp(R2); fprintf('B(s) = -1 *'); disp(R1);
            end
        end
    end
end

% Conferência: g~_ij calculado como função de transferência x resposta em
% frequência da seção 4 (o erro deve ser ~0)
erro_g = zeros(2);
for i = 1:2
    for j = 1:2
        gw_tf = squeeze(freqresp(gtil_sys{i,j}, w));
        erro_g(i,j) = max(abs(gw_tf - gtil{i,j}) ./ abs(gtil{i,j}));
    end
end
fprintf('\nErro relativo máximo entre g~_ij (FT) e g~_ij (freq.):\n'); disp(erro_g);

%% 5) m_ji E AJUSTE DE BAIXA ORDEM
Kp = zeros(2); Ki = zeros(2); Kd = zeros(2);
m_real = cell(2); m_aprox = cell(2);

for i = 1:2
    lbar_w = squeeze(freqresp(lbar{i}, w));
    for j = 1:2
        m = gtil{i,j} .* exp(jw*tau_til(i,j)) ./ lbar_w;

        if ordem == 1
            A = [ones(N_w,1), -m];
        else
            A = [ones(N_w,1), -m, -m.*jw.^2];
        end
        b = m.*jw;
        x = [real(A); imag(A)] \ [real(b); imag(b)];

        b0 = x(1); a0 = x(2); a1 = 1;
        if ordem == 2, a2 = x(3); else, a2 = 0; end

        % Ganhos do PID do elemento k_ji (sem o k_i): eq. (34)
        Kp(j,i) = a1/b0;
        Ki(j,i) = a0/b0;
        Kd(j,i) = a2/b0;

        m_real{j,i}  = m;
        m_aprox{j,i} = b0 ./ (a2*jw.^2 + a1*jw + a0);
    end
end

%% 6) GANHOS k_i
if ischar(criterio), criterio = {criterio, criterio}; end   % mesmo critério nas duas malhas
k = zeros(1,2);
for i = 1:2
    switch criterio{i}
        case 'Am'                       % eq. (21)
            k(i) = pi/(2*espec(i)*tau_i(i));
        case 'phim'                     % eq. (20)
            k(i) = pi*(90 - espec(i))/(180*tau_i(i));
        case 'T'                        % eq. (22): T_i = 1/k_i
            k(i) = 1/espec(i);
        case 'delta'
            if caso(i) == 2             % eq. (26): delta = (z-k)/(2 sqrt(k z))
                k(i) = z_l(i)*(sqrt(1 + espec(i)^2) - espec(i))^2;
            elseif caso(i) == 3         % eq. (30): delta = sqrt(k/(4 z))
                k(i) = 4*z_l(i)*espec(i)^2;
            else
                error('criterio ''delta'' vale só para os Casos 2 e 3.');
            end
    end
    if tau_i(i) == 0 && any(strcmp(criterio{i}, {'Am','phim'}))
        error('tau_%d = 0: use criterio ''T'' (Caso 1) ou ''delta'' (Casos 2 e 3).', i);
    end
end
fprintf('\n GANHOS\n');
fprintf('k1 = %.4f   k2 = %.4f\n', k);

%% 6b) FUNÇÕES DE TRANSFERÊNCIA DAS MALHAS SISO PROJETADAS
l_esp = cell(1,2);  h_esp = cell(1,2);
fprintf('\n MALHAS SISO PROJETADAS\n');
for i = 1:2
    l_esp{i} = minreal(k(i)*lbar{i}/s, 1e-6)*exp(-tau_i(i)*s);
    h_esp{i} = feedback(ss(l_esp{i}), 1);       % ss: aceita atraso na malha
    fprintf('\nl_%d(s) =', i); disp(l_esp{i});
    if tau_i(i) == 0
        hi = minreal(tf(h_esp{i}), 1e-6);
        [wn_i, zeta_i] = damp(hi);
        fprintf('h_%d(s) =', i); disp(hi);
        fprintf('polos de h_%d: wn = %s, delta = %s\n', i, ...
                mat2str(unique(round(wn_i, 4)).'), mat2str(unique(round(zeta_i, 4)).'));
    else
        fprintf('h_%d(s) = l_%d/(1 + l_%d)  (tem atraso na malha: não é racional)\n', i, i, i);
    end
end

%% 7) DESACOPLADOR D^PID(s) E CONTROLADOR
D  = tf(zeros(2));
Dc = cell(2);                       % D_ji sem atraso (parte racional) -> Simulink
fprintf('\n DESACOPLADOR D^PID(s)\n');
for j = 1:2
    for i = 1:2
        if zerar(j,i)
            D(j,i) = 0;
            Dc{j,i} = tf(0);
            fprintf('D%d%d(s) ~ 0  (desprezado)\n', j, i);
            continue
        end
        Dji = Kp(j,i) + Ki(j,i)/s;
        if Kd(j,i) ~= 0
            Tf  = abs(Kd(j,i)/(Kp(j,i)*N_filtro));
            Dji = Dji + Kd(j,i)*s/(Tf*s + 1);
        end
        Dc{j,i} = minreal(tf(Dji), 1e-9);
        D(j,i) = Dji*exp(-tauK(j,i)*s);

        if ordem == 1
            % (s + Ki/Kp)/(s/Kp) * e^{-tau s}
            fprintf('D%d%d(s) = (s + %.4f)/(%.4f s) * e^{-%.3g s}\n', ...
                    j, i, Ki(j,i)/Kp(j,i), 1/Kp(j,i), tauK(j,i));
        else
            % (s^2 + Kp/Kd s + Ki/Kd)/(s/Kd) * e^{-tau s}  (formato da eq. 46)
            fprintf('D%d%d(s) = (s^2 + %.4g s + %.4g)/(%.4g s) * e^{-%.3g s}\n', ...
                    j, i, Kp(j,i)/Kd(j,i), Ki(j,i)/Kd(j,i), 1/Kd(j,i), tauK(j,i));
        end
    end
end

%  EXEMPLO 1
% fprintf(['\nArtigo, eq. (41):\n', ...
%  'D11 = (s+0.13)/(0.81s)          D12 = -(s+0.23)/(1.50s) e^{-2s}\n', ...
%  'D21 = (s+0.085)/(1.57s) e^{-4s} D22 = -(s+0.165)/(1.57s)\n']);

% EXEMPLO 2: eq. (44)
% fprintf(['\nArtigo, eq. (44):\n', ...
%  'D11 = -(s+1.09)/(5.62s)         D12 = (s+1.58)/(4.15s)\n', ...
%  'D21 = (s+4.01)/(11.56s)         D22 = -(s+0.74)/(2.88s)\n']);

% EXEMPLO 3: eq. (46)
fprintf(['\nArtigo, eq. (46):\n', ...
 'D11 = (s^2+1.82s+1.24)/(30.79s)    D12 = -(s^2+0.037s+2.66e-11)/(656s) ~ 0\n', ...
 'D21 = -(s^2+1.36s+0.72)/(454s)     D22 = (s^2+0.206s+0.0079)/(54.1s)\n']);

% Verificação em regime: para baixas frequências K = G^-1 L leva a
%   Ki (matriz dos ganhos integrais) = G(0)^-1 * diag(lbar_1(0), lbar_2(0))
% (só quando G(0) e lbar_i(0) são finitos -> não se aplica ao Exemplo 3)
G0 = dcgain(G);
lbar0 = [dcgain(lbar{1}), dcgain(lbar{2})];
if all(isfinite(G0(:))) && all(isfinite(lbar0))
    Ki_esperado = G0 \ diag(lbar0);
    fprintf('\nVerificação em regime (Ki do ajuste x Ki exato):\n');
    disp(table(Ki(:), Ki_esperado(:), 'VariableNames', {'Ki_ajuste','Ki_exato'}, ...
         'RowNames', {'D11','D21','D12','D22'}));
end

K = D*diag(k);                      % controlador centralizado completo

%% 8a) GRÁFICO: aproximação de m_ji
figure('Name', 'Aproximação de m_{ji}');
for i = 1:2
    for j = 1:2
        subplot(2, 2, (j-1)*2 + i);
        plot(real(m_real{j,i}),  imag(m_real{j,i}),  'b-',  'LineWidth', 1.4); hold on;
        plot(real(m_aprox{j,i}), imag(m_aprox{j,i}), 'r--', 'LineWidth', 1.4);
        grid on; axis equal;
        xlabel('Real'); ylabel('Imaginário');
        title(sprintf('m_{%d%d}(j\\omega)', j, i));
        legend('real', 'aproximado', 'Location', 'best');
    end
end

%% 8b) SIMULAÇÃO EM MALHA FECHADA (Fig. 4 no Ex. 1; Fig. 6 no Ex. 2; Fig. 8 no Ex. 3)
L  = ss(G)*ss(K);                   % ss suporta os atrasos internos
CL = feedback(L, eye(2));

t = (0:dt:t_final).';
r = [ones(size(t)), double(t >= t_r2)];
y = lsim(CL, r, t);

% resposta especificada de cada malha SISO, h_i (seção 6b)
y_esp = zeros(length(t), 2);
for p = 1:2
    y_esp(:,p) = lsim(h_esp{p}, r(:,p), t);
end

figure('Name', 'Resposta ao degrau - malha fechada');
for p = 1:2
    subplot(2, 1, p);
    plot(t, r(:,p), 'k:', 'LineWidth', 1.0); hold on;
    plot(t, y_esp(:,p), 'r--', 'LineWidth', 1.2);
    plot(t, y(:,p), 'b-', 'LineWidth', 1.4);
    grid on;
    xlabel(['Tempo (' unid_tempo ')']); ylabel(sprintf('y_%d', p));
    title(sprintf('Saída y_%d', p));
    legend(sprintf('r_%d', p), sprintf('especificado h_%d', p), ...
           sprintf('y_%d (TITO)', p), 'Location', 'best');
end

%% 8c) Malhas abertas resultantes
figure('Name', 'Malha aberta');
bodemag(L, {w_min, 10}); grid on;
title('|L_{ij}(j\omega)| - termos fora da diagonal = acoplamento residual');
%% 9) SIMULINK

modelo  = 'morilla_test_simu';

k_i     = k;                        % ganhos proporcionais k1, k2
tau_pid = tauK;                     % atrasos dos elementos do desacoplador
Tau_tot = theta;                    % atrasos totais de cada g_ij
% Dc e Gr já foram criados nas seções 7 e 4b

load_system(modelo);
set_param(modelo, 'StopTime', num2str(t_final));
sim_out = sim(modelo, 'SaveOutput', 'on', 'OutputSaveName', 'yout', ...
              'SaveFormat', 'Dataset');

% Comparação Simulink x lsim (precisa de blocos Outport: Out1 = y1, Out2 = y2)
yout = [];
try, yout = sim_out.yout; catch, end
if isempty(yout) || yout.numElements < 2
    fprintf(['\nSimulink rodou. Para comparar com o lsim, adicione Out1 = y1 e ', ...
             'Out2 = y2 no modelo (ou veja os Scopes).\n']);
else
    figure('Name', 'Simulink x lsim');
    for p = 1:2
        ts = yout.getElement(p).Values;
        subplot(2, 1, p);
        plot(t, y(:,p), 'b-', 'LineWidth', 1.6); hold on;
        plot(ts.Time, squeeze(ts.Data), 'r--', 'LineWidth', 1.3);
        grid on;
        xlabel(['Tempo (' unid_tempo ')']); ylabel(sprintf('y_%d', p));
        title(sprintf('Saída y_%d', p));
        legend('lsim (script)', 'Simulink', 'Location', 'best');
    end
end