% lattice network -- superinfect

output_folder = "C:\Users\xiaoy\Work\primes\point_source\superinfect\graphs_layer";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows = 40;
cols = 40;
N = rows * cols;

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);

D = 0.005; r = 0.1; A = 0.1; K = 1;
beta_1 = 0.5; beta_2 = 0.4;
gamma_1 = 0.2; gamma_2 = 0.1;
alpha_1 = 0.01; alpha_2 = 0.05;
sigma = 0.9;

d11 = 0.5; d12 = 0.1; d13 = 0.1; d22 = 0.3; d33 = 0.1;

all_k1 = [1,2,2,2,2,3,3,3];
all_k2 = [1,2,1,1,2,2,1,3];
all_k3 = [1,1,1,2,2,1,2,3];
colors = lines(length(all_k1));

figI = figure('Visible', 'off'); hold on;
figJ = figure('Visible', 'off'); hold on;

num_t =500;
time_bet = 1;

for d_idx = 1:length(all_k1)
    d_idx

    k_1 = all_k1(d_idx);
    k_2 = all_k2(d_idx);
    k_3 = all_k3(d_idx);

    A1 = generate_LA2D(rows, cols, k_1);
    A2 = generate_LA2D(rows, cols, k_2);
    A3 = generate_LA2D(rows, cols, k_3);

    L1 = A1 - diag(sum(A1, 2)); 
    L2 = A2 - diag(sum(A2, 2)); 
    L3 = A3 - diag(sum(A3, 2));

    S0 = ones(N, 1);
    I0 = zeros(N, 1); J0 = zeros(N, 1);
    I0(1)=0.05; J0(1600)=0.05;

    y0 = [S0;I0;J0];

    nodes_I = zeros(1, num_t+1);
    nodes_J = zeros(1, num_t+1);
    time_vals = zeros(1, num_t+1);

    nodes_I(1) = sum(I0 > 0.01) / N;
    nodes_J(1) = sum(J0 > 0.01) / N;
    time_vals(1) = 0;

    for i = 1:num_t
 
        tspan = [0 time_bet];
        ode_wrap = @(t, y) superinfect_diff(t, y, L1, L2, L3, r,K,A,D,beta_1, beta_2, ...
            gamma_1, gamma_2, alpha_1, alpha_2, sigma, d11, d12, d13, d22, d33, N);
        [t, y] = ode113(ode_wrap, tspan, y0);

        S = y(end, 1:N)';
        I = y(end, N+1:2*N)';
        J = y(end, 2*N+1:3*N)';
    
        S = max(S, 0);
        I = max(I, 0);
        J = max(J, 0);

        nodes_I(i+1) = sum(I > 0.01) / N;
        nodes_J(i+1) = sum(J > 0.01) / N;

        time_vals(i+1) = i * time_bet;

        y0 = [S;I;J];
    end

    figure(figI); 
    plot(time_vals, nodes_I, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figJ); 
    plot(time_vals, nodes_J, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);
end

figure(figI); xlabel('Time'); ylabel('I-Spread Index');
legend(compose('(LA%d,LA%d,LA%d)', [all_k1(:), all_k2(:), all_k3(:)]), 'Location', 'best');
saveas(figI, fullfile(output_folder, 'la_I_v1'));
close(figI);

figure(figJ); xlabel('Time'); ylabel('J-Spread Index');
legend(compose('(LA%d,LA%d,LA%d)', [all_k1(:), all_k2(:), all_k3(:)]), 'Location', 'best');
saveas(figJ, fullfile(output_folder, 'la_J_v1'));
close(figJ);

function dydt = superinfect_diff(t,y,L1,L2,L3,r,K,A,D,beta_1,beta_2,gamma_1,gamma_2,alpha_1,alpha_2,sigma,d11,d12,d13,d22,d33,N)
    S = y(1:N);
    I = y(N+1:2*N);
    J = y(2*N+1:end);

    f = r*S.*(1-S/K).*(S/A-1)-(beta_1*I+beta_2*J).*S./(S+I+J)+gamma_1*I+gamma_2*J-D*S;
    g = I.*(beta_1*S./(S+I+J)-alpha_1-gamma_1-sigma*beta_2*J./(S+I+J))-D*I;
    h = J.*(beta_2*S./(S+I+J)-alpha_2-gamma_2+sigma*beta_2*I./(S+I+J))-D*J;

    dSdt = f+d11*(L1*S)+d12*(L2*I)+d13*(L3*J);
    dIdt = g+d22*(L2*I);
    dJdt = h+d33*(L3*J);

    dydt = [dSdt; dIdt; dJdt];
end


function A = generate_LA2D(rows, cols, range)
    N = rows * cols;
    A = zeros(N);

    for idx = 1:N
        yi = ceil(idx / cols);
        xi = idx - (yi - 1) * cols;

        neighbors = [];
        for dx = -range:range
            for dy = -range:range
                if dx == 0 && dy == 0
                    continue;
                end
                d2 = dx^2 + dy^2;
                if d2 <= range^2 && d2 ~= 8
                    xj = xi + dx;
                    yj = yi + dy;
                    if xj >= 1 && xj <= cols && yj >= 1 && yj <= rows
                        jdx = (yj - 1) * cols + xj;
                        neighbors(end+1) = jdx;
                    end
                end
            end
        end

        neighbors = sort(neighbors);
        A(idx, neighbors) = 1;
    end
end