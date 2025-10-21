output_folder = "C:\Users\xiaoy\Work\primes\point_source\coinfect\graphs_raw";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows = 40;
cols = 40;
N = rows * cols;

k_1=1; k_2=1; k_3=1;

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);

L1 = A1 - diag(sum(A1, 2)); 
L2 = A2 - diag(sum(A2, 2)); 
L3 = A3 - diag(sum(A3, 2));

D=0.005; r=0.1; A=0.1; K=1;
beta_1 = 0.3; beta_2 = 0.4; beta_10 = 0.2; beta_02 = 0.3;
gamma_1 = 0.1; gamma_2 = 0.05;
alpha_1 = 0.05; alpha_2 = 0.15;
d11 = 0.3; d12 = 0.1; d13 = 0.1; d22 = 0.3; d33 = 0.1;

params.L1 = L1; params.L2 = L2; params.L3 = L3;
params.r = r; params.K = K; params.A = A; params.D = D;
params.beta_1 = beta_1; params.beta_2 = beta_2;
params.beta_10 = beta_10; params.beta_02 = beta_02;
params.gamma_1 = gamma_1; params.gamma_2 = gamma_2;
params.alpha_1 = alpha_1; params.alpha_2 = alpha_2;
params.d11 = d11; params.d12 = d12; params.d13 = d13;
params.d22 = d22; params.d33 = d33;
params.N = N;
params.num_t = 1000;
params.time_bet = 1;

S0 = ones(N, 1); I0 = zeros(N, 1); J0 = zeros(N, 1); C0 = zeros(N, 1);
I0(1)=0.05; J0(1600)=0.05;
params.y0 = [S0; I0; J0; C0];

alpha12_vals = 0.36:0.02:0.66;  
b12_thresholds = zeros(size(alpha12_vals));

for idx = 1:length(alpha12_vals)
    fprintf('Computing for alpha_12 = %.3f...\n', alpha12_vals(idx));
    params.alpha_12 = alpha12_vals(idx);   
    b12_thresholds(idx) = find_b12_threshold([0, 0.5], 1e-3, 20, params);
end

figure;
scatter(alpha12_vals, b12_thresholds, 40, 'filled');
xlabel('\alpha_{12}');
ylabel('b_{12} Threshold');
filename = sprintf('b12thres_vs_alpha_new');
savefig(fullfile(output_folder, filename));

function b12_thresh = find_b12_threshold(b12_range, tol, max_iter, params)
    low = b12_range(1);
    high = b12_range(2);
    b12_thresh = low;

    for iter = 1:max_iter
        mid = (low + high) / 2;
        [C_peak, C_final] = simulate_coinfect_C(mid, params);
        fprintf('Iter %d: b12=%.4f, C_peak=%.3f, C_final=%.3f\n', iter, mid, C_peak, C_final);

        %if C_peak < 1 && C_final < 0.1
        if C_final ==0
            b12_thresh = mid;
            low = mid;
        else
            high = mid;
        end

        if abs(high - low) < tol
            break;
        end
    end
end

function [C_peak, C_final] = simulate_coinfect_C(beta_12, params)
    % Run the coinfection model and return C peak and final values

    y0 = params.y0;
    num_t = params.num_t;
    time_bet = params.time_bet;

    nodes_C = zeros(1, num_t+1);
    nodes_C(1) = sum(y0(3*params.N+1:end) > 0.005) / params.N;

    for i = 1:num_t
        tspan = [0 time_bet];
        ode_wrap = @(t, y) coinfect_diff(t,y,params.L1,params.L2,params.L3, ...
            params.r, params.K, params.A, params.D, params.beta_1, params.beta_2, ...
            params.beta_10, params.beta_02, beta_12, params.gamma_1, params.gamma_2, ...
            params.alpha_1, params.alpha_2, params.alpha_12, ...
            params.d11, params.d12, params.d13, params.d22, params.d33, params.N);
        [~, y] = ode113(ode_wrap, tspan, y0);

        C = y(end, 3*params.N+1:end)';
        C = max(C, 0);
        nodes_C(i+1) = sum(C > 0.005) / params.N;

        y0 = y(end,:)';
    end

    C_peak = max(nodes_C);
    C_final = nodes_C(end);
end

function dydt = coinfect_diff(t,y,L1,L2,L3,r, K, A,D,beta_1,beta_2,beta_10,beta_02,beta_12,gamma_1,gamma_2,alpha_1,alpha_2,alpha_12,d11,d12,d13,d22,d33,N)
    S = y(1:N);
    I = y(N+1:2*N);
    J = y(2*N+1:3*N);
    C = y(3*N+1:end);

    f = r*S.*(1-S/K).*(S/A-1)-(beta_1*I + beta_2*J).*S./(S+I+J-C) ...  
    -(beta_10+beta_02+beta_12-beta_1-beta_2)*C.*S./(S+I+J-C)+gamma_1*I+gamma_2*J-(gamma_1+gamma_2)*C-D*S;
    g = (beta_1*I+(beta_10+beta_12-beta_1)*C).*(S+J-C)./(S+I+J-C)-gamma_1*I-alpha_1*(I-C)-alpha_12*C-D*I;
    h = (beta_2*J+(beta_02+beta_12-beta_2)*C).*(S+I-C)./(S+I+J-C)-gamma_2*J-alpha_2*(J-C)-alpha_12*C-D*J;
    l = beta_12*S.*C./(S+I+J-C)+(beta_2*J+(beta_02+beta_12-beta_2)*C).*(I-C)./(S+I+J-C)...
    +(beta_1*I+(beta_10+beta_12-beta_1)*C).*(J-C)./(S+I+J-C)-(gamma_1+gamma_2+alpha_12)*C-D*C;

    dSdt = f+d11*(L1*S)+d12*(L2*I)+d13*(L3*J);
    dIdt = g+d22*(L2*I);
    dJdt = h+d33*(L3*J);
    dCdt = l;

    dydt = [dSdt; dIdt; dJdt; dCdt];
end


function A = generate_LA2D(rows, cols, range)
    N = rows * cols;
    A = zeros(N);
    [X, Y] = meshgrid(1:cols, 1:rows);
    coords = [X(:), Y(:)];
    for idx = 1:N
        xi = coords(idx, 1);
        yi = coords(idx, 2);
        for dx = -range:range
            for dy = -range:range
                if dx == 0 && dy == 0, continue; end
                if sqrt(dx^2 + dy^2) <= range
                    xj = xi + dx; yj = yi + dy;
                    if xj >= 1 && xj <= cols && yj >= 1 && yj <= rows
                        jdx = (yj - 1) * cols + xj;
                        A(idx, jdx) = 1;
                    end
                end
            end
        end
    end
end