output_folder = ""; % fill in folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows = 40; cols = 40; N = rows * cols;
k_1=1; k_2=1; k_3=1;

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);

D = 0.005; r = 0.1; A = 0.1; K = 1;
beta_1 = 0.5; beta_2 = 0.4;
gamma_1 = 0.2; gamma_2 = 0.1;
alpha_1 = 0.01; alpha_2 = 0.05;

d11 = 0.3; d12 = 0.1; d13 = 0.1; d22 = 0.3; d33 = 0.1;

L1 = A1 - diag(sum(A1, 2)); 
L2 = A2 - diag(sum(A2, 2)); 
L3 = A3 - diag(sum(A3, 2));

S0 = ones(N, 1);
I0 = zeros(N, 1); J0 = zeros(N, 1);
I0(1)=0.05; J0(1600)=0.05;

all_sigma = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4];

num_t = 500;
time_bet = 1;

I_peaks = zeros(1, length(all_sigma));
I_peak_times = zeros(1, length(all_sigma));
J_sat_times = zeros(1, length(all_sigma));

for d_idx = 1:length(all_sigma)

    d_idx

    sigma = all_sigma(d_idx);

    y0 = [S0;I0;J0];
    nodes_I = zeros(1, num_t+1);
    nodes_J = zeros(1, num_t+1);
    time_vals = zeros(1, num_t+1);

    nodes_I(1) = sum(I0 > 0.01) / N;
    nodes_J(1) = sum(J0 > 0.01) / N;
    time_vals(1) = 0;

    sat_time = NaN;

    for i = 1:num_t
        tspan = [0 time_bet];
        ode_wrap = @(t, y) superinfect_diff(t, y, L1, L2, L3, r,K,A,D,beta_1, beta_2, ...
            gamma_1, gamma_2, alpha_1, alpha_2, all_sigma(d_idx), d11, d12, d13, d22, d33, N);
        [t, y] = ode113(ode_wrap, tspan, y0);

        S = max(y(end, 1:N)', 0);
        I = max(y(end, N+1:2*N)', 0);
        J = max(y(end, 2*N+1:end)', 0);

        nodes_I(i+1) = sum(I > 0.01) / N;
        nodes_J(i+1) = sum(J > 0.01) / N;
        time_vals(i+1) = i * time_bet;

        if isnan(sat_time) && nodes_J(i+1) >= 1.0
            sat_time = time_vals(i+1);
        end

        y0 = [S; I; J];
    end

    [I_peaks(d_idx), idx_peak] = max(nodes_I);
    I_peak_times(d_idx) = time_vals(idx_peak);

    J_sat_times(d_idx) = sat_time;
end

figure;
scatter(all_sigma, I_peaks, 40, 'filled');
xlabel('\sigma');
ylabel('I_1 Peak Value');
filename = sprintf('Ipeak_vs_sigma_superinfect');
savefig(fullfile(output_folder, filename));

figure;
scatter(all_sigma, I_peak_times, 40, 'filled');
xlabel('\sigma');
ylabel('I_1 Peak Time');
filename = sprintf('Ipeaktime_vs_sigma_superinfect');
savefig(fullfile(output_folder, filename));

figure;
scatter(all_sigma, J_sat_times, 40, 'filled');
xlabel('\sigma');
ylabel('I_2 Saturation Time');
filename = sprintf('Jsat_vs_sigma_superinfect');
savefig(fullfile(output_folder, filename));

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
