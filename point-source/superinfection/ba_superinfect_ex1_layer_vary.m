% ba network -- superinfect

output_folder = "C:\Users\xiaoy\Work\primes\point_source\superinfect\graphs_layer";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

N=1600;

D = 0.005; r = 0.1; A = 0.1; K = 1;
beta_1 = 0.5; beta_2 = 0.4;
gamma_1 = 0.2; gamma_2 = 0.1;
alpha_1 = 0.01; alpha_2 = 0.05;
sigma = 0.9;

d11 = 0.3; d12 = 0.1; d13 = 0.1; d22 = 0.3; d33 = 0.1;

all_m0_1 = [5,13,13,13,13,25,25,25];
all_m_1 = [2,6,6,6,6,12,12,12];
all_m0_2 = [5,13,5,5,13,13,5,25];
all_m_2 = [2,6,2,2,6,6,2,12];
all_m0_3 = [5,5,5,13,13,5,13,25];
all_m_3 = [2,2,2,6,6,2,6,12];
colors = lines(length(all_k1));

num_t =200;
time_bet = 1;
run_time = 1000;

avg_I_peak = zeros(length(all_k1),1);
avg_I_peak_time = zeros(length(all_k1),1);
avg_J_peak = zeros(length(all_k1),1);
avg_J_peak_time = zeros(length(all_k1),1);

for d_idx = 1:length(all_k1)
    d_idx

    I_peak = zeros(1,run_time);
    I_peak_time = zeros(1,run_time);
    J_peak = zeros(1,run_time);
    J_peak_time = zeros(1,run_time);

    for run_idx = 1:run_time
    
        m0_1 = all_m0_1(d_idx);
        m0_2 = all_m0_2(d_idx);
        m0_3 = all_m0_3(d_idx);
        m_1 = all_m_1(d_idx);
        m_2 = all_m_2(d_idx);
        m_3 = all_m_3(d_idx);

        G_ba1 = ba_network(N, m0_1,m_1);
        G_ba2 = ba_network(N, m0_2,m_2);
        G_ba3 = ba_network(N, m0_3,m_3);

        A1 = adjacency(G_ba1);
        A2 = adjacency(G_ba2);
        A3 = adjacency(G_ba3);
    
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

        [peak_I, idx_peak_I] = max(nodes_I);
        I_peak(run_idx) = peak_I;
        I_peak_time(run_idx) = time_vals(idx_peak_I);
        
        [peak_J, idx_peak_J] = max(nodes_J);
        J_peak(run_idx) = peak_J;
        J_peak_time(run_idx) = time_vals(idx_peak_J);
    end
    avg_I_peak(d_idx) = mean(I_peak);
    avg_I_peak_time(d_idx) = mean(I_peak_time);
    avg_J_peak(d_idx) = mean(J_peak);
    avg_J_peak_time(d_idx) = mean(J_peak_time);
end

results_table = table(all_k1(:), all_k2(:), all_k3(:), ...
avg_I_peak, avg_I_peak_time, avg_J_peak, avg_J_peak_time, ...
'VariableNames', {'k1','k2','k3','AvgPeak_I','AvgTimePeak_I','AvgPeak_J','AvgTimePeak_J'});

disp(results_table);

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

function G = ba_network(N, m0, m)
    A = zeros(N);
    for i = 1:m0
        for j = i+1:m0
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end

    for newNode = m0+1:N
        degrees = sum(A(1:newNode-1, 1:newNode-1), 2);
        prob = degrees / sum(degrees);
        targets = zeros(1, m);
        selected = false(1, newNode-1);

        k = 0;
        while k < m
            r = rand;
            cum_prob = cumsum(prob);
            idx = find(r <= cum_prob, 1);
            if ~selected(idx)
                k = k + 1;
                targets(k) = idx;
                selected(idx) = true;
            end
        end

        for i = 1:m
            A(newNode, targets(i)) = 1;
            A(targets(i), newNode) = 1;
        end
    end

    G = graph(A);
end