% ws network -- superinfect

output_folder = ""; % fill in folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

N=1600;

D=0.005; r=0.1; A=0.1; K=1;
beta_1 = 0.3; beta_2 = 0.4; beta_10 = 0.2; beta_02 = 0.3; beta_12 = 0.05;
gamma_1 = 0.1; gamma_2 = 0.05;
alpha_1 = 0.05; alpha_2 = 0.15; alpha_12=0.25;

d11 = 0.3; d12 = 0.1; d13 = 0.1; d22 = 0.3; d33 = 0.1;

all_m0_1 = [5,13,13,13,13,25,25,25];
all_m_1 = [2,6,6,6,6,12,12,12];
all_m0_2 = [5,13,5,5,13,13,5,25];
all_m_2 = [2,6,2,2,6,6,2,12];
all_m0_3 = [5,5,5,13,13,5,13,25];
all_m_3 = [2,2,2,6,6,2,6,12];
colors = lines(length(all_m0_1));

num_t =500;
time_bet = 1;
run_time = 1000;

avg_I_peak = zeros(length(all_m0_1),1);
avg_I_peak_time = zeros(length(all_m0_1),1);
avg_J_peak = zeros(length(all_m0_1),1);
avg_J_peak_time = zeros(length(all_m0_1),1);
avg_C_peak = zeros(length(all_m0_1),1);
avg_C_peak_time = zeros(length(all_m0_1),1);

for d_idx = 1:length(all_m0_1)
    d_idx

    I_peak = zeros(1,run_time);
    I_peak_time = zeros(1,run_time);
    J_peak = zeros(1,run_time);
    J_peak_time = zeros(1,run_time);
    C_peak = zeros(1,run_time);
    C_peak_time = zeros(1,run_time);

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
        I0 = zeros(N, 1); J0 = zeros(N, 1); C0 = zeros(N,1);
        I0(1)=0.05; J0(1600)=0.05;
    
        y0 = [S0;I0;J0;C0];
    
        nodes_I = zeros(1, num_t+1);
        nodes_J = zeros(1, num_t+1);
        nodes_C = zeros(1, num_t+1);
        nodes_I_mono = zeros(1, num_t+1);
        nodes_J_mono = zeros(1, num_t+1);
        time_vals = zeros(1, num_t+1);
    
        nodes_I(1) = sum(I0 > 0.01) / N;
        nodes_J(1) = sum(J0 > 0.01) / N;
        nodes_C(1)=sum(C0>0.005)/N;
        time_vals(1) = 0;
    
        for i = 1:num_t
     
            tspan = [0 time_bet];
            ode_wrap = @(t, y) coinfect_diff(t,y,L1,L2,L3,r, K, A,D,beta_1,beta_2,beta_10,beta_02,beta_12,...
                gamma_1,gamma_2,alpha_1,alpha_2,alpha_12,d11,d12,d13,d22,d33,N);
            [t, y] = ode113(ode_wrap, tspan, y0);

            S = y(end, 1:N)';
            I = y(end, N+1:2*N)';
            J = y(end, 2*N+1:3*N)';
            C = y(end, 3*N+1:end)';
        
            S = max(S, 0);
            I = max(I, 0);
            J = max(J, 0);
            C = max(C,0);
    
            I_mono = I-C;
            J_mono = J-C;
    
            nodes_I(i+1) = sum(I > 0.01) / N;
            nodes_J(i+1) = sum(J > 0.01) / N;
            nodes_C(i+1) = sum(C > 0.005) / N;
    
            nodes_I_mono(i+1) = sum(I_mono > 0.01) / N;
            nodes_J_mono(i+1) = sum(J_mono > 0.01) / N;
    
            time_vals(i+1) = i * time_bet;
    
            y0 = [S;I;J;C];
        end

        [peak_I, idx_peak_I] = max(nodes_I_mono);
        I_peak(run_idx) = peak_I;
        I_peak_time(run_idx) = time_vals(idx_peak_I);
        
        [peak_J, idx_peak_J] = max(nodes_J_mono);
        J_peak(run_idx) = peak_J;
        J_peak_time(run_idx) = time_vals(idx_peak_J);

        [peak_C, idx_peak_C] = max(nodes_C);
        C_peak(run_idx) = peak_C;
        C_peak_time(run_idx) = time_vals(idx_peak_C);
    end
    avg_I_peak(d_idx) = mean(I_peak);
    avg_I_peak_time(d_idx) = mean(I_peak_time);
    avg_J_peak(d_idx) = mean(J_peak);
    avg_J_peak_time(d_idx) = mean(J_peak_time);
    avg_C_peak(d_idx) = mean(C_peak);
    avg_C_peak_time(d_idx) = mean(C_peak_time);
end

results_table = table(all_m0_1(:), all_m0_2(:), all_m0_3(:), ...
avg_I_peak, avg_I_peak_time, avg_J_peak, avg_J_peak_time, avg_C_peak, avg_C_peak_time,...
'VariableNames', {'k1','k2','k3','AvgPeak_I','AvgTimePeak_I','AvgPeak_J','AvgTimePeak_J', 'AvgPeak_C','AvgTimePeak_C'});

disp(results_table);

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
