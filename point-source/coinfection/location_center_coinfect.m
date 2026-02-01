% i starts at the center, j varies -- coinfection

output_folder=""; % fill in folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows=40; cols=40; N=rows*cols;

k_1=1; k_2=1; k_3=1;

A1=generate_LA2D(rows, cols, k_1);
A2=generate_LA2D(rows, cols, k_2);
A3=generate_LA2D(rows, cols, k_3);

L1=A1-diag(sum(A1, 2)); 
L2=A2-diag(sum(A2, 2)); 
L3=A3-diag(sum(A3, 2));

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.4; beta_10=0.2; beta_02=0.3;
gamma_1=0.1; gamma_2=0.05;
alpha_1=0.05; alpha_2=0.15;
d11=0.3; d12=0.1; d13=0.1; d22=0.3; d33=0.1;


all_dist = [1190,9*40+10,206,1600];
% 19, 21, 29, 39

G_ws1=graph(A1);

disp(all_dist);

colors = lines(length(all_dist));

figI=figure('Visible', 'off'); hold on;
figJ=figure('Visible', 'off'); hold on;
figC=figure('Visible', 'off'); hold on;

figI_mono=figure('Visible', 'off'); hold on;
figJ_mono=figure('Visible', 'off'); hold on;

num_t=250;
time_bet=1;

for d_idx=1:length(all_dist)
    
    d=all_dist(d_idx);

    S0=ones(N, 1); I0=zeros(N, 1); J0=zeros(N, 1); C0=zeros(N,1);

    I0_mono=zeros(N,1); J0_mono=zeros(N,1);
    I0(820)=0.05; J0(d)=0.05;
    I0_mono(820)=0.05; J0_mono(d)=0.05;

    y0=[S0;I0;J0;C0];

    nodes_I=zeros(1, num_t+1);
    nodes_I_mono=zeros(1, num_t+1);
    nodes_J=zeros(1, num_t+1);
    nodes_J_mono=zeros(1, num_t+1);
    nodes_C=zeros(1, num_t+1);
    time_vals=zeros(1, num_t+1);

    nodes_I(1)=sum(I0>0.01)/N;
    nodes_I_mono(1)=sum(I0_mono>0.01)/N;
    nodes_J(1)=sum(J0>0.01)/N;
    nodes_J_mono(1)=sum(J0_mono>0.01)/N;
    nodes_C(1)=sum(C0>0.005)/N;
    time_vals(1)=0;

    for i=1:num_t
 
        tspan=[0 time_bet];
        ode_wrap=@(t, y) coinfect_diff(t,y,L1,L2,L3,r,K,A,D,beta_1,beta_2,beta_10,beta_02,beta_12,gamma_1,gamma_2,alpha_1,alpha_2,alpha_12,d11,d12,d13,d22,d33,N);
        [t, y]=ode113(ode_wrap, tspan, y0);

        S=y(end, 1:N)'; I=y(end, N+1:2*N)'; J=y(end, 2*N+1:3*N)'; C=y(end, 3*N+1:end)';
    
        S=max(S, 0); I=max(I, 0); J=max(J, 0); C=max(C,0);

        I_mono=I-C; J_mono=J-C;

        nodes_I(i+1)=sum(I>0.01)/N;
        nodes_J(i+1)=sum(J>0.01)/N;
        nodes_C(i+1)=sum(C>0.005)/N;

        nodes_I_mono(i+1)=sum(I_mono>0.01)/N;
        nodes_J_mono(i+1)=sum(J_mono>0.01)/N;

        time_vals(i+1)=i*time_bet;

        y0=[S;I;J;C];
    end

    figure(figI); 
    plot(time_vals, nodes_I, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figJ); 
    plot(time_vals, nodes_J, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figC); 
    plot(time_vals, nodes_C, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figI_mono); 
    plot(time_vals, nodes_I_mono, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figJ_mono); 
    plot(time_vals, nodes_J_mono, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);
end

figure(figI); xlabel('Time'); ylabel('I-Spread Index');
legend(compose('d_{ij} = %.3f', all_dist), 'Location', 'best'); grid on;
saveas(figI, fullfile(output_folder, 'location_I_center'));
close(figI);

figure(figJ); xlabel('Time'); ylabel('J-Spread Index');
legend(compose('d_{ij} = %.3f', all_dist), 'Location', 'best'); grid on;
saveas(figJ, fullfile(output_folder, 'location_J_center'));
close(figJ);

figure(figC); xlabel('Time'); ylabel('C-Spread Index');
legend(compose('d_{ij} = %.3f', all_dist), 'Location', 'best'); grid on;
saveas(figC, fullfile(output_folder, 'location_C_center'));
close(figC);

figure(figI_mono); xlabel('Time'); ylabel('I_1-Spread Index');
legend(compose('d_{ij} = %.3f', all_dist), 'Location', 'best'); grid on;
saveas(figI_mono, fullfile(output_folder, 'location_I_mono_center'));
close(figI_mono);

figure(figJ_mono); xlabel('Time'); ylabel('I_2-Spread Index');
legend(compose('d_{ij} = %.3f', all_dist), 'Location', 'best'); grid on;
saveas(figJ_mono, fullfile(output_folder, 'location_J_mono_center'));
close(figJ_mono);

function dydt=coinfect_diff(t,y,L1,L2,L3,r,K,A,D,beta_1,beta_2,beta_10,beta_02,beta_12,gamma_1,gamma_2,alpha_1,alpha_2,alpha_12,d11,d12,d13,d22,d33,N)
    S=y(1:N); I=y(N+1:2*N); J=y(2*N+1:3*N); C=y(3*N+1:end);

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

    dydt = [dSdt;dIdt;dJdt;dCdt];
end

function A=generate_LA2D(rows, cols, range)
    N=rows*cols;
    A=zeros(N);
    [X, Y]=meshgrid(1:cols, 1:rows);
    coords=[X(:), Y(:)];
    for idx=1:N
        xi=coords(idx, 1);
        yi=coords(idx, 2);
        for dx=-range:range
            for dy=-range:range
                if dx==0 && dy==0, continue; end
                if sqrt(dx^2+dy^2)<=range
                    xj=xi+dx; yj=yi+dy;
                    if xj>=1 && xj<=cols && yj>=1 && yj<=rows
                        jdx=(yj-1)*cols+xj;
                        A(idx, jdx)=1;
                    end
                end
            end
        end
    end

end
