rows=40; cols=40; N=rows * cols;
k_1=3; k_2=3; k_3=3;

A1=generate_LA2D(rows,cols,k_1);
A2=generate_LA2D(rows,cols,k_2);
A3=generate_LA2D(rows,cols,k_3);

D=0.005; r=0.1; A=0.1; K=1;
beta_10.3; beta_02=0.15; beta_10=0.1; beta_02=0.1; beta_12=0.05;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.02; alpha_2=0.15; alpha_12=0.1;

d11=0.4; d12=-0.2; d13=-0.2; d22=0.01; d33=4.8;

S_eq=0.99128; I_eq=0.010155; J_eq=0.032663; C_eq=0.010403

L1=A1-diag(sum(A1, 2)); 
L2=A2-diag(sum(A2, 2)); 
L3=A3-diag(sum(A3, 2));

S0=S_eq+1e-5*randn(N,1);
I0=I_eq+1e-5*randn(N,1);
J0=J_eq+1e-5*randn(N,1);
C0=C_eq+1e-5*randn(N,1);

y0 = [S0;I0;J0;C0];

output_folder = ""; % update folder name here
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

num_t=2000; time_bet=1;

amplitude_overall=zeros(1, num_t+1);
time_vals=zeros(1, num_t+1);

for i=1:num_t
    tspan=[0 time_bet];
    
    ode_wrap=@(t, y) coinfect_diff(t,y,L1,L2,L3,r,K,A,D,beta_1,beta_2,beta_10,beta_02,beta_12,gamma_1,gamma_2,alpha_1,alpha_2,alpha_12,d11,d12,d13,d22,d33,N);
    [t, y]=ode113(ode_wrap, tspan, y0);

    S=y(end, 1:N)'; I=y(end, N+1:2*N)'; J=y(end, 2*N+1:3*N)'; C=y(end, 3*N+1:end)';

    S=max(S, 0); I=max(I, 0); J=max(J, 0); C=max(C,0);

    S_plot = reshape(S, [cols, rows])';
    I_plot = reshape(I, [cols, rows])';
    J_plot = reshape(J, [cols, rows])';
    C_plot = reshape(C, [cols, rows])';

    if mod(i, 50) == 0
        cbSize=14;

        figS=figure('Visible','off');
        imagesc(S_plot);
        colormap(figS,'turbo');
        c=colorbar; 
        c.TickLabelInterpreter='latex';
        c.FontSize=cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figS, fullfile(output_folder, sprintf('S_t%03d_40.png', i*time_bet)), 'Resolution', 600);
        close(figS);
        
        figI=figure('Visible','off');
        imagesc(I_plot);
        colormap(figI,'turbo');
        c=colorbar; 
        c.TickLabelInterpreter='latex';
        c.FontSize=cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figI, fullfile(output_folder, sprintf('I_t%03d_40.png', i*time_bet)), 'Resolution', 600);
        close(figI);
        
        figJ=figure('Visible','off');
        imagesc(J_plot);
        colormap(figJ,'turbo');
        c=colorbar; 
        c.TickLabelInterpreter='latex';
        c.FontSize=cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figJ, fullfile(output_folder, sprintf('J_t%03d_40.png', i*time_bet)), 'Resolution', 600);
        close(figJ);

        figC=figure('Visible','off');
        imagesc(C_plot);
        colormap(figC,'turbo');
        c=colorbar; 
        c.TickLabelInterpreter='latex';
        c.FontSize=cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figC, fullfile(output_folder, sprintf('C_t%03d_40.png', i*time_bet)), 'Resolution', 600);
        close(figC);
    end

    amplitude_overall(i+1)=sqrt(sum((S-S_eq).^2+(I-I_eq).^2+(J-J_eq).^2+(C-C_eq).^2));
    time_vals(i+1)=i*time_bet;

    y0 = [S;I;J;C];
end

close all;

function dydt = coinfect_diff(t,y,L1,L2,L3,r, K, A,D,beta_1,beta_2,beta_10,beta_02,beta_12,gamma_1,gamma_2,alpha_1,alpha_2,alpha_12,d11,d12,d13,d22,d33,N)
    S = y(1:N); I = y(N+1:2*N); J = y(2*N+1:3*N); C = y(3*N+1:end);

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


function A = generate_LA2D(rows, cols, range)
    N=rows*cols;
    A=zeros(N);
    [X, Y]=meshgrid(1:cols, 1:rows);
    coords = [X(:), Y(:)];
    for idx=1:N
        xi=coords(idx, 1);
        yi=coords(idx, 2);
        for dx=-range:range
            for dy=-range:range
                if dx==0 && dy==0, continue; end
                if sqrt(dx^2+dy^2) <= range
                    xj=xi+dx; yj=yi+dy;
                    if xj>=1 && xj<=cols && yj>=1 && yj<=rows
                        jdx = (yj-1)*cols+xj;
                        A(idx, jdx)=1;
                    end
                end
            end
        end
    end
end