rows=40; cols=40; N=rows*cols;
k_1=2; k_2=2; k_3=2;

A1=generate_LA2D(rows,cols,k_1);
A2=generate_LA2D(rows,cols,k_2);
A3=generate_LA2D(rows,cols,k_3);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2= 0.15;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.02; alpha_2= 0.15;
sigma=3;

d11= 0.1; d12= -0.2; d13=-0.2; d22= =0.01; d33=4.8;

S_eq=0.82354; I_eq =1.3971; J_eq=0.42647

L1=A1-diag(sum(A1,2)); 
L2=A2-diag(sum(A2,2)); 
L3=A3-diag(sum(A3,2));

S0=S_eq+1e-3*randn(N,1);
I0=I_eq+1e-3*randn(N,1);
J0=J_eq+1e-3*randn(N,1);
y0=[S0;I0;J0];

output_folder = ""; % enter folder name here
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

num_t=4000; time_bet=1;

amplitude_overall=zeros(1, num_t+1);
time_vals=zeros(1, num_t+1);

for i=1:num_t
    tspan=[0 time_bet];
    
    ode_wrap=@(t, y) superinfect_diff(t, y, L1, L2, L3, r,K,A,D,beta_1, beta_2, gamma_1, gamma_2, alpha_1, alpha_2, sigma, d11, d12, d13, d22, d33, N);
    [t, y]=ode113(ode_wrap, tspan, y0);

    S=y(end, 1:N)'; I=y(end, N+1:2*N)'; J=y(end, 2*N+1:end)';

    S=max(S,0); I=max(I,0); J=max(J,0);

    S_plot=reshape(S, [cols, rows])';
    I_plot=reshape(I, [cols, rows])';
    J_plot=reshape(J, [cols, rows])';

    if mod(i,50) == 0
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
        
        figJ = figure('Visible','off');
        imagesc(J_plot);
        colormap(figJ,'turbo');
        c=colorbar; 
        c.TickLabelInterpreter='latex';
        c.FontSize=cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figJ, fullfile(output_folder, sprintf('J_t%03d_40.png', i*time_bet)), 'Resolution', 600);
        close(figJ);
    end

    amplitude_overall(i+1)=sqrt(sum((S-S_eq).^2+(I-I_eq).^2+(J-J_eq).^2));
    time_vals(i+1)=i*time_bet;

    y0=[S;I;J];
end

close all;

figAmpOverall=figure('Visible', 'on');
clf; hold on;

plot(time_vals, amplitude_overall, '-', 'LineWidth', 1);
xlabel('Time');
ylabel('Amplitude');
grid on;

savefig(figAmpOverall, fullfile(output_folder, 'amp_overall'));

close(figAmpOverall);

function dydt = superinfect_diff(t,y,L1,L2,L3,r,K,A,D,beta_1,beta_2,gamma_1,gamma_2,alpha_1,alpha_2,sigma,d11,d12,d13,d22,d33,N)
    S=y(1:N); I=y(N+1:2*N); J=y(2*N+1:end);

    f = r*S.*(1-S/K).*(S/A-1)-(beta_1*I+beta_2*J).*S./(S+I+J)+gamma_1*I+gamma_2*J-D*S;
    g = I.*(beta_1*S./(S+I+J)-gamma_1-sigma*beta_2*J./(S+I+J))-D*I-alpha_1*I;
    h = J.*(beta_2*S./(S+I+J)-gamma_2+sigma*beta_2*I./(S+I+J))-D*J-alpha_2*J;

    dSdt = f+d11*(L1*S)+d12*(L2*I)+d13*(L3*J);
    dIdt = g+d22*(L2*I);
    dJdt = h+d33*(L3*J);

    dydt=[dSdt;dIdt;dJdt];
end


function A = generate_LA2D(rows, cols, range)
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
                if sqrt(dx^2+dy^2) <= range
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