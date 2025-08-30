function burgers_equation_implicit()
% burger的幅值不能过大
    clear all; close all; clc
    % 参数设置
    L = 20;           % 空间域长度
    T = 28;             % 总时间
    Nx = 100;           % 空间网格点数
    Nt = 500;          % 时间步数
    nu = 0.1;          % 粘性系数
    
    dx = L / (Nx - 1);  % 空间步长
    dt = T / Nt;        % 时间步长
    r = nu * dt / (dx^2); % 扩散系数
    
    % 初始化网格和初始条件
    x = linspace(0, L, Nx);
    % u0 = sin(2*pi*x) + 0.5*sin(pi*x); % 初始条件
    u0 = 2*sin(2*pi/20*x); % 初始条件
    
    % 初始化解矩阵
    u = zeros(Nx, Nt);
    u(:,1) = u0';
    
    % 构建三对角矩阵
    main_diag = (1 + r) * ones(Nx, 1);
    off_diag = (-r/2) * ones(Nx-1, 1);
    
    A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
    
    % 边界条件（Dirichlet边界条件）
    A(1, :) = 0; A(1, 1) = 1;
    A(Nx, :) = 0; A(Nx, Nx) = 1;
    
    % 时间迭代
    for n = 1:Nt-1
        % 当前时间步的解
        u_current = u(:, n);
        
        % 构建右侧向量（考虑非线性项）
        b = zeros(Nx, 1);
        for i = 2:Nx-1
            % Crank-Nicolson格式
            nonlinear_term = -0.25 * (u_current(i+1)^2 - u_current(i-1)^2) / dx;
            diffusion_term = 0.5 * r * (u_current(i+1) - 2*u_current(i) + u_current(i-1));
            
            b(i) = u_current(i) + dt * nonlinear_term + diffusion_term;
        end
        
        % 边界条件
        b(1) = 0;   % u(0,t) = 0
        b(Nx) = 0;  % u(L,t) = 0
        
        % 求解线性系统
        u(:, n+1) = A \ b;
    end
    
    % 可视化结果
    figure;
    surf(linspace(0, T, Nt), x, u, 'EdgeColor', 'none');
    xlabel('时间');
    ylabel('空间');
    zlabel('u(x,t)');
    title('一维Burgers方程数值解（隐式方法）');
    
    % 绘制特定时间点的解
    figure;
    plot(x, u(:,1), 'r-', 'LineWidth', 2, 'DisplayName', 't=0');
    hold on;
    plot(x, u(:,round(Nt/4)), 'g-', 'LineWidth', 2, 'DisplayName', ['t=' num2str(T/4)]);
    plot(x, u(:,round(Nt/2)), 'b-', 'LineWidth', 2, 'DisplayName', ['t=' num2str(T/2)]);
    plot(x, u(:,end), 'k-', 'LineWidth', 2, 'DisplayName', ['t=' num2str(T)]);
    xlabel('x');
    ylabel('u(x,t)');
    title('不同时间点的解');
    legend;
    grid on;
end