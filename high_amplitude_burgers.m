function high_amplitude_burgers()
    % 高幅值Burgers方程的隐式求解
    % 针对幅值约为1000的情况进行优化
    clear all; close all; clc
    
    % 参数设置
    L = 30.0;           % 空间域长度
    T = 28.0;           % 总时间（减小时间以避免发散）
    Nx = 200;          % 增加空间网格点数
    Nt = 1000;         % 增加时间步数
    Re = 1000;          % 雷诺数
    amplitude = 1054;  % 信号幅值
    
    dx = L / (Nx - 1); % 空间步长
    dt = T / Nt;       % 时间步长
    x = linspace(0, L, Nx); % 空间网格
    
    % 初始条件 - 高幅值信号
    u0 = amplitude * (sin(2*pi/30*x));
    
    % 初始化解矩阵
    u = zeros(Nx, Nt);
    u(:,1) = u0';
    
    % 使用自适应时间步长的牛顿迭代法
    for n = 1:Nt-1
        % 当前时间步的解
        u_current = u(:, n);
        
        % 尝试求解下一个时间步
        max_iter = 10;      % 最大迭代次数
        tol = 1e-6;         % 收敛容差
        u_guess = u_current; % 初始猜测
        
        for iter = 1:max_iter
            % 计算残差和雅可比矩阵
            [F, J] = burgers_jacobian(u_guess, u_current, dx, dt, Re, Nx);
            
            % 求解线性系统
            delta = -J \ F;
            
            % 更新解
            u_next = u_guess + delta;
            
            % 检查收敛
            if norm(delta) < tol
                break;
            end
            
            u_guess = u_next;
            
            % 如果迭代次数过多，减小时间步长并重试
            if iter == max_iter
                warning('在时间步 %d 未收敛，减小时间步长', n);
                dt = dt / 2;
                [F, J] = burgers_jacobian(u_guess, u_current, dx, dt, Re, Nx);
                delta = -J \ F;
                u_next = u_guess + delta;
            end
        end
        
        u(:, n+1) = u_next;
        
        % 检查解是否发散
        if any(isnan(u_next)) || any(isinf(u_next))
            error('解在时间步 %d 发散', n);
        end
    end
    
    % 可视化结果
    plot_high_amplitude_results(x, u, T, Nt, amplitude);
end

function [F, J] = burgers_jacobian(u_next, u_current, dx, dt, Re, Nx)
    % 计算Burgers方程的残差和雅可比矩阵
    
    % 初始化
    F = zeros(Nx, 1);
    J = zeros(Nx, Nx);
    
    % 内部点
    for i = 2:Nx-1
        % 对流项 (使用迎风格式以提高稳定性)
        if u_next(i) >= 0
            conv_term = u_next(i) * (u_next(i) - u_next(i-1)) / dx;
        else
            conv_term = u_next(i) * (u_next(i+1) - u_next(i)) / dx;
        end
        
        % 扩散项
        diff_term = (1/Re) * (u_next(i+1) - 2*u_next(i) + u_next(i-1)) / dx^2;
        
        % 时间导数
        time_deriv = (u_next(i) - u_current(i)) / dt;
        
        % 残差
        F(i) = time_deriv + conv_term - diff_term;
        
        % 雅可比矩阵 (对流项的导数)
        if u_next(i) >= 0
            % 对流项对u_i的导数
            J(i, i) = J(i, i) + (2*u_next(i) - u_next(i-1)) / dx;
            % 对流项对u_{i-1}的导数
            J(i, i-1) = J(i, i-1) - u_next(i) / dx;
        else
            % 对流项对u_i的导数
            J(i, i) = J(i, i) + (u_next(i+1) - 2*u_next(i)) / dx;
            % 对流项对u_{i+1}的导数
            J(i, i+1) = J(i, i+1) + u_next(i) / dx;
        end
        
        % 扩散项的导数
        J(i, i-1) = J(i, i-1) - (1/Re) / dx^2;
        J(i, i) = J(i, i) + 2*(1/Re) / dx^2;
        J(i, i+1) = J(i, i+1) - (1/Re) / dx^2;
        
        % 时间导数的导数
        J(i, i) = J(i, i) + 1/dt;
    end
    
    % 边界条件 (Dirichlet边界条件)
    F(1) = u_next(1);    % u(0,t) = 0
    F(Nx) = u_next(Nx);  % u(L,t) = 0
    
    J(1, :) = 0; J(1, 1) = 1;
    J(Nx, :) = 0; J(Nx, Nx) = 1;
end

function plot_high_amplitude_results(x, u, T, Nt, amplitude)
    % 高幅值结果可视化函数
    
    % 创建时间向量
    t = linspace(0, T, Nt);
    
    % 3D表面图
    figure;
    surf(t, x, u, 'EdgeColor', 'none');
    xlabel('时间 (t)');
    ylabel('空间 (x)');
    zlabel('u(x,t)');
    title(sprintf('高幅值Burgers方程数值解 (幅值=%d)', amplitude));
    colormap('jet');
    colorbar;
    
    % 不同时间点的解
    figure;
    % plot(x, u(:,1), 'r-', 'LineWidth', 2, 'DisplayName', 't=0');
    hold on;
    % plot(x, u(:,round(Nt/4)), 'g-', 'LineWidth', 2, 'DisplayName', ['t=' num2str(T/4)]);
    % plot(x, u(:,round(Nt/2)), 'b-', 'LineWidth', 2, 'DisplayName', ['t=' num2str(T/2)]);
    plot(x, u(:,end), 'k-', 'LineWidth', 2, 'DisplayName', ['t=' num2str(T)]);
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('不同时间点的解 (幅值=%d)', amplitude));
    legend;
    grid on;
    
    % 计算并显示最大最小值
    fprintf('初始最大值: %.2f\n', max(u(:,1)));
    fprintf('最终最大值: %.2f\n', max(u(:,end)));
    fprintf('初始最小值: %.2f\n', min(u(:,1)));
    fprintf('最终最小值: %.2f\n', min(u(:,end)));
end