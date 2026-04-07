% PID Controller Example with Step Perturbation
clear; close all; clc;
% figure;
% System Parameters
Kp = 1; % Proportional gain
Ki = 0; % Integral gain
Kd = 0; % Derivative gain
tau = 1; % Time constant of the system

% Simulation Parameters
dt = 0.01; % Time step
t_end = 100; % End time
time = 0:dt:t_end; % Time vector
n = length(time);

% Desired Set Point
set_point = 1;

% Perturbation Parameters
step_time_1 = 35; % Time at which the step goes from 1 to 0
step_time_2 = 65; % Time at which the step goes from 0 to 1

step_time_1_C = 50; % Time at which the step goes from 1 to 0
step_time_2_C = 70; % Time at which the step goes from 0 to 1

% PID Controller Initialization
err = 0;
err_sum = 0;
err_prev = 0;
u = 0;

% System Initialization
y = 0;
y_prev = 0;

% Simulation Loop
for i = 1:n
    % Calculate error
    err = set_point - y;
    err_sum = err_sum + err*dt;
    err_diff = (err - err_prev)/dt;

    % PID Controller
    u = Kp*err + Ki*err_sum + Kd*err_diff;

    % Perturbation (Step Disturbance)
    if time(i) < step_time_1
        N = 1;
    elseif time(i) < step_time_2
        N = 2;
    else
        N = 0.5;
    end
    
    if time(i) < step_time_1_C
        C = 1;
    elseif time(i) < step_time_2_C
        C = 2;
    else
        C = 1;
    end
    
    

    % System Dynamics (First-Order System with Disturbance)
%     dydt = (- y*C/N  + u ) / tau;
    dydt = (- y  + N +u) / tau;
    y = y_prev + dydt*dt;
    
    
    % Update variables for the next iteration
    err_prev = err;
    y_prev = y;

    % Save data for plotting
    y_data(i) = y;
    out_data(i) = 5-y;
    u_data(i) = u;
    disturbance_data(i) = N;
    disturbance_C_data(i) = C;
end

fig=figure('Renderer', 'painters', 'Position', [10 20 1000 1000]);

pos1 = [0.1, 0.4+0.4, 0.3, 0.1];
pos2 = [0.1, 0.32+0.4, 0.3, 0.04];

subplot1 = axes('Parent', fig, 'Position', pos1);
subplot2 = axes('Parent', fig, 'Position', pos2);

% Plotting the results

plot(subplot1,time, out_data,time, 5-ones(size(time))*set_point,'k--', 'LineWidth', 2);
xticks(subplot1,[])
yticks(subplot1,[])
box(subplot1,'off');
set(subplot1,'FontName','Times new roman','FontSize',15,'linewidth',1)
% plot(subplot2,time, u_data, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Control Signal (u)');
% title('PID Control Signal');
% grid on;
% 

plot(subplot2,time, disturbance_data, 'LineWidth', 2);
xlabel('Time');
xticks(subplot2,[])
yticks(subplot2,[])
box(subplot2,'off');
set(subplot2,'FontName','Times new roman','FontSize',15,'linewidth',1)
% ylabel('Disturbance');
% title('Nitrogen Concentration');
% grid on;
% 
% subplot(4,1,4);
% plot(time, disturbance_C_data, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Disturbance');
% title('Carbon Concentration');
% grid on;


%%
Kp = 0; % Proportional gain
Ki = 0.5; % Integral gain
Kd = 0; % Derivative gain
tau = 1; % Time constant of the system

% Simulation Parameters
dt = 0.01; % Time step
t_end = 100; % End time
time = 0:dt:t_end; % Time vector
n = length(time);

% Desired Set Point
set_point = 1;

% Perturbation Parameters
step_time_1 = 35; % Time at which the step goes from 1 to 0
step_time_2 = 65; % Time at which the step goes from 0 to 1

step_time_1_C = 50; % Time at which the step goes from 1 to 0
step_time_2_C = 70; % Time at which the step goes from 0 to 1

% PID Controller Initialization
err = 0;
err_sum = 0;
err_prev = 0;
u = 0;

% System Initialization
y = 0;
y_prev = 0;

% Simulation Loop
for i = 1:n
    % Calculate error
    err = set_point - y;
    err_sum = err_sum + err*dt;
    err_diff = (err - err_prev)/dt;

    % PID Controller
    u = Kp*err + Ki*err_sum + Kd*err_diff;

    % Perturbation (Step Disturbance)
    if time(i) < step_time_1
        N = 1;
    elseif time(i) < step_time_2
        N = 2;
    else
        N = 1;
    end
    
    if time(i) < step_time_1_C
        C = 1;
    elseif time(i) < step_time_2_C
        C = 2;
    else
        C = 1;
    end
    
    

    % System Dynamics (First-Order System with Disturbance)
%     dydt = (- y*C/N  + u ) / tau;
    dydt = (- y  + N + u) / tau;
    y = y_prev + dydt*dt;
    
    
    % Update variables for the next iteration
    err_prev = err;
    y_prev = y;

    % Save data for plotting
    y_data(i) = y;
    out_data(i) = 5-y;
    u_data(i) = u;
    disturbance_data(i) = N;
    disturbance_C_data(i) = C;
end

% fig=figure('Renderer', 'painters', 'Position', [10 20 1000 1000]);
pos3 = [0.1, 0.4+0.1, 0.3, 0.1];
pos4 = [0.1, 0.32+0.1, 0.3, 0.04];

subplot1 = axes('Parent', fig, 'Position', pos3);
subplot2 = axes('Parent', fig, 'Position', pos4);


plot(subplot1,time, out_data, time, 5-ones(size(time))*set_point,'k--','LineWidth', 2);
xticks(subplot1,[])
yticks(subplot1,[])
box(subplot1,'off');
set(subplot1,'FontName','Times new roman','FontSize',15,'linewidth',1)



plot(subplot2,time, disturbance_data, 'LineWidth', 2);
xlabel('Time');
xticks(subplot2,[])
yticks(subplot2,[])
box(subplot2,'off');
set(subplot2,'FontName','Times new roman','FontSize',15,'linewidth',1)



%%
Kp = 0; % Proportional gain
Ki = 0.5; % Integral gain
Kd = 0; % Derivative gain
tau = 1; % Time constant of the system

% Simulation Parameters
dt = 0.01; % Time step
t_end = 100; % End time
time = 0:dt:t_end; % Time vector
n = length(time);

% Desired Set Point
set_point = 1;

% Perturbation Parameters
step_time_1 = 35; % Time at which the step goes from 1 to 0
step_time_2 = 65; % Time at which the step goes from 0 to 1

step_time_1_C = 50; % Time at which the step goes from 1 to 0
step_time_2_C = 70; % Time at which the step goes from 0 to 1

% PID Controller Initialization
err = 0;
err_sum = 0;
err_prev = 0;
u = 0;

% System Initialization
y = 0;
y_prev = 0;

% Simulation Loop
for i = 1:n
    % Calculate error
    err = set_point - y;
    err_sum = err_sum + err*dt;
    err_diff = (err - err_prev)/dt;

    % PID Controller
    u = Kp*err + Ki*err_sum + Kd*err_diff;

    % Perturbation (Step Disturbance)
    if time(i) < step_time_1
        N = 1;
    elseif time(i) < step_time_2
        N = 2;
    else
        N = 1;
    end
    
    if time(i) < step_time_1_C
        C = 1;
    elseif time(i) < step_time_2_C
        C = 2;
    else
        C = 1;
    end
    
    

    % System Dynamics (First-Order System with Disturbance)
%     dydt = (- y*C/N  + u ) / tau;
    dydt = (- y/N + u) / tau;
    y = y_prev + dydt*dt;
    
    
    % Update variables for the next iteration
    err_prev = err;
    y_prev = y;

    % Save data for plotting
    y_data(i) = y;
    out_data(i) = 5-y;
    u_data(i) = u;
    disturbance_data(i) = N;
    disturbance_C_data(i) = C;
end

pos3 = [0.1, 0.4-0.2, 0.3, 0.1];
pos4 = [0.1, 0.32-0.2, 0.3, 0.04];

subplot1 = axes('Parent', fig, 'Position', pos3);
subplot2 = axes('Parent', fig, 'Position', pos4);



plot(subplot1,time, out_data,time, 5-ones(size(time))*set_point,'k--','LineWidth', 2);
xticks(subplot1,[])
yticks(subplot1,[])
box(subplot1,'off');
set(subplot1,'FontName','Times new roman','FontSize',15,'linewidth',1)



plot(subplot2,time, disturbance_data, 'LineWidth', 2);
xlabel('Time');
xticks(subplot2,[])
yticks(subplot2,[])
box(subplot2,'off');
set(subplot2,'FontName','Times new roman','FontSize',15,'linewidth',1)


