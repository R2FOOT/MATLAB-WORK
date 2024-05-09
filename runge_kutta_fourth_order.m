function runge_kutta_fourth_order()
    
    fig = figure('Position', [300 200 700 500], 'Name', 'Runge-Kutta Solver');
    
    uicontrol('Style', 'text', 'String', 'Enter ODE (in terms of x & y):', ...
        'Position', [0 470 200 20]);
    ODE = uicontrol('style', 'edit', 'string', '', 'position', [28 445 350 30]);
    uicontrol('Style', 'text', 'String', 'y(x0)','Position', [385, 470, 50, 20]);
    y0_ui = uicontrol('Style', 'edit', 'String', '','Position', [385, 445, 50, 30]);
    uicontrol('Style', 'text', 'String', 'x0','Position', [450, 470, 50, 20]);
    x0_ui = uicontrol('Style', 'edit', 'String', '','Position', [450, 445, 50, 30]);
    uicontrol('Style', 'text', 'String', 'xn','Position', [515, 470, 50, 20]);
    xn_ui = uicontrol('Style', 'edit', 'String', '','Position', [515, 445, 50, 30]);
    uicontrol('Style', 'text', 'String', 'h','Position', [580, 470, 50, 20]);
    step_size = uicontrol('Style', 'edit', 'String', '','Position', [580, 445, 50, 30]);
    push_btn = uicontrol('style', 'pushbutton', 'string', 'SUBMIT','Position', ...
        [100, 370, 100, 30], 'Callback', @compute);
    ax = axes('Parent', fig, 'Position', [0.1, 0.1, 0.8, 0.7]);
    
    function compute(~,~)
        ODE_string = get(ODE, 'String');
        f = str2func(['@(x, y) ', ODE_string]);
        y0 = str2double(get(y0_ui, 'String'));
        h = str2double(get(step_size, 'String'));
        x0 = str2double(get(x0_ui, 'String'));
        xn = str2double(get(xn_ui, 'String'));
        
        x_range = [x0 xn+2];
        
        n_steps = ceil((x_range(2) - x_range(1))/h);
    
        x_values = zeros(1, n_steps+1);
        y_values_RK4 = zeros(1, n_steps+1);
        
        x_values(1) = x0;
        y_values_RK4(1) = y0;
        
        for i = 1:n_steps
            k1 = h * f(x_values(i), y_values_RK4(i));
            k2 = h * f(x_values(i) + h/2, y_values_RK4(i) + k1/2);
            k3 = h * f(x_values(i) + h/2, y_values_RK4(i) + k2/2);
            k4 = h * f(x_values(i) + h, y_values_RK4(i) + k3);

            y_values_RK4(i+1) = y_values_RK4(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
            x_values(i+1) = x_values(i) + h;
        end 
        %syms y(x);
        %exact_solution_sys = dsolve(@(x,y) diff(y) == x-y, y(x0) == y0);
        %exact_solution = matlabFunction(exact_solution_sys);
        %y_values_exact = exact_solution(x_values);
        
        plot(ax, x_values, y_values_RK4, '-o', 'LineWidth', 1.5, 'MarkerSize',...
            6, 'DisplayName', 'Runge-Kutta 4th order');
        hold(ax, 'off');
        %plot(ax, x_values, y_values_exact, '-s', 'LineWidth', 1.5, 'MarkerSize',...
        %    6, 'DisplayName', 'Exact Solution');
        %hold(ax, 'off');
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        title(ax, 'Comparison of Runge-Kutta 4th order and Exact Solution');
        legend(ax, 'Location', 'best');
        grid(ax, 'on');
    end
    
end
