% File: MAE162A.m

classdef MAE162A
     methods (Static)
    % Solving Vector Loops
    function sol = NRsolve(f1, f2, angles, guess, depth)
        J = jacobian([f1 f2], sym(angles));
        F = -inv(J)*[f1; f2];
        sol.angles(1) = 0;
        sol.angles(2) = 0;

        for i = 1:depth
            guess = (F(guess(1), guess(2)) + guess);
        end
        sol.(angles(1)) = mod(guess(1), 2*pi);
        sol.(angles(2)) = mod(guess(2), 2*pi);
    end

    function [s_out, exists] = solveLoop(f1, f2, s_in,optimize)
        s_out = s_in;
        x = s_out.parameters(1);
        y = s_out.parameters(2);
        s_out.(x) = 0;
        s_out.(y) = 0;
        exists = false;
        guess = [0 1; 0 5];
        NRguess = [0 0];
        eqns = [f1 == 0, f2 == 0];

        if x == "a3" && y == "a4"
            guess = [2, 5.5; 0 2*pi];
            NRguess = [2, 0];
        end

        if optimize
            sol = vpasolve(eqns, sym([x y]), guess, "Random", true);
            if ~isempty(sol.(x)) && ~isempty(sol.(y))
                s_out.(x) = double(sol.(x)(1));
                s_out.(y) = double(sol.(y)(1));
                exists = true;
            end
        else
            sol = MAE162A.NRsolve(f1, f2, [x y], NRguess, 4);
            if (sol.(x) ~= 0) && (sol.(y) ~= 0)
                s_out.(x) = sol.(x)(1);
                s_out.(y) = sol.(y)(1);
                exists = true;
            end
        end
    end

    function updated = velocityLoop(data, links)
        a = [string(['a' num2str(links(1))]), string(['a' num2str(links(2))])];
        r = [string(['r' num2str(links(1))]), string(['r' num2str(links(2))])];
        
        J = zeros(2, 2);
        f = zeros(2, 1);
        
        for i = 1:length(data)
            s = data(i);
            f(1) = s.r2*sin(s.a2);
            f(2) = -s.r2*cos(s.a2);
            J(1,1) = -s.(a(1))*sin(s.(a(1)));
            J(1,2) = -s.(a(2))*cos(s.(a(2)));
            J(2,1) = s.(a(1))*sin(s.(a(1)));
            J(2,2) = s.(a(2))*cos(s.(a(2)));
            kc = J\f;
            
            s.(['kc1_', num2str(links(1))])= kc(1);
            s.(['kc1_', num2str(links(2))])= kc(2);
            updated(i) = s;
        end

    end
    
    % Building Vector Loops
    function [f1, f2] = buildLoop(s, vectors)
        x = sym(s.parameters(1));
        y = sym(s.parameters(2));
        f1(x, y) = sym(0);
        f2(x, y) = sym(0);
        for i = 1:length(vectors)
            fields = s.variables(i, :);
            R = vectors{i}(s.(fields(1)), s.(fields(2)));
            f1 = f1 + R(1);
            f2 = f2 + R(2);
        end
    end

    % Parsing Vector Loops
    function vectors = parseLoop(vars)
        vectors = cell(1, length(vars));
        for i = 1:length(vars)
            v = sym(vars(i, :));
            R(v(1), v(2)) = [v(1)*cos(v(2)); v(1)*sin(v(2))];
            vectors{i} = R;
        end
    end
    
    % Showing Vector Loops
    function showLoop(s, varargin)
        n = length(s.variables);
       [x, y] = MAE162A.getPositions(s);

        hold on
        q = quiver(0, 0, x(1), y(1), "off");
        q.MaxHeadSize = 0.5;
        if ~isempty(varargin)
            q.Color = varargin{1};
        end
        
        for i = 2:n  
            q = quiver(sum(x(1:i-1)), sum(y(1:i-1)), x(i), y(i), "off");
            q.MaxHeadSize = 0.5;
            if ~isempty(varargin)
                q.Color = varargin{1};
            end
        end
        hold off
    end
    
    function showSimulation(fig, scale, movie, name, varargin)
        if movie
            v = VideoWriter([name(find(~isspace(name))), '.avi']);
            open(v);
        end

        if length(varargin) == 1
            data = varargin{1};
            for i = 1:length(data)
                clf(fig);
                axis(scale);
                title(gca, name);

                MAE162A.showLoop(data(i));
                MAE162A.showBox();
                drawnow

                if movie
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end
            end
        end

        if length(varargin) == 2
            data1 = varargin{1};
            data2 = varargin{2};
            for i = 1:length(data1)
                clf(fig);
                axis(scale);
                title(gca, name);

                MAE162A.showLoop(data1(i), "red");
                MAE162A.showLoop(data2(i), "blue");
                MAE162A.showTop(data2(i))
                MAE162A.showBox();
                drawnow

                if movie
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end
            end
        end

        if movie
            close(v);
        end
    end

    function showTop(s)
        [x, y] = MAE162A.getPositions(s);
        x1 = sum(x(1:1));
        x2 = sum(x(1:2));
        y1 = sum(y(1:1));
        y2 = sum(y(1:2));
        
        h = 200;

        hold on
        plot([x1 x1 x2 x2],[y1 y1+h y2+h y2])
        plot(0.5*(x1+x2), 0.5*(y1+y2)+h, 'o')
        hold off
    end

    function showBox()
        hold on
        line([0, 850], [1050, 1050],'Color','k')
        line([0, 850], [0, 0],'Color','k')
        line([0 0],[0 1050],'Color','k');
        line([850 850],[0 1050],'Color','k');
        hold off
    end

    % Initializing Data
    function s = initData(variables, values, parameters)
        for i = 1:length(variables)
            vals = values(i, :);
            vars = variables(i, :);
            s.(vars(1)) = vals(1);
            s.(vars(2)) = vals(2);
        end
        x = parameters(1);
        y = parameters(2);
        s.(x) = sym(x);
        s.(y) = sym(y);
        s.parameters = parameters;
        s.variables = variables;
    end

    % Simulating Vector Loops
    function [data, idx] = simulate1(input, s, optimize)
        n = length(input);
        R = MAE162A.parseLoop(s.variables); 
        data(1, n) = s; % pre-allocate struct array
        idx(1, n) = false; 

        for i = 1:length(input)
            s.a2 = input(i);      
            [data(i), idx(i)] = MAE162A.step(s, R, optimize);
        end    
    end

    function [data2, idx2] = simulate2(data1, s2, optimize)
        n = length(data1);
        data2(1, n) = s2; 
        idx2(1, n) = false;
        
        R = MAE162A.parseLoop(s2.variables); 

        for i = 1:n     
            s = data1(i);

            % finding r16 and a16 (constrained by R2 and R3)
            x = s.r2*cos(s.a2)+0.5*s.r3*cos(s.a3); 
            y = s.r2*sin(s.a2)+0.5*s.r3*sin(s.a3);
            
            s2.r16 = double(sqrt(x^2+y^2));
            s2.a16 = double(atan2(y, x));
            s2.a15 = s.a4;
            
            [data2(i), idx2(i)] = MAE162A.step(s2, R, optimize);
        end
    end

    function [s_out, exists] = step(s_in, R, optimize)
        [f1, f2] = MAE162A.buildLoop(s_in, R);
        [s_out, exists] = MAE162A.solveLoop(f1, f2, s_in, optimize);
    end
    
    % Analysis
    function [x, y] = getPositions(s)
        n = length(s.variables);
        x = zeros(1, n);
        y = zeros(1, n);
        for i = 1:n
            r = s.variables(i, 1);
            a = s.variables(i, 2);
            x(i) =  s.(r) * cos(s.(a));
            y(i) =  s.(r) * sin(s.(a));
        end
    end

    function positionAnalysis(a, maxima, minima, fignum1, fignum2, movie)
        theta2_max = deg2rad(maxima);
        theta2_min = deg2rad(minima);
        theta2 = linspace(theta2_min, theta2_max, 1000)';
        R1 = 4*a;
        R2 = 5*a;
        R3 = 2*a; 
        R4 = 5*a; 
        R6 = 2*a; 
        R5 = 1.25001*R6;
        R15 = 0.5*R4;
        R26 = 0.5*R3;
        
        theta3 = ones(length(theta2),1);
        theta4 = ones(length(theta2),1);
        theta5 = ones(length(theta2),1);
        theta6 = ones(length(theta2),1);
        theta3 = theta3 .* (80 * pi/180);
        theta4 = theta4 .* (300 * pi/180);
        theta5 = theta5 .* (200 * pi/180);
        theta6 = theta6 .* (5 * pi/180);
    
        for j = 1:length(theta2)
            n = 0;
            D = [100; 100];
            theta3_current = theta3(j);
            theta4_current = theta4(j);
            while norm(D,1) > 10^-6 
                n = n + 1;
                func1 = R2*cos(theta2(j)) + R3*cos(theta3_current) + R4*cos(theta4_current) - R1;
                func2 = R2*sin(theta2(j)) + R3*sin(theta3_current) + R4*sin(theta4_current);
                F = [func1; func2];
                J11 = -R3*sin(theta3_current);
                J12 = -R4*sin(theta4_current);
                J21 = R3*cos(theta3_current);
                J22 = R4*cos(theta4_current);
                Jacob1 = [J11 J12; J21 J22];
                D = -Jacob1\F;
                dt3 = D(1);
                dt4 = D(2);
                theta3_current = theta3_current + dt3;
                theta4_current = theta4_current + dt4;
                if n > 100
                    break;
                end
            end
          
            theta3(j) = theta3_current;
            theta4(j) = theta4_current;
             if j < length(theta2)
                 theta3(j+1) = theta3(j);
                theta4(j+1) = theta4(j);
             end
        end
        
        for j = 1: length(theta2)
             n = 0;
             X = [100; 100];
             theta5_current = theta5(j);
             theta6_current = theta6(j);
             while norm(X,1) > 10^(-6)
                   n = n + 1;
                   func3 = R2*cos(theta2(j)) + R26*cos(theta3(j)) + R6*cos(theta6_current) + R5*cos(theta5_current) + R15*cos(theta4(j)) - R1;
                   func4 = R2*sin(theta2(j)) + R26*sin(theta3(j)) + R6*sin(theta6_current) + R5*sin(theta5_current) + R15*sin(theta4(j));
                   f = [func3;func4];
                   J11 = -R5*sin(theta5_current);
                   J12 = -R6*sin(theta6_current);
                   J21 = R5*cos(theta5_current);
                   J22 = R6*cos(theta6_current);
                   Jacob2 =  [J11 J12; J21 J22];
                   X = -Jacob2\f;
                   theta5_current = theta5_current + X(1);
                   theta6_current = theta6_current + X(2);
                   if n > 100
                        break;
                   end
             end
                
             theta5(j) = theta5_current;
             theta6(j) = theta6_current;
             
             if j < length(theta2)
                 theta5(j+1) = theta5(j);
                 theta6(j+1) = theta6(j);
             end
        end
           
        R6_xcom = 202;
        R6_ycom = 150;
        R6mag = (R6_xcom^2 + R6_ycom^2)^(0.5);
        theta6com = 37.2 * pi/180;
        Xcom = zeros(length(theta2), 1);
        Xpos = zeros(length(theta2), 1);
        Ycom = zeros(length(theta2), 1);
        Ypos = zeros(length(theta2), 1);
        
        for i = 1:length(theta2)
            Xcom(i) = R6mag * cos(theta6com + theta6(i)) + R2* cos(theta2(i)) + (R6-R26)*cos(theta3(i));
            Ycom(i) = R6mag * sin(theta6com + theta6(i)) + R2* sin(theta2(i)) + (R6-R26)*sin(theta3(i));
            Xpos(i) = abs(Xcom(i) - Xcom(1));
            Ypos(i) = abs(Ycom(i) - Ycom(1));
        end

        fig1 = figure(fignum1);
        plot(rad2deg(theta2), Xpos);
        title('Horizontal Displacement of Tabletop COM')
        xlabel('θ2 (degress)')
        ylabel('Displacement (mm)')
        
        fig2 = figure(fignum2);
        plot(rad2deg(theta2), Ypos);
        title('Vertical Displacement of Tabletop COM')
        xlabel('θ2 (degrees)')
        ylabel('Displacement (mm)')

        if movie 
            saveas(fig1, 'xpos.jpg')
            saveas(fig2, 'ypos.jpg')
        end
    end
    
    function angularAnalysis(a, maxima, minima, fignum, movie)
        theta2_max = deg2rad(maxima);
        theta2_min = deg2rad(minima);
        theta2 = linspace(theta2_min, theta2_max, 1000)';
        R1 = 4*a;
        R2 = 5*a;
        R3 = 2*a; 
        R4 = 5*a; 
        R6 = 2*a; 
        R5 = 1.25001*R6;
        R15 = 0.5*R4;
        R26 = 0.5*R3;
        
        theta3 = ones(length(theta2),1);
        theta4 = ones(length(theta2),1);
        theta5 = ones(length(theta2),1);
        theta6 = ones(length(theta2),1);
        theta3 = theta3 .* (80 * pi/180);
        theta4 = theta4 .* (300 * pi/180);
        theta5 = theta5 .* (200 * pi/180);
        theta6 = theta6 .* (5 * pi/180);
    
        for j = 1:length(theta2)
            n = 0;
            D = [100; 100];
            theta3_current = theta3(j);
            theta4_current = theta4(j);
            while norm(D,1) > 10^-6 
                n = n + 1;
                func1 = R2*cos(theta2(j)) + R3*cos(theta3_current) + R4*cos(theta4_current) - R1;
                func2 = R2*sin(theta2(j)) + R3*sin(theta3_current) + R4*sin(theta4_current);
                F = [func1; func2];
                J11 = -R3*sin(theta3_current);
                J12 = -R4*sin(theta4_current);
                J21 = R3*cos(theta3_current);
                J22 = R4*cos(theta4_current);
                Jacob1 = [J11 J12; J21 J22];
                D = -Jacob1\F;
                dt3 = D(1);
                dt4 = D(2);
                theta3_current = theta3_current + dt3;
                theta4_current = theta4_current + dt4;
                if n > 100
                    break;
                end
            end
          
            theta3(j) = theta3_current;
            theta4(j) = theta4_current;
             if j < length(theta2)
                 theta3(j+1) = theta3(j);
                theta4(j+1) = theta4(j);
             end
        end
        
        for j = 1: length(theta2)
             n = 0;
             X = [100; 100];
             theta5_current = theta5(j);
             theta6_current = theta6(j);
             while norm(X,1) > 10^(-6)
                   n = n + 1;
                   func3 = R2*cos(theta2(j)) + R26*cos(theta3(j)) + R6*cos(theta6_current) + R5*cos(theta5_current) + R15*cos(theta4(j)) - R1;
                   func4 = R2*sin(theta2(j)) + R26*sin(theta3(j)) + R6*sin(theta6_current) + R5*sin(theta5_current) + R15*sin(theta4(j));
                   f = [func3;func4];
                   J11 = -R5*sin(theta5_current);
                   J12 = -R6*sin(theta6_current);
                   J21 = R5*cos(theta5_current);
                   J22 = R6*cos(theta6_current);
                   Jacob2 =  [J11 J12; J21 J22];
                   X = -Jacob2\f;
                   theta5_current = theta5_current + X(1);
                   theta6_current = theta6_current + X(2);
                   if n > 100
                        break;
                   end
             end
                
             theta5(j) = theta5_current;
             theta6(j) = theta6_current;
             
             if j < length(theta2)
                 theta5(j+1) = theta5(j);
                 theta6(j+1) = theta6(j);
             end
        end
           
        R6_xcom = 202;
        R6_ycom = 150;
        R6mag = (R6_xcom^2 + R6_ycom^2)^(0.5);
        theta6com = 37.2 * pi/180;
        Xcom = zeros(length(theta2), 1);
        Xpos = zeros(length(theta2), 1);
        Ycom = zeros(length(theta2), 1);
        Ypos = zeros(length(theta2), 1);
        
        for i = 1:length(theta2)
            Xcom(i) = R6mag * cos(theta6com + theta6(i)) + R2* cos(theta2(i)) + (R6-R26)*cos(theta3(i));
            Ycom(i) = R6mag * sin(theta6com + theta6(i)) + R2* sin(theta2(i)) + (R6-R26)*sin(theta3(i));
            Xpos(i) = abs(Xcom(i) - Xcom(1));
            Ypos(i) = abs(Ycom(i) - Ycom(1));
        end
        fig = figure(fignum);
        plot(rad2deg(theta2), rad2deg(theta6));
        title('Angular Rotation of Tabletop');
        xlabel('θ2 (degrees)');
        ylabel('Angle of Tabletop (degrees)');

        if movie
            saveas(fig, 'rad.jpg')
        end
    end

    function torqueAnalysis(a, maxima, minima, mass, fignum, movie)
        theta2_max = deg2rad(maxima);
        theta2_min = deg2rad(minima);
        theta2 = linspace(theta2_min, theta2_max, 1000);
        R1 = 4*a;
        R2 = 5*a;
        R3 = 2*a; 
        R4 = 5*a; 
        R6 = 2*a; 
        R5 = 1.25001*R6;
        R15 = 0.5*R4;
        R26 = 0.5*R3;
        
        links = [R1; R2; R3; R4; R5; R6];

        theta2 = transpose(theta2);
        theta3 = ones(length(theta2),1);
        theta4 = ones(length(theta2),1);
        theta5 = ones(length(theta2),1);
        theta6 = ones(length(theta2),1);
        theta3 = theta3 .* (80 * pi/180);
        theta4 = theta4 .* (300 * pi/180);
        theta5 = theta5 .* (200 * pi/180);
        theta6 = theta6 .* (5 * pi/180);
    
        for j = 1:length(theta2)
            n = 0;
            D = [100; 100];
            theta3_current = theta3(j);
            theta4_current = theta4(j);
            while norm(D,1) > 10^-6 
                n = n + 1;
                func1 = R2*cos(theta2(j)) + R3*cos(theta3_current) + R4*cos(theta4_current) - R1;
                func2 = R2*sin(theta2(j)) + R3*sin(theta3_current) + R4*sin(theta4_current);
                F = [func1; func2];
                J11 = -R3*sin(theta3_current);
                J12 = -R4*sin(theta4_current);
                J21 = R3*cos(theta3_current);
                J22 = R4*cos(theta4_current);
                Jacob1 = [J11 J12; J21 J22];
                D = -Jacob1\F;
                dt3 = D(1);
                dt4 = D(2);
                theta3_current = theta3_current + dt3;
                theta4_current = theta4_current + dt4;
                if n > 100
                    break;
                end
            end
          
            theta3(j) = theta3_current;
            theta4(j) = theta4_current;
             if j < length(theta2)
                 theta3(j+1) = theta3(j);
                theta4(j+1) = theta4(j);
             end
        end
        
        for j = 1: length(theta2)
             n = 0;
             X = [100; 100];
             theta5_current = theta5(j);
             theta6_current = theta6(j);
             while norm(X,1) > 10^(-6)
                   n = n + 1;
                   func3 = R2*cos(theta2(j)) + R26*cos(theta3(j)) + R6*cos(theta6_current) + R5*cos(theta5_current) + R15*cos(theta4(j)) - R1;
                   func4 = R2*sin(theta2(j)) + R26*sin(theta3(j)) + R6*sin(theta6_current) + R5*sin(theta5_current) + R15*sin(theta4(j));
                   f = [func3;func4];
                   J11 = -R5*sin(theta5_current);
                   J12 = -R6*sin(theta6_current);
                   J21 = R5*cos(theta5_current);
                   J22 = R6*cos(theta6_current);
                   Jacob2 =  [J11 J12; J21 J22];
                   X = -Jacob2\f;
                   theta5_current = theta5_current + X(1);
                   theta6_current = theta6_current + X(2);
                   if n > 100
                        break;
                   end
             end
                
             theta5(j) = theta5_current;
             theta6(j) = theta6_current;
             
             if j < length(theta2)
                 theta5(j+1) = theta5(j);
                 theta6(j+1) = theta6(j);
             end
        end
           
        R6_xcom = 202;
        R6_ycom = 150;
        R6mag = (R6_xcom^2 + R6_ycom^2)^(0.5);
        theta6com = 37.2 * pi/180;
        Xcom = zeros(length(theta2), 1);
        Xpos = zeros(length(theta2), 1);
        Ycom = zeros(length(theta2), 1);
        Ypos = zeros(length(theta2), 1);
        
        for i = 1:length(theta2)
            Xcom(i) = R6mag * cos(theta6com + theta6(i)) + R2* cos(theta2(i)) + (R6-R26)*cos(theta3(i));
            Ycom(i) = R6mag * sin(theta6com + theta6(i)) + R2* sin(theta2(i)) + (R6-R26)*sin(theta3(i));
            Xpos(i) = abs(Xcom(i) - Xcom(1));
            Ypos(i) = abs(Ycom(i) - Ycom(1));
        end
          
        I = (mass/12).*(30^2+links.^2);
        time = linspace(0,1000,1000)';
        omega = pi/180;
        motor_theta2 = 0.5*(theta2_max-theta2_min)*sin(omega*time) + 0.5*(theta2_max+theta2_min);
        motor_dtheta2 = 0.5*omega*(theta2_max-theta2_min)*cos(omega*time);
        motor_ddtheta2 = -0.5*omega^2*(theta2_max-theta2_min)*sin(omega*time);
        
        motor_theta3 = ones(length(motor_theta2),1);
        motor_theta4 = ones(length(motor_theta2),1);
        motor_theta5 = ones(length(motor_theta2),1);
        motor_theta6 = ones(length(motor_theta2),1);
        motor_theta3 = motor_theta3 .* (80 * pi/180);
        motor_theta4 = motor_theta4 .* (300 * pi/180);
        motor_theta5 = motor_theta5 .* (200 * pi/180);
        motor_theta6 = motor_theta6 .* (5 * pi/180);
    
        for j = 1:length(time)
            n = 0;
            D = [100; 100];
            theta3_current = motor_theta3(j);
            theta4_current = motor_theta4(j);
            while norm(D,1) > 10^-6 
                n = n + 1;
                func1 = R2*cos(motor_theta2(j)) + R3*cos(theta3_current) + R4*cos(theta4_current) - R1;
                func2 = R2*sin(motor_theta2(j)) + R3*sin(theta3_current) + R4*sin(theta4_current);
                F = [func1; func2];
                J11 = -R3*sin(theta3_current);
                J12 = -R4*sin(theta4_current);
                J21 = R3*cos(theta3_current);
                J22 = R4*cos(theta4_current);
                Jacob1 = [J11 J12; J21 J22];
                D = -Jacob1\F;
                dt3 = D(1);
                dt4 = D(2);
                theta3_current = theta3_current + dt3;
                theta4_current = theta4_current + dt4;
                if n > 100
                    break;
                end
            end
          
            motor_theta3(j) = theta3_current;
            motor_theta4(j) = theta4_current;
             if j < length(motor_theta2)
                motor_theta3(j+1) = motor_theta3(j);
                motor_theta4(j+1) = motor_theta4(j);
             end
        end
        
        for z = 1: length(time)
         n = 0;
         X = [100; 100];
         theta5_current = motor_theta5(z);
         theta6_current = motor_theta6(z);
         while norm(X,1) > 10^(-6)
               n = n + 1;
               func3 = R2*cos(motor_theta2(z)) + R26*cos(motor_theta3(z)) + R6*cos(theta6_current) + R5*cos(theta5_current) + R15*cos(motor_theta4(z)) - R1;
               func4 = R2*sin(motor_theta2(z)) + R26*sin(motor_theta3(z)) + R6*sin(theta6_current) + R5*sin(theta5_current) + R15*sin(motor_theta4(z));
               f = [func3;func4];
               J11 = -R5*sin(theta5_current);
               J12 = -R6*sin(theta6_current);
               J21 = R5*cos(theta5_current);
               J22 = R6*cos(theta6_current);
               Jacob2 =  [J11 J12; J21 J22];
               X = -Jacob2\f;
               theta5_current = theta5_current + X(1);
               theta6_current = theta6_current + X(2);
               if n > 100
                    break;
               end
         end
        
         motor_theta5(z) = theta5_current;
         motor_theta6(z) = theta6_current;
         
         if z < length(theta2)
         motor_theta5(z+1) = motor_theta5(z);
         motor_theta6(z+1) = motor_theta6(z);
         end
        end
    
        
        iterations = length(time);
        motor_theta3_prime = ones(length(iterations),1);
        motor_theta4_prime = ones(length(iterations),1);
        motor_theta5_prime = ones(length(iterations),1);
        motor_theta6_prime = ones(length(iterations),1);
        
        for i = 1:iterations
            Jacob1(1,1) = -R3*sin(motor_theta3(i,1));
            Jacob1(1,2) = -R4*sin(motor_theta4(i,1));
            Jacob1(2,1) = R3*cos(motor_theta3(i,1));
            Jacob1(2,2) = R4*cos(motor_theta4(i,1));
            c1 = R2*sin(motor_theta2(i,1));
            c2 = -R2*cos(motor_theta2(i,1));
            C1 = [c1;c2];
            KC1 = Jacob1\C1;
            motor_theta3_prime(i,1) = KC1(1);
            motor_theta4_prime(i,1) = KC1(2);
        end
       
        for i = 1:iterations
            Jacob2(1,1) = -R5*sin(motor_theta5(i,1));
            Jacob2(1,2) = -R6*sin(motor_theta6(i,1));
            Jacob2(2,1) = R5*cos(motor_theta5(i,1));
            Jacob2(2,2) = R6*cos(motor_theta6(i,1));
            c3 = R15*motor_theta4_prime(i,1)*sin(motor_theta4(i,1))+R2*sin(motor_theta2(i,1))+R26*motor_theta3_prime(i,1)*sin(motor_theta3(i,1));
            c4 = -R15*motor_theta4_prime(i,1)*cos(motor_theta4(i,1))-R2*cos(motor_theta2(i,1))-R26*motor_theta3_prime(i,1)*cos(motor_theta3(i,1));
            C2 = [c3;c4];
            SecondOrder_KC2 = Jacob2\C2;
            motor_theta5_prime(i,1) = SecondOrder_KC2(1);
            motor_theta6_prime(i,1) = SecondOrder_KC2(2);
        end
        
        motor_theta3_Dprime = zeros(length(iterations),1);
        motor_theta4_Dprime = zeros(length(iterations),1);
        motor_theta5_Dprime = zeros(length(iterations),1);
        motor_theta6_Dprime = zeros(length(iterations),1);
    
        for i = 1:iterations
            Jacob1(1,1) = -R3*sin(motor_theta3(i,1));
            Jacob1(1,2) = -R4*sin(motor_theta4(i,1));
            Jacob1(2,1) = R3*cos(motor_theta3(i,1));
            Jacob1(2,2) = R4*cos(motor_theta4(i,1));
            c1 = R2*cos(motor_theta2(i,1))+R3*motor_theta3_prime(i,1)^2*cos(motor_theta3(i,1))+R4*motor_theta4_prime(i,1)^2*cos(motor_theta4(i,1));
            c2 = R2*sin(motor_theta2(i,1))+R3*motor_theta3_prime(i,1)^2*sin(motor_theta3(i,1))+R4*motor_theta4_prime(i,1)^2*sin(motor_theta4(i,1));
            C1 = [c1;c2];
            SecondOrder_KC1 = Jacob1\C1;
            motor_theta3_Dprime(i,1) = SecondOrder_KC1(1);
            motor_theta4_Dprime(i,1) = SecondOrder_KC1(2);
        end
        
        for i = 1:iterations
            Jacob2(1,1) = -R5*sin(motor_theta5(i,1));
            Jacob2(1,2) = -R6*sin(motor_theta6(i,1));
            Jacob2(2,1) = R5*cos(motor_theta5(i,1));
            Jacob2(2,2) = R6*cos(motor_theta6(i,1));
            c3 = R5*motor_theta5_prime(i,1)^2*cos(motor_theta5(i,1))+R15*motor_theta4_Dprime(i,1)*sin(motor_theta4(i,1))+R15*motor_theta4_prime(i,1)^2*cos(motor_theta4(i,1)) +R2*cos(motor_theta2(i,1))+R26*motor_theta3_Dprime(i,1)*sin(motor_theta3(i,1))+R26*motor_theta3_prime(i,1)^2*cos(motor_theta3(i,1))+R6*motor_theta6_prime(i,1)^2*cos(motor_theta6(i,1));
            c4 = R5*motor_theta5_prime(i,1)^2*sin(motor_theta5(i,1))-R15*motor_theta4_Dprime(i,1)*cos(motor_theta4(i,1))+R15*motor_theta4_prime(i,1)^2*sin(motor_theta4(i,1)) +R2*sin(motor_theta2(i,1))-R26*motor_theta3_Dprime(i,1)*cos(motor_theta3(i,1))+R26*motor_theta3_prime(i,1)^2*sin(motor_theta3(i,1))+R6*motor_theta6_prime(i,1)^2*sin(motor_theta6(i,1));
            C2 = [c3;c4];
            SecondOrder_KC2 = Jacob2\C2;
            motor_theta5_Dprime(i,1) = SecondOrder_KC2(1);
            motor_theta6_Dprime(i,1) = SecondOrder_KC2(2);
        end
        xprime_G2 = -0.5*R2*sin(motor_theta2)/1000;
        xDprime_G2 = -0.5*R2*cos(motor_theta2)/1000;
        yprime_G2 = 0.5*R2*cos(motor_theta2)/1000;
        yDprime_G2 = -0.5*R2*sin(motor_theta2)/1000;
        xprime_G3 = (-R2*sin(motor_theta2)-0.5*R3*motor_theta3_prime.*sin(motor_theta3))/1000;
        xDprime_G3 = (-R2*cos(motor_theta2)-0.5*R3*motor_theta3_Dprime.*sin(motor_theta3)-0.5*R3*motor_theta3_prime.^2.*cos(motor_theta3))/1000;
        yprime_G3 = (R2*cos(motor_theta3)+0.5*R3*motor_theta3_prime.*cos(motor_theta3))/1000;
        yDprime_G3 = (-R2*sin(motor_theta3)+R3*motor_theta3_Dprime.*cos(motor_theta3)-0.5*R3*motor_theta3_prime.^2.*sin(motor_theta3))/1000;
        xprime_G4 = (-R2.*sin(motor_theta2)-R3*motor_theta3_prime.*sin(motor_theta3)-0.5.*R4*motor_theta4_Dprime.*sin(motor_theta4))/1000;
        xDprime_G4 = (-R2.*cos(motor_theta2)-R3*motor_theta3_Dprime.*sin(motor_theta3)-R3*motor_theta3_prime.^2.*cos(motor_theta3)- 0.5*R4*motor_theta4_Dprime.*sin(motor_theta4)-0.5*R4*motor_theta4_Dprime.^2.*cos(motor_theta4))/1000;
        yprime_G4 = (R2*cos(motor_theta2)+R3*motor_theta3_prime.*cos(motor_theta3)+ 0.5.*R4*motor_theta4_prime.*cos(motor_theta4))/1000;
        yDprime_G4 = (-R2*sin(motor_theta2)+R3*motor_theta3_Dprime.*cos(motor_theta3)-R3*motor_theta3_prime.^2.*sin(motor_theta3)+0.5*R4*motor_theta4_Dprime.*cos(motor_theta4)-0.5*R4*motor_theta4_Dprime.^2.*sin(motor_theta4))/1000;
        xprime_G5 = xprime_G3 +(-R6*motor_theta6_prime.*sin(motor_theta6)-0.5*R5*motor_theta5_prime.*sin(theta5))/1000;
        xDprime_G5 = xDprime_G3 +(-R6*motor_theta6_Dprime.*sin(motor_theta6)-R6*motor_theta6_prime.^2.*cos(motor_theta6)-0.5*R5*motor_theta5_Dprime.*sin(motor_theta5)-0.5*R5*motor_theta5_prime.^2.*cos(motor_theta5))/1000;
        yprime_G5 = yprime_G3 +(R6*motor_theta6_prime.*cos(motor_theta6)+0.5*R5*motor_theta5_prime.*cos(motor_theta5))/1000;
        yDprime_G5 = yDprime_G3 +(R6*motor_theta6_Dprime.*cos(motor_theta6)-R6*motor_theta6_prime.^2.*sin(motor_theta6)+0.5*R5*motor_theta5_Dprime.*cos(motor_theta5) - 0.5*R5*motor_theta5_prime.^2.*sin(motor_theta5))/1000;
        xprime_G6 = (-R2*sin(motor_theta2)-R26*motor_theta3_prime.*sin(motor_theta3)-R6mag*motor_theta6_prime.*sin(motor_theta6))/1000;
        xDprime_G6 = (-R2*cos(motor_theta2)-R26*motor_theta3_Dprime.*sin(motor_theta3)-R26*motor_theta3_prime.^2.*cos(motor_theta3)-R6mag*motor_theta6_Dprime.*cos(motor_theta6)-R6mag*motor_theta6_prime.^2.*sin(motor_theta6))/1000;
        yprime_G6 = (R2*cos(motor_theta2)+R26*motor_theta3_prime.*cos(motor_theta3)+R6mag*motor_theta6_prime.*cos(motor_theta6))/1000;
        yDprime_G6 = (-R2*sin(motor_theta2)+R26*motor_theta3_Dprime.*cos(motor_theta3)-R26*motor_theta3_prime.^2.*sin(motor_theta3)+R6mag*motor_theta6_Dprime.*cos(motor_theta6)-R6mag*motor_theta6_prime.^2.*sin(motor_theta6))/1000;
        xprime_G = [zeros(length(xprime_G2),1), xprime_G2, xprime_G3, xprime_G4, xprime_G5, xprime_G6];
        yprime_G = [zeros(length(yprime_G2),1), yprime_G2, yprime_G3, yprime_G4, yprime_G5, yprime_G6];    
        xDprime_G = [zeros(length(xDprime_G2),1), xDprime_G2, xDprime_G3, xDprime_G4, xDprime_G5, xDprime_G6];
        yDprime_G = [zeros(length(yDprime_G2),1), yDprime_G2, yDprime_G3, yDprime_G4, yDprime_G5, yDprime_G6];
        motor_theta2_prime = ones(length(motor_theta3_prime),1);
        motor_theta2_Dprime = zeros(length(motor_theta3_prime),1);
        motor_thetaj_prime = [zeros(length(motor_theta3_prime)), motor_theta2_prime, motor_theta3_prime, motor_theta4_prime, motor_theta5_prime, motor_theta6_prime];
        motor_thetaj_Dprime = [zeros(length(motor_theta3_Dprime)), motor_theta2_Dprime, motor_theta3_Dprime, motor_theta4_Dprime, motor_theta5_Dprime, motor_theta6_Dprime];
        
        A = zeros(length(motor_theta2),6);
        Asum = zeros(6,1);
        B = zeros(length(motor_theta2),6);
        Bsum = zeros(6,1);
        for i = 1: length(motor_theta2)
            for j = 2: 6
                A(i,j) = mass(j) .* (xprime_G(i,j).^2 + yprime_G(i,j).^2) + I(j) .* motor_thetaj_prime(i,j)^2;
                Asum(j) = norm(A(j,:),1);
                B(i,j) = mass(j) .* (xprime_G(i,j).*xDprime_G(i,j) + yprime_G(i,j).*yDprime_G(i,j)) + I(j).*motor_thetaj_prime(i,j).*motor_thetaj_Dprime(i,j);
                Bsum(j) = norm(B(j,:),1);
            end
        end
    
        Bsum(1,1) = norm(Bsum(:,1),1);
        
        g = 9810; 
        g2 = mass(2)*g*yprime_G2;
        g3 = mass(3)*g*yprime_G3;
        g4 = mass(4)*g*yprime_G4;
        g5 = mass(5)*g*yprime_G5;
        g6 = mass(6)*g*yprime_G6;
        g_tot = (g2 + g3 + g4 + g5 + g6)/1000;
        
        T = Asum(1,1).*motor_ddtheta2 + Bsum(1,1).*motor_dtheta2.^2;
        
        torque = (T + g_tot)/1000;
        
        fig = figure(fignum);
        plot(time, torque);
        title('Torque over Three Osciliitory Periods');
        xlabel('Time (s)');
        ylabel('Torque (Nm)');

        if movie
            saveas(fig, 'torque.jpg')
        end
end

    % End Methods
    end
end


