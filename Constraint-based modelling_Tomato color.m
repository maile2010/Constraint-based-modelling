% Define the stoichiometric matrix (example placeholder) 
% Rows correspond to metabolites, columns to reactions 
S = [ 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      0 1 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      0 0 3 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      0 0 0 1 1 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0; 
      0 0 0 0 0 0 1 0 0 -1 -1 -1 0 0 0 0 0 0 0 0 0; 
      0 0 0 0 0 0 0 0 0 0 1 0 -1 0 0 0 0 0 0 0 0; 
      0 0 0 0 0 0 0 0 0 0 0 0 1 -1 -1 0 0 0 0 0 0; 
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0; 
      0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 -1 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1];
% Define the objective function coefficients 
% We maximize v7 + v13 + v18 + v20 -> minimize -v7 - v13 - v18 - v20 
f = zeros(21, 1); % 21 reactions 
f(7) = -1; % v7 (red) 
f(13) = -1; % v13 (orange) 
f(18) = -1; % v18 (orange) 
f(20) = -1; % v20 (yellow)
% Define the bounds for each reaction flux 
lb = zeros(21, 1); % Lower bounds: irreversible reactions have flux >= 0 
ub = ones(21, 1)*1000; % Upper bounds: initialize as 1000
%set bound for reversible reaction 6
lb(6) = -1000;
% Apply given upper bounds 
ub(1) = 10; % v1 <= 10 mmol/(gDW*h)
% Set fixed flux values 
lb(10) = 5; % v10 = 5 mmol/(gDW*h) 
ub(10) = 5; % Fixed upper and lower bounds for v10 
lb(17) = 3; % v17 = 3 mmol/(gDW*h) 
ub(17) = 3; % Fixed upper and lower bounds for v17 
%knock-out v11
%ub(11) = 0;
%lb(11) = 0;
% Enforce v11 = v12 using additional constraints 
A_eq = S; % Start with the stoichiometric constraints 
b_eq = zeros(size(S, 1), 1); % Steady-state constraints 
% Add v11 = v12 constraint 
A_eq = [A_eq; zeros(1, 21)]; % Add a new row for v11 = v12 
A_eq(end, 11) = 1; % Coefficient for v11 
A_eq(end, 12) = -1; % Coefficient for -v12 
b_eq = [b_eq; 0]; % Right-hand side for v11 = v12 constraint 
% Solver options
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'iter');
% Solve the linear program
[x, fval, exitflag, output] = linprog(f, [], [], A_eq, b_eq, lb, ub, options);
    % Output results
disp('Optimized flux distribution:');
disp(x);
disp('Maximized carotenoid production (total flux):');
disp(-fval); % Negate because linprog minimizes by default



