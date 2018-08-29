%SIMULATED ANNEALING ALGORITHM

%Main Parameters
Dimension = 5;
Maximum_evaluations = 10000;
Maximum_runs = 300;
%Stepsize parameters
Alpha = 0.1; %Suggested value
Omega = 2.1; %Suggested value
D_initial = 12;
D_min = 40;
D_max = 280;
%ECS parameters
Cooling_constant = 0.932;
%Markov chain length
L_k = 130;

%Initialise value storage for each run
Best_objective_storage = zeros(1,Maximum_runs);
X_storage = cell(1,Maximum_runs);
Current_objective = cell(1,Maximum_runs);
Infeasible_storage = zeros(1,Maximum_runs);
Average_obj_population = zeros(Maximum_runs, 10066);

for Run = 1:Maximum_runs
    %Ensure same random seeds are used at each run
    rand('seed',Run); randn('seed',Run)
    
    % Initialise x within -512 and 512
    x = 1024*(rand(Dimension,1)-0.5);
    %Initialise stepsize
    D = diag(D_initial*rand(Dimension,1).*ones(Dimension,1));
    %Initialise Temperature
    T0 = InitialiseTemperature(Dimension,x,D_initial);
    
    X_storage{Run} = x;
    Current_objective{Run} = objective(x, Dimension);
    
    T = T0;
    Iteration = 1;
    Evaluations = 1;
    
    while Evaluations <= Maximum_evaluations
        Last_objective_store = objective(x, Dimension);
        for L = 1:L_k
            %Generate random u vector in range [-1,1] and update x
            u = 2*(rand(Dimension,1)-0.5);
            x_new = x+D*u;
            
            %If solution is outside of constraints, generate new x
            Infeasible_meter = 0;
            while sum(x_new<-512)>0 || sum(x_new>512)>0
                %Stop search if new solutions consistently infeasible
                if Infeasible_meter > 1000
                    disp('No feasible x found')
                    return
                end
                u = 2*(rand(Dimension,1)-0.5);
                x_new = x+D*u;
                Infeasible_meter = Infeasible_meter + 1;
            end
            
            %Update actual step size
            D_actual = sqrt(sum(D*u.^2));
            
            %If objective decreases always accept
            %If increases, accept with acceptance probability
            f_new = objective(x_new, Dimension);
            f_old = Last_objective_store(L)
            Evaluations = Evaluations + 1;
            if f_new<f_old
                Acceptance_probability = 1;
            else
                Acceptance_probability = exp(-(f_new-f_old)/(T*D_actual));
            end
            
            %Compare acceptance probability with random number
            Random_probability = rand;
            if Random_probability < Acceptance_probability
                %Accept new x and update step size
                x = x_new;
                R = diag(abs(D*u));
                D = (1-Alpha)*D+Alpha*Omega*R;
                Last_objective_store = [Last_objective_store f_new];
                %Check if D is within upper and lower limits
                for k = 1:Dimension
                    if abs(D(k,k)) < D_min
                        D(k,k) = D_min;
                    elseif abs(D(k,k)) > D_max
                        D(k,k) = D_max;
                    end
                end
            else
                Last_objective_store = [Last_objective_store f_old];
            end
            X_storage{Run} = [X_storage{Run} x];
            Current_objective{Run} = [Current_objective{Run} Last_objective_store(end)];
            Average_obj_population(Run,Evaluations) = Last_objective_store(end);
        end
        Iteration = Iteration + 1;
        
        % Find new temperature with Exponential Cooling Scheme (ECS)
        T = Cooling_constant*T;
        
    end
    Best_objective_storage(Run) = min(Current_objective{Run}(:));
end

% Final Statistics
Objective_run_mean = mean(Best_objective_storage)
Objective_run_std = std(Best_objective_storage)
Infeasible_mean = mean(Infeasible_storage)
%Identify the run with the best objective function and the best final value
Best_objective_index = find(Best_objective_storage==min(Best_objective_storage));
Best_final_value = (Current_objective{Best_objective_index}(end));

% Final Plots
%Plot objective function value for the best run
figure(1);
plot(1:length(Current_objective{Best_objective_index}),Current_objective{Best_objective_index})
xlim([0 10000])
xlabel('Evaluation Number')
ylabel('Objective Function Value')

%Plot the search pattern of the 2D-EF (2D Eggholder Function)
if Dimension == 2
    figure(2);
    hold on
    x1_scatter = X_storage{Best_objective_index}(1,:);
    x2_scatter = X_storage{Best_objective_index}(2,:);
    xlabel('x_{1}')
    ylabel('x_{2}')
    title('Eggholder 2D Function')
    scatter(x1_scatter,x2_scatter,'r')
    [X1,X2] = meshgrid(linspace(-512,512,100),linspace(-512,512,100));
    Z = zeros(100,100);
    for i = 1:100
        for j = 1:100
            Z(i,j) = objective([X1(i,j) X2(i,j)], Dimension);
        end
    end
    contour(X1,X2,Z,10)
end

%Plot average and minimum objective function values
figure(3)
Average_obj = zeros(1,10066);
Minimum_obj = zeros(1,10066);
for i = 1:10066
    Average_obj(i) = mean(Average_obj_population(:,i));
    Minimum_obj(i) = min(Average_obj_population(:,i));
end
hold on
yyaxis left
ylabel('Average Objective Function')
plot(1:10066, Average_obj)
yyaxis right
plot (1:10066, Minimum_obj)
ylabel('Minimum Objective Function')
xlim([0 10000])
xlabel('Generation')
