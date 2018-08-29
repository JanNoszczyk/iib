%EVOLUTION STRATEGIES ALGORITHM

close all

%Main Parameters
Dimension = 5;
Mu_Lambda_ratio = 7; %Recommended by Schwefel 1:7
Parent_size = 68;
Offspring_size = Parent_size*Mu_Lambda_ratio;
Maximum_runs = 300;
Maximum_evaluations = 10000;
Optimal_objective = -3800; %Desired value of objective function
%Calculated maximum number of generations based on the limit of evaluations
Max_generations = ceil((Maximum_evaluations-Parent_size)/(Parent_size+Offspring_size));
%1 for (mu, lambda), 2 for (mu + lambda)
Selection_type = 2;
%Mutation Parameters
Tau   = 1/(sqrt(2*sqrt(Dimension)));
Tau_prime  = 1/(sqrt(2*Dimension));
Minimum_sd = 0.1;
Mutation_update = 2; %1 for uncorralated, 2 for correlated
Beta  = 0.0873; %Suggested value

%Recombination Methods
%Control variable recombination operator, 1 for discrete, 2 for global
%discrete
Control_recombination = 2;
%Strategy parameters recombination operator, 1 for intermediate, 2 for
%global intermediate
Strategy_recombination = 2;
Intermediate_weight = 0.5; %Suggested value

%Initialise value storage for each run
X_storage = cell(1,Maximum_runs); 
Current_objective = cell(1,Maximum_runs); 
Best_objective = zeros(Parent_size,Maximum_runs); 
Average_obj_population = zeros(Maximum_runs, Max_generations);
Minimum_obj_population = zeros(Maximum_runs, Max_generations);

for Run = 1:Maximum_runs
    %Ensure same random seeds are used at each run
    rand('seed',Run); randn('seed',Run);
    
    % Initialise x within -512 and 512
    x = 1024*(rand([Dimension,Parent_size])-0.5);
    
    %Initialise covariance matrix
    Cov = cell(1,Parent_size); %(sigma)
    for i = 1:Parent_size
        Random_matrix = rand(Dimension);
        %Ensures covariance matrix is symmetric
        Cov{i} = Random_matrix+Random_matrix'; 
    end
    
    %Calculate corresponding rotation angles
    Angles = cell(1,Parent_size);
    for i = 1:Parent_size
        Angles{i} = zeros(Dimension);
        for j = 1:Dimension
            for k = 1:Dimension
                %Use atan2 for rotation angle to be between [-pi,pi]
                Angles{i}(j,k) = 0.5*atan2(2*Cov{i}(j,k),(Cov{i}(j,j) - Cov{i}(k,k)));
            end
        end
    end
    
    %Inialise Parent Error
    Parent_objective = zeros(Parent_size,1);
    for i = 1:Parent_size
        Parent_objective(i) = objective(x(:,i), Dimension);
    end
    Parent_error = abs(Optimal_objective - Parent_objective);
    
    Current_objective{Run} = Parent_objective;
    
    Parent_x = x;
    Parent_cov = Cov;
    Parent_angles = Angles;
    
    Evaluations = Parent_size;
    Generation = 1;
    Infeasible_meter = 0;
    
    %Initialise temporary storage tensor for x values
    X_temporary_storage = zeros(Dimension, Parent_size, Max_generations);
    X_temporary_storage(:,:,1) = x;
    
    
    %START EVOLUTION STRATEGY
    while Evaluations <= Maximum_evaluations
        
        %RECOMBINATION
        
        %CONTROL VARIABLES
        %Discrete Recombination
        if Control_recombination == 1
            for i = 1:Offspring_size
                %Select two parents randomly
                Random_parents = randsample(1:Parent_size,2);
                for j = 1:Dimension
                    x(j,i) = x(j,(randsample(Random_parents,1)));
                end
            end
        end
        %Global Discrete Recombination
        if Control_recombination == 2
            for i = 1:Offspring_size
                Random_parents = randsample(1:Parent_size,Parent_size);
                for j = 1:Dimension
                    x(j,i) = x(j,(randsample(Random_parents,1)));
                end
            end
        end
        
        %STRATEGY PARAMETERS
        %Intermediate Recombination
        if Strategy_recombination == 1
            for i = 1:Offspring_size
                Random_parents = randsample(1:Parent_size,2);
                Cov{i} = Intermediate_weight*Cov{Random_parents(1)}+(1-Intermediate_weight)*Cov{Random_parents(2)};
                Angles{i} = Intermediate_weight*Angles{Random_parents(1)}+(1-Intermediate_weight)*Angles{Random_parents(2)};
            end
        end
        %Global Intermediate Recombination
        if Strategy_recombination == 2
            for i = 1:Offspring_size
                Random_parents = randsample(1:Parent_size,Parent_size);
                Cov{i} = Intermediate_weight*Cov{Random_parents(1)}+(1-Intermediate_weight)*Cov{Random_parents(2)};
                Angles{i} = Intermediate_weight*Angles{Random_parents(1)}+(1-Intermediate_weight)*Angles{Random_parents(2)};
            end
        end
        
        
        %MUTATION
        for i = 1:Offspring_size
            Cov{i} = Cov{i}.*exp(Tau_prime*randn(Dimension) + Tau*randn(Dimension));
            Angles{i} = Angles{i} + Beta*randn(Dimension);
            
            %Uncorellated variable update
            if Mutation_update == 1
                x(:,i) = x(:,i) +randn(Dimension,1).*sqrt(diag(Cov{i}));
            end
            
            %Cottelated variable update (Rotation)
            if Mutation_update == 2
                R = eye(Dimension);
                for j = 1:Dimension-1
                    for k = j+1:Dimension
                        R_temp = eye(Dimension);
                        R_temp(j,j) = cos(Angles{i}(j,k));
                        R_temp(k,k) = cos(Angles{i}(j,k));
                        R_temp(j,k) = -sin(Angles{i}(j,k));
                        R_temp(k,j) = sin(Angles{i}(j,k));
                        R = R*R_temp;
                    end
                end
                
                x(:,i) = x(:,i) + R*sqrt(diag(diag(Cov{i})))*randn(Dimension,1);
            end
        end
        
        %Asses rotation angles and covariances
        for i = 1:Offspring_size
            %Adjust angles outside [-pi,pi]
            [A,B] = find(abs(Angles{i})>pi);
            Angles{i}(A,B) = Angles{i}(A,B) - 2*pi*(Angles{i}(A,B)/abs(Angles{i}(A,B)));
            %Reset negative covariances to constant value
            [A,B] = find(Cov{i}<=0);
            Cov{i}(A,B) = Minimum_sd;
        end
        
        %SELECTION
        %Evaluate Offspring Error
        Offspring_objective = zeros(Offspring_size,1);
        for i = 1:Offspring_size
            Offspring_objective(i) = objective(x(:,i), Dimension);
        end
        Offspring_error = abs(Optimal_objective - Offspring_objective);
        
        %(Mu, Lambda)-Selection
        if Selection_type == 1
            %Obtain indices for minimum Offpsring error
            [X_sorted, Index] = sort(Offspring_error);
            x = x(:,Index(1:Parent_size));
            Cov = Cov(Index(1:Parent_size));
            Angles = Angles(Index(1:Parent_size));
        end
        
        
        %(Mu + Lambda)-Selection
        if Selection_type == 2
            %Store parent and offspring errors and variables together
            Combined_x = [x Parent_x];
            Combined_cov = [Cov Parent_cov];
            Combined_angles = [Angles Parent_angles];
            Combined_error = [Offspring_error' Parent_error'];
            
            [X_sorted, Index] = sort(Combined_error);
            x = Combined_x(:,Index(1:Parent_size));
            Cov = Combined_cov(Index(1:Parent_size));
            Angles = Combined_angles(Index(1:Parent_size));
        end
        
        %Evaluate Parent Error
        Parent_objective = zeros(Parent_size,1);
        for i = 1:Parent_size
            Parent_objective(i) = objective(x(:,i), Dimension);
        end
        Parent_error = abs(Optimal_objective - Parent_objective);
        
        %Handle constraints
        %Reject solutions outside [-512,512], keeping previous ones
        for i = 1:Parent_size
            for j = 1:Dimension
                if x(j,i) >= 512 || x(j,i) <= -512
                    x = Parent_x;
                    Cov = Parent_cov;
                    Angles = Parent_angles;
                    Infeasible_meter = Infeasible_meter+1;
                end
            end
        end
        
        if Infeasible_meter > 1000
            error('No feasible x found')
        end
        
        %Update storage arrays
        Average_obj_population(Run, Generation) = mean(Parent_objective);
        Minimum_obj_population(Run, Generation) = min(Parent_objective);
        
        Evaluations = Evaluations + Parent_size + Offspring_size;
        Generation = Generation + 1;
        
        X_temporary_storage(:,:,Generation) = x;
        Current_objective{Run} = [Current_objective{Run} Parent_objective];
        
        Parent_x = x;
        Parent_cov = Cov;
        Parent_angles = Angles;
        
        
    end
    
    %Store control variables and best objective for each run
    X_storage{Run} = X_temporary_storage;
    for i = 1:Parent_size
        Best_objective(i,Run) = min(Current_objective{Run}(i,:));
    end
    
end



% Final Statistics
Best_parent_obj = zeros(1,Maximum_runs);
for i = 1:Maximum_runs
    Best_parent_obj(i) = mean(Best_objective(:,i));
end
Objective_run_mean = mean(Best_parent_obj)
Objective_run_std = std(Best_parent_obj)


% Final Plots
%Identify the run with the best objective function
[c, ix] = min(Best_objective(:));
[Best_parent, Best_run] = ind2sub(size(Best_objective), ix)
figure(1);
plot(1:length(Current_objective{Best_run}(Best_parent,:)),Current_objective{Best_run}(Best_parent,:))
xlabel('Generation')
ylabel('Objective Function Value')

Final_guy = (Current_objective{Best_run}(Best_parent,end));

%Plot average and minimum objective function values
Average_obj = zeros(1,Max_generations);
Minimum_obj = zeros(1,Max_generations);
for i = 1:Max_generations
    Average_obj(i) = mean(Average_obj_population(:,i));
    Minimum_obj(i) = min(Minimum_obj_population(:,i));
end
figure(2)
hold on
yyaxis left
ylabel('Average Objective Function')
plot(1:Max_generations, Average_obj)
yyaxis right
plot (1:Max_generations, Minimum_obj)
ylabel('Minimum Objective Function')
xlim([1 Max_generations])
xlabel('Generation')

%Plot the search pattern of the 2D-EF (2D Eggholder Function)
if Dimension == 2
    figure(3);
    hold on
    for i = 1:Parent_size
        x1_scatter = X_storage{Best_run}(1,i,:);
        x2_scatter = X_storage{Best_run}(2,i,:);
        xlabel('x_{1}')
        ylabel('x_{2}')
        title('Eggholder 2D Function')
        scatter(x1_scatter,x2_scatter,'r')
    end
    [X1,X2] = meshgrid(linspace(-512,512,100),linspace(-512,512,100));
    Z = zeros(100,100);
    for i = 1:100
        for j = 1:100
            Z(i,j) = objective([X1(i,j) X2(i,j)], Dimension);
        end
    end
    contour(X1,X2,Z,10)
end


