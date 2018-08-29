%INITIALISE TEMPERATURE USING KIRKPATRICK METHOD
function T0 = InitialiseTemperature(Dimension,x,D_initial)

%Initial temperature paremeters
Inital_acceptance_probability = 0.80; %Suggested value
Inital_temperature_iterations = 1000;
Number_of_increases = 0;
Df_storage = zeros(Inital_temperature_iterations,1);

for search = 1:Inital_temperature_iterations
    
    %Initialise stepsize matrix and find search point
    D = diag(D_initial.*ones(Dimension,1));
    u = 2.*(rand(Dimension,1)-0.50);
    x_search = x+D*u;
    %Keep generating new samples if they violate constraints
    while sum(x_search<-512)>0 || sum(x_search>512)>0
        u = 2.*(rand(Dimension,1)-0.50);
        x_search = x+D*u;
    end
    
    %Store values of objective function increases
    df = objective(x_search, Dimension)-objective(x, Dimension);
    if df > 0
        Number_of_increases = Number_of_increases + 1;
        Df_storage(Number_of_increases) = df;
    end
end

%Calculate initial temperature by accepting increases with the acceptence
%probability
if Number_of_increases > 0
    Df_avg = mean(Df_storage);
    T0 = -Df_avg/log(Inital_acceptance_probability);
else
    disp('No objective increases found, T0 = 1')
    T0 = 1;
end

end

