% Joshua Guthrie, Charlebois Laboratory, University of Alberta, Department of Physics, jdguthri@ualberta.ca
% Deterministic ODE simulation code for the manuscript "Phenotypic Heterogeneity Facilitates Survival WhileHindering the Evolution of Drug Resistance"


clear;

% define constant simulation parameters
relative_to_N0 = false; % true if using N_o (initial total pop. size) as the base measure for the fixation/establishment calculations
log_bool = true; % true if plotting time series on a loglog scale
save_time_series = true;
plot_heatmaps = true;
linspace_n = 15; % heatmap grid size
save_heatmaps = true;
deathrate_sweep = false; % experimental feature, ignore
save_to = ""; % directory to save output files to

% simulation parameters
t_end = 150; dt = 0.001; %hours
S_i = 5.5e+5; G1_i = 0.0; G2_i = 0.0; %G2 is an experimental feature, ignore

% define the globals
global N0 k n kS kN kG1 kG2 rSN rNS rG1S rG1N rG2G1 dS dN dG1 dG2
k = 1e+7; n = 2; % Baryani-Hill function terms for model 2
kS = 0.0; 
rSN = 0.0035; rNS = 0.0625; rG1S = 0.0; rG1N = (1e-6)/3; rG2G1 = (1e-6)/3; 
dS = 0.0; dN = 0.0; dG1 = 0.0; dG2 = 0.0;

% script code to run multiple simulations, not needed for main functionality
models = [2]; % 1 for exponential growth, 2 for logistic growth
scenarios = [1]; % 1 for G1 only, 2 for G1 and G2
cidal_list = [true];
deathfactor_list = [1];
G1_susceptible_list = [false]; 
G2_susceptible_list = [false];
Ni_list = [5.5e+4]; 

for a = 1:length(Ni_list)
    N_i = Ni_list(a);
    N0 = S_i  + N_i + G1_i + G2_i;
    for b = 1:length(models)
        model = models(b);
        for c = 1:length(scenarios)
            scenario = scenarios(c);
            for d = 1:length(cidal_list)
                cidal = cidal_list(d);
                if cidal == false
                    deathfactor = 0.0;
                    G1_susceptible = false;
                    G2_susceptible = false;
                    deterministic_simulations(model, scenario, relative_to_N0, cidal, deathfactor, ... 
                                              G1_susceptible, G2_susceptible, log_bool, save_time_series, ... 
                                              plot_heatmaps, linspace_n, save_heatmaps, deathrate_sweep, ...
                                              t_end, dt, S_i, N_i, G1_i, G2_i, save_to)
                else
                    for e = 1:length(deathfactor_list)
                        deathfactor = deathfactor_list(e);
                        for f = 1:length(G1_susceptible_list)
                            G1_susceptible = G1_susceptible_list(f);
                            for g = 1:length(G2_susceptible_list)
                                G2_susceptible = G2_susceptible_list(g);
                                deterministic_simulations(model, scenario, relative_to_N0, cidal, deathfactor, ... 
                                                          G1_susceptible, G2_susceptible, log_bool, save_time_series, ... 
                                                          plot_heatmaps, linspace_n, save_heatmaps, deathrate_sweep, ...
                                                          t_end, dt, S_i, N_i, G1_i, G2_i, save_to)
                            end
                        end
                    end
                end
            end
        end
    end
end

