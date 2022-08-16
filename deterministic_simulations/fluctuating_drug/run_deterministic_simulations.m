% Joshua Guthrie, University of Alberta
% Script to run simulations for different situations
clear;

% define constant simulation parameters
relative_to_N0 = false; %true if using N_o (initial total pop. size) as the base measure for the fixation/establishment calculations
log_bool = true;
save_time_series = true;
plot_heatmaps = false;
linspace_n = 15;
deathrate_sweep = false;
save_heatmaps = false;
save_to = "C:\Users\JoshuaGuthrie\cLab\SNG_July_2022\fluct_drug\";

t_end = 4800; dt = 0.01; %hours
S_i = 5.5e+5; G1_i = 0; G2_i = 0.0; 

% define the globals
global N0 k n kS kN kG1 kG2 rSN rNS rG1S rG1N rG2G1 dS dN dG1 dG2 kN_dsweep kG_dsweep
k = 1e+7; n = 2; % Baryani-Hill function terms for model 2
kS = 0; 
rSN = 0.0035; rNS = 0.0625; rG1S = 0; rG1N = (1e-6)/3; rG2G1 = (1e-6)/3; 
% basal death used for static drug case was the mean chronological life span of
% budding yeast in 2% gluclose solution (6.5 days = 156 hours in
% https://doi.org/10.1016/j.cmet.2012.06.002)
dS = 1/156; dN = 1/156; dG1 = 1/156; dG2 = 0.0;
%N0 = S_i + N_i + G1_i + G2_i;
kN_dsweep = 0.3466; kG_dsweep = 0.3466;

% define the simulation parameters
models = [2]; % 1 for exponential growth, 2 for logistic growth
scenarios = [1]; % 1 for G1 only, 2 for G1 and G2
cidal_list = [true];
deathfactor_list = [1];
G1_susceptible_list = [false]; 
G2_susceptible_list = [false];
Ni_list = [5.5e+4]; 

% deterministic_simulations(1, 1, true, false, 0, ... 
%                           false, false, true, true, ... 
%                           false, 15, false, false, ...
%                           t_end, dt, S_i, N_i, G1_i, G2_i)

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

