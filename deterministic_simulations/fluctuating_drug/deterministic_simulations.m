function deterministic_simulations(model, scenario, relative_to_N0, cidal, deathfactor, ... 
                                   G1_susceptible, G2_susceptible, log_bool, save_time_series, ... 
                                   plot_heatmaps, linspace_n, save_heatmaps, deathrate_sweep, ...
                                   t_end, dt, S_i, N_i, G1_i, G2_i, save_to)
                               
    % Joshua Guthrie, Charlebois Laboratory, University of Alberta
    % Simulation code using deterministic ODE models for quantitative
    % non-genetic/genetic drug resistance/evolution study. Covers numerical simulations
    % for multiple models, scenarios, and drug types, as well as result plotting.

    global N0 k n kS kN kG1 kG2 rSN rNS rG1S rG1N rG2G1 dS dN dG1 dG2 kN_dsweep kG_dsweep
    kN_list = [0.1733 0.2600 0.3466]; kG_list = [0.3466 0.3466 0.3466];
    dN_list = [0.1 0.5 1.0];
    dS_list = [1 1 1];
    

   
    %% Scenario 1
    % Run Scenario 1 (no G2) for different values of kG1 and kN
    if scenario == 1
        t_list = [];
        S_list = [];
        N_list = [];
        G1_list = [];
        dN_dt_list = [];
        t_est_list = [];
        t_fix_list = [];
        legend_list = [];

        % run simulations for each kN and kG
        for i = 1:length(kN_list)
            kN = kN_list(i); kG1 = kG_list(i);
            % if using cidal drugs, specify death rates
            if cidal
                dN = dN_list(i); 
                dS = dS_list(i); 
                if G1_susceptible
                    dG1 = 0.5;
                else
                    dG1 = 1/156;
                end
            else
                dS = 1/156; dN = 1/156; dG1 = 1/156; dG2 = 0.0;         
            end
            % run simulations, save results to lists
            [t, S, N, G1, ~, dN_dt, t_est, t_fix] = simulate(model,scenario,relative_to_N0,t_end,dt,S_i,N_i,G1_i,G2_i);
            if cidal
                legend_text = sprintf('\\delta_N = %.1f /hr', dN);
            else
                legend_text = sprintf('k_N = %.4f /hr', kN);
            end
            t_list = [t_list, t];
            S_list = [S_list, S];
            N_list = [N_list, N];
            G1_list = [G1_list, G1];
            dN_dt_list = [dN_dt_list, dN_dt];
            t_est_list = [t_est_list, t_est];
            t_fix_list = [t_fix_list, t_fix];
            legend_list = [legend_list, {legend_text}];
        end

        % plot simulation results
        % create the time series file name
        if save_time_series
           directory = save_to;
           if cidal
               drug_type = 'cidal';
               deathfactor_str = sprintf('_deathfactor%d',deathfactor);
           else
               drug_type = 'static';
               deathfactor_str ='';
           end
           if G1_susceptible
               G1_susceptible_str = '_G1susceptible';
           else
               G1_susceptible_str = '';
           end

           if G2_susceptible
               G2_susceptible_str = '_G2susceptible';
           else
               G2_susceptible_str = '';
           end

           if log_bool
               log_str = '_loglog';
           else
               log_str = '';
           end

           TS_filename = sprintf("%smodel%d_scenario%d_%s%s%s%s_Ni%0.1e%s_TimeSeries",directory,model,scenario,...
                                 drug_type,deathfactor_str,G1_susceptible_str,G2_susceptible_str, N_i, log_str);
        else
            TS_filename = 'n/a';
        end
        
        plot_simulation_results(t_list,S_list,N_list,G1_list,'n/a',dN_dt_list,legend_list,log_bool,save_time_series,TS_filename)

        % run large amount of simulations and plot t_est and t_fix heatmaps 
        if plot_heatmaps
            [kG_used_est, kG_used_fix, kN_used_est, kN_used_fix, t_est_list, t_fix_list] = population_times(model, scenario, relative_to_N0, t_end, dt,  ...
                                                                                                            S_i, N_i, G1_i, G2_i, cidal, G1_susceptible, ... 
                                                                                                            G2_susceptible, deathfactor, linspace_n);
            if save_heatmaps
               directory = save_to;
               if cidal
                   drug_type = 'cidal';
                   deathfactor_str = sprintf('_deathfactor%d',deathfactor);
               else
                   drug_type = 'static';
                   deathfactor_str ='';
               end
               if G1_susceptible
                   G1_susceptible_str = '_G1susceptible';
               else
                   G1_susceptible_str = '';
               end

               if G2_susceptible
                   G2_susceptible_str = '_G2susceptible';
               else
                   G2_susceptible_str = '';
               end

               HM_filename_est = sprintf("%smodel%d_scenario%d_%s%s%s%s_Ni%0.1e_establishment_heatmap",directory,model,scenario,...
                                     drug_type,deathfactor_str,G1_susceptible_str,G2_susceptible_str, N_i);
               HM_filename_fix = sprintf("%smodel%d_scenario%d_%s%s%s%s_Ni%0.1e_fixation_heatmap",directory,model,scenario,...
                                     drug_type,deathfactor_str,G1_susceptible_str,G2_susceptible_str, N_i);
            else
                HM_filename_est = 'n/a';
                HM_filename_fix = 'n/a';
            end            
            plot_population_times(kG_used_est, kN_used_est, t_est_list, 'Establishment', scenario, save_heatmaps, HM_filename_est)
            plot_population_times(kG_used_fix, kN_used_fix, t_fix_list, 'Fixation', scenario, save_heatmaps, HM_filename_fix)
        end

        % run a sweep of varying cidal drug death rate values
        if deathrate_sweep
            [dG_used_est, dG_used_fix, dN_used_est, dN_used_fix, t_est_list, t_fix_list] = population_times_delta(model, scenario, relative_to_N0, t_end, dt,  ...
                                                                                                            S_i, N_i, G1_i, G2_i, G1_susceptible, ... 
                                                                                                            G2_susceptible, kN_dsweep, kG_dsweep, linspace_n);
            if save_heatmaps
               directory = save_to;
               if G1_susceptible
                   G1_susceptible_str = '_G1susceptible';
               else
                   G1_susceptible_str = '';
               end

               if G2_susceptible
                   G2_susceptible_str = '_G2susceptible';
               else
                   G2_susceptible_str = '';
               end

               HM_filename_est = sprintf("%smodel%d_scenario%d_%s%s_Ni%0.1e_kN%0.4f_kG%0.4f_deathrate_establishment_heatmap",directory,model,scenario,...
                                     G1_susceptible_str,G2_susceptible_str, N_i, kN_dsweep, kG_dsweep);
               HM_filename_fix = sprintf("%smodel%d_scenario%d_%s%s_Ni%0.1e_kN%0.4f_kG%0.4f_deathrate_fixation_heatmap",directory,model,scenario,...
                                     G1_susceptible_str,G2_susceptible_str, N_i, kN_dsweep, kG_dsweep);
            else
                HM_filename_est = 'n/a';
                HM_filename_fix = 'n/a';
            end            
            plot_population_times_delta(dG_used_est, dN_used_est, t_est_list, 'Establishment', kN_dsweep, kG_dsweep, save_heatmaps, HM_filename_est)
            plot_population_times_delta(dG_used_fix, dN_used_fix, t_fix_list, 'Fixation', kN_dsweep, kG_dsweep, save_heatmaps, HM_filename_fix)
        end
    end
end

%% functions
function [time, S, N, G1, G2, dN_dt, t_est, t_fix] = simulate(model,scenario,relative_to_N0,t_end,dt,S_i,N_i,G1_i,G2_i)
    % This fuction runs the simulations for specified simulation parameters
    % and numerical values.
    % Input: 
    % [int or float] simulation time (t_end), [int or float] time step (dt), [int or float] initial population concentrations,
    % [int] simulation scenario
    % Output:
    % time, population concentrations as functions of time, conservation dN/dt as a
    % function of time, establishment time (t_est) and fixation time (t_fix) for the given scenario
    
    % bring globals into function space
    global N0 k n kS kN kG1 kG2 rSN rNS rG1S rG1N rG2G1 dS dN dG1 dG2
    
    G2 = 'n/a';
    
    % fluctuate the environment, run sims, and combine results (only for model 2, scenario 1)
    % distinquish between drug and no-drug environments
    fluc_interval = 12; %hours
    kS_drug = kS; kN_drug = kN; kG1_drug = kG1; rG1S_drug = rG1S; dS_drug = dS; dN_drug = dN; dG1_drug = dG1;
    kS_no = 0.3466; kN_no = 0.1733; kG1_no = 0; rG1S_no = rG1N; dS_no = 1/156; dN_no = 1/156; dG1_no = 1/156;


    t_fluc = 0:fluc_interval:t_end;
    S_i_fluc = S_i; N_i_fluc = N_i; G1_i_fluc = G1_i;
    
    for i = 1:(length(t_fluc)-1)
        % start with drug, remove drug on even intervals
        if rem(i,2)==0
            kS = kS_no; kN = kN_no; kG1 = kG1_no; rG1S = rG1S_no; dS = dS_no; dN = dN_no; dG1 = dG1_no;
        else
            kS = kS_drug; kN = kN_drug; kG1 = kG1_drug; rG1S = rG1S_drug; dS = dS_drug; dN = dN_drug; dG1 = dG1_drug;
        end

        t_i = t_fluc(i);
        t_end_fluc = t_fluc(i+1);

        [t, X] = ODE_solver_model2_scenario1(t_i,t_end_fluc,dt,S_i_fluc,N_i_fluc,G1_i_fluc);
        
        % calculate dT/dt
        dN_dt_temp = [];
        for j = 1:length(X(:,1))
            dN_dt_calc = X(j,1)*kS*(k^n)/(k^n + (X(j,1)+X(j,2)+X(j,3))^n) ...
                         + X(j,2)*kN*(k^n)/(k^n + (X(j,1)+X(j,2)+X(j,3))^n) ...
                         + X(j,3)*kG1*(k^n)/(k^n + (X(j,1)+X(j,2)+X(j,3))^n) ...
                         - X(j,1)*dS - X(j,2)*dN - X(j,3)*dG1;

            dN_dt_temp(j,:) = dN_dt_calc;
        end
        
        % update lists
        if i == 1
            time = t;
            S = X(:,1);
            N = X(:,2);
            G1 = X(:,3);
            dN_dt = dN_dt_temp;
        else
            time = [time; t];
            S = [S; X(:,1)];
            N = [N; X(:,2)];
            G1 = [G1; X(:,3)];
            dN_dt = [dN_dt; dN_dt_temp];
        end
        
        S_i_fluc = S(end);
        N_i_fluc = N(end);
        G1_i_fluc = G1(end);
    end
    
    % reset for safety
    kS = kS_drug; kN = kN_drug; kG1 = kG1_drug; rG1S = rG1S_drug; dS = dS_drug; dN = dN_drug; dG1 = dG1_drug;
    

    % calculate conservation equation results based on scenario and model
    
%     for i = 1:length(S)
%         if model == 1
%             if scenario == 2
%                 dN_dt_calc = S(i)*kS + N(i)*kN + G1(i)*kG1 +G2(i)*kG2 - S(i)*dS - N(i)*dN - G1(i)*dG1 - G2(i)*dG2;
%             elseif scenario == 1
%                 dN_dt_calc = S(i)*kS + N(i)*kN + G1(i)*kG1 - S(i)*dS - N(i)*dN - G1(i)*dG1 ;
%             end            
%         elseif model == 2
%             if scenario == 2
%                 dN_dt_calc =  S(i)*kS*(k^n)/(k^n + (S(i)+N(i)+G1(i)+G2(i))^n) ...
%                     + N(i)*kN*(k^n)/(k^n + (S(i)+N(i)+G1(i)+G2(i))^n) ...
%                     + G1(i)*kG1*(k^n)/(k^n + (S(i)+N(i)+G1(i)+G2(i))^n) ...
%                     + G2(i)*kG2*(k^n)/(k^n + (S(i)+N(i)+G1(i)+G2(i))^n) ...
%                     - S(i)*dS - N(i)*dN - G1(i)*dG1 - G2(i)*dG2;
%             elseif scenario == 1
%                 dN_dt_calc = S(i)*kS*(k^n)/(k^n + (S(i)+N(i)+G1(i))^n) ...
%                     + N(i)*kN*(k^n)/(k^n + (S(i)+N(i)+G1(i))^n) ...
%                     + G1(i)*kG1*(k^n)/(k^n + (S(i)+N(i)+G1(i))^n) ...
%                     - S(i)*dS - N(i)*dN - G1(i)*dG1;
%             end
%         end
%         dN_dt_list(i,:) = dN_dt_calc;
%     end
%     dN_dt = dN_dt_list;

    % establishment and fixation time calculation
    establishment = false;
    fixation = false;
    t_est = "> t_end";
    t_fix = "> t_end";
    index = 1;
    pop_fraction = 0;
    if scenario == 2
        population = G2;
    elseif scenario == 1
        population = G1;
    end
    % calculate relative to either N_0 or N_tot
    if relative_to_N0
        % loop through the results until fixation is found
        while fixation == false && index <= length(population)
            pop_t = time(index);
            pop_fraction = population(index)/N0;
            % check for establishment
            if establishment == false && pop_fraction > 0.05
                t_est = pop_t;
                establishment = true;
            end
            % check for fixation (exit loop if found)
            if pop_fraction > 0.95
                t_fix = pop_t;
                fixation = true;
            end            
            index = index + 1;
        end
    else
        while fixation == false && index <= length(population)
            pop_t = time(index);
            if scenario == 2
                pop_fraction = population(index) / (S(index) + N(index) + G1(index) + population(index));
            elseif scenario == 1
                pop_fraction = population(index) / (S(index) + N(index) + population(index));
            end
            % check for establishment
            if establishment == false && pop_fraction > 0.05
                t_est = pop_t;
                establishment = true; % to make sure t_est isn't overwritten
            end
            % check for fixation (exits loop if found)
            if pop_fraction > 0.95
                t_fix = pop_t;
                fixation = true;
            end   
            index = index + 1;
        end
    end
    fprintf("t_est = %0.2f, t_fix = %0.2f\n", t_est, t_fix) % print est/fix time results to command window
end


function [kG_used_est, kG_used_fix, kN_used_est, kN_used_fix, t_est_list, t_fix_list] = population_times(model, scenario, relative_to_N0, t_end, dt, S_i, N_i, G1_i, G2_i, ...
                                                                                                         cidal, G1_susceptible, G2_susceptible, death_factor, linspace_n)
    % Create simulations for large variations of kG1 and kN
    % Input:
    % [int] number of combinations to consider (linspace_n),
    % all parameters for simulate() function, [int or float] cidal drug
    % death factor (death_factor)
    % Output:
    % lists of the kG and kN values used for both the establishment and
    % fixation times, lists for the corresponding establishment and
    % fixation time results
    
    % bring globals into function space

    global N0 k n kS kN kG1 kG2 rSN rNS rG1S rG1N rG2G1 dS dN dG1 dG2
    % create intial lists of kG and kN values and create array of all
    % combinations of these values
    kG_list = linspace(0.1733,0.3466,linspace_n);
    kN_list = linspace(0.1733,0.3466,linspace_n);
    kG_kN_combinations = combvec(kG_list,kN_list);

    kG_used_est = [];
    kG_used_fix = [];
    kN_used_est = [];
    kN_used_fix = [];
    t_est_list = [];
    t_fix_list = [];
    
    % Go through all combinations of kG and kN and run simulations
    for i = 1:length(kG_kN_combinations)
        combination = kG_kN_combinations(:,i); % gets a combination
        kG = combination(1);
        kN = combination(2);
        
        % change values based on scenario
        if scenario == 1
            kG1 = kG;
        elseif scenario == 2
            kG2 = kG;
            kG1 = 0.1733;            
        end
        
        % use death factor if one is given
        if cidal
            dN = 0.1733*death_factor;
            dS = 0.3466*death_factor;
            if G1_susceptible
                dG1 = 0.2600*death_factor;
            else
                dG1 = 1/156;
            end
            if G2_susceptible
                dG2 = 0.2600*death_factor;
            else
                dG2 = 0.0;
            end
        else
            dS = 1/156; dN = 1/156; dG1 = 1/156; dG2 = 0.0;
        end
        
        % check if kN <= kG (as required)
        if kN <= kG
            [~,~,~,~,~,~, t_est, t_fix] = simulate(model,scenario,relative_to_N0,t_end, dt, S_i, N_i, G1_i, G2_i); % run simulation
            % if t_est is found, save it and the kG and kN used to the result lists
            if isnumeric(t_est)
                kG_used_est = [kG_used_est,kG];
                kN_used_est = [kN_used_est,kN];
                t_est_list = [t_est_list,t_est];
            end
            % same for t_fix
            if isnumeric(t_fix)
                kG_used_fix = [kG_used_fix,kG];
                kN_used_fix = [kN_used_fix,kN];
                t_fix_list = [t_fix_list,t_fix];
            end
        end
    end
end

function [dG_used_est, dG_used_fix, dN_used_est, dN_used_fix, t_est_list, t_fix_list] = population_times_delta(model, scenario, relative_to_N0, t_end, dt, S_i, N_i, G1_i, G2_i, ...
                                                                                                               G1_susceptible, G2_susceptible, kN_dsweep, kG_dsweep, linspace_n)
    % Create simulations for large variations of kG1 and kN
    % Input:
    % [int] number of combinations to consider (linspace_n),
    % all parameters for simulate() function, [int or float] cidal drug
    % death factor (death_factor)
    % Output:
    % lists of the kG and kN values used for both the establishment and
    % fixation times, lists for the corresponding establishment and
    % fixation time results
    
    % bring globals into function space

    global N0 k n kS kN kG1 kG2 rSN rNS rG1S rG1N rG2G1 dS dN dG1 dG2
    % create intial lists of kG and kN values and create array of all
    % combinations of these values

    kN = kN_dsweep;
    if scenario == 1
        kG1 = kG_dsweep;
    elseif scenario == 2
        kG1 = 0.1733;
        kG2 = kG_dsweep;
    end

    dN_list = linspace(0.1,1.0,linspace_n);
    dG_list = linspace(1/156,0.05,linspace_n);
    dG_dN_combinations = combvec(dG_list,dN_list);

    dG_used_est = [];
    dG_used_fix = [];
    dN_used_est = [];
    dN_used_fix = [];
    t_est_list = [];
    t_fix_list = [];
    
    % Go through all combinations of dG and dN and run simulations
    for i = 1:length(dG_dN_combinations)
        combination = dG_dN_combinations(:,i); % gets a combination
        dG1 = combination(1);
        dN = combination(2);
        dS = 1.0;

        [~,~,~,~,~,~, t_est, t_fix] = simulate(model,scenario,relative_to_N0,t_end, dt, S_i, N_i, G1_i, G2_i); % run simulation
        % if t_est is found, save it and the kG and kN used to the result lists
        if isnumeric(t_est)
            dG_used_est = [dG_used_est,dG1];
            dN_used_est = [dN_used_est,dN];
            t_est_list = [t_est_list,t_est];
        end
        % same for t_fix
        if isnumeric(t_fix)
            dG_used_fix = [dG_used_fix,dG1];
            dN_used_fix = [dN_used_fix,dN];
            t_fix_list = [t_fix_list,t_fix];
        end
    end
end

function plot_simulation_results(t,S,N,G1,G2,dN_dt,legend_set,log_bool, savefig_bool, file_name)
    % Plots the results of a given simulation
    % Inputs:
    % [list] 1D lists containing results for different simulations, [set]
    % legend string for different kN kG values (legend_set), bool to use
    % loglog instead of regular plots (log_bool), bool to save figure to a
    % file (savefig_bool), [str] name of file with or without extension
    % (default is .fig)
    
    % Create figure and subplots, for each population plot all results on the same subplot
    fig = figure; 
    figure('DefaultAxesFontSize',28);
    
    subplot(2,3,1);
    % plot each of the simulation results in the input lists
    for i = 1:length(S(1,:))
        % plot in loglog 
        if log_bool
            % for different simulation sets, use different line formats
            if i == 1
                loglog(t(:,i),S(:,i),'LineWidth',4);
            elseif i == 2
                loglog(t(:,i),S(:,i), '-.','LineWidth',4);
            elseif i == 3
                loglog(t(:,i),S(:,i),'--','LineWidth',4);
            end
        else
            if i == 1
                plot(t(:,i),S(:,i),'LineWidth',4);
            elseif i == 2
                plot(t(:,i),S(:,i), '-.','LineWidth',4);
            elseif i == 3
                plot(t(:,i),S(:,i),'--','LineWidth',4);
            end
        end
        hold on
    end
    hold off
    % axis labels and title
    xlabel('t (hr)'); ylabel('S(t) (cells/mL)');
    
    xlim([10^-2 10^3]);
    ylim([10^-5 10^9]);
    xticks([ 10^-2 10^0 10^3])
    xticklabels({'10^{-2}','10^{0}','10^{3}'})
    
    ttl = title('A)');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.36; 
    ttl.HorizontalAlignment = 'left';
    
%     subplot(2,3,2);
% %     %plot(0,0,  0,0,  0,0, 'LineWidth', 2);
% %     plot(0,0,'LineWidth',3);
% %     hold on
% %     plot(0,0, '-.','LineWidth',3);
% %     hold on
% %     plot(0,0,'--','LineWidth',3);
% %     hold off
% %     axis off;
%     ttl = title('A)');
%     ttl.Units = 'Normalize'; 
%     ttl.Position(1) = -0.33; 
%     ttl.HorizontalAlignment = 'left';
    % plot legend 
    l = legend(legend_set, 'Location', 'southwest');
    l.FontSize = 18;
    
    subplot(2,3,2);
    % same comments as above
    for i = 1:length(N(1,:))
        if log_bool
            if i == 1
                loglog(t(:,i),N(:,i),'LineWidth',4);
            elseif i == 2
                loglog(t(:,i),N(:,i), '-.','LineWidth',4);
            elseif i == 3
                loglog(t(:,i),N(:,i),'--','LineWidth',4);
            end
        else
            if i == 1
                plot(t(:,i),N(:,i),'LineWidth',4);
            elseif i == 2
                plot(t(:,i),N(:,i), '-.','LineWidth',4);
            elseif i == 3
                plot(t(:,i),N(:,i),'--','LineWidth',4);
            end
        end
        hold on
    end
    hold off
    xlabel('t (hr)'); ylabel('N(t) (cells/mL)');
    
    xlim([10^-2 10^3]);
    ylim([10^-5 10^9]);

    xticks([ 10^-2 10^0 10^3])
    xticklabels({'10^{-2}','10^{0}','10^{3}'})
    
    ttl = title('B)');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.35; 
    ttl.HorizontalAlignment = 'left';
    %title('N')
    
    subplot(2,3,3);
    for i = 1:length(G1(1,:))
        if log_bool
            if i == 1
                loglog(t(:,i),G1(:,i),'LineWidth',4);
            elseif i == 2
                loglog(t(:,i),G1(:,i), '-.','LineWidth',4);
            elseif i == 3
                loglog(t(:,i),G1(:,i),'--','LineWidth',4);
            end
        else
            if i == 1
                plot(t(:,i),G1(:,i),'LineWidth',4);
            elseif i == 2
                plot(t(:,i),G1(:,i), '-.','LineWidth',4);
            elseif i == 3
                plot(t(:,i),G1(:,i),'--','LineWidth',4);
            end
        end
        hold on
    end
    hold off
    xlabel('t (hr)'); ylabel('G(t) (cells/mL)');
    
    xlim([10^-2 10^3]);
    ylim([10^-5 10^9]);
% 
    xticks([ 10^-2 10^0 10^3])
    xticklabels({'10^{-2}','10^{0}','10^{3}'})
    
    ttl = title('C)');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.35; 
    ttl.HorizontalAlignment = 'left';
    %title('G')
    
    if isnumeric(G2)
        subplot(2,3,5);
        for i = 1:length(G2(1,:))
            if log_bool
                if i == 1
                    loglog(t(:,i),G2(:,i),'LineWidth',4);
                elseif i == 2
                    loglog(t(:,i),G2(:,i), '-.','LineWidth',4);
                elseif i == 3
                    loglog(t(:,i),G2(:,i),'--','LineWidth',4);
                end
            else
                if i == 1
                    plot(t(:,i),G2(:,i),'LineWidth',4);
                elseif i == 2
                    plot(t(:,i),G2(:,i), '-.','LineWidth',4);
                elseif i == 3
                    plot(t(:,i),G2(:,i),'--','LineWidth',4);
                end
            end
            hold on
        end
        hold off
        xlabel('t (hr)'); ylabel('G_2(t) (cells/mL)');
        %title('G_2') 
    end
    
    subplot(2,3,4);
    for i = 1:length(G1(1,:))
        T = S(:,i) + N(:,i) + G1(:,i);
        if log_bool
            if i == 1
                loglog(t(:,i),T,'LineWidth',4);
            elseif i == 2
                loglog(t(:,i),T, '-.','LineWidth',4);
            elseif i == 3
                loglog(t(:,i),T,'--','LineWidth',4);
            end
        else
            if i == 1
                plot(t(:,i),T,'LineWidth',4);
            elseif i == 2
                plot(t(:,i),T, '-.','LineWidth',4);
            elseif i == 3
                plot(t(:,i),T,'--','LineWidth',4);
            end
        end
        hold on
    end
    hold off
    xlabel('t (hr)'); ylabel('T(t) (cells/mL)');
    
    xlim([10^-2 10^3]);
    ylim([10^-5 10^9]);

    xticks([ 10^-2 10^0 10^3])
    xticklabels({'10^{-2}','10^{0}','10^{3}'})
    
    ttl = title('D)');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.35; 
    ttl.HorizontalAlignment = 'left';
    
    subplot(2,3,5);
    for i = 1:length(G1(1,:))
        if isnumeric(G2)
            ratio = G2(:,i)./N(:,i);
            ratio_string = 'G_2(t)/N(t)';
            ratio_title_string = 'G_2/N';
        else
            ratio = G1(:,i)./(S(:,i) + N(:,i) + G1(:,i));
            ratio_string = 'G(t)/T(t)';
            ratio_title_string = 'G/T Ratio';

        end
        if log_bool
            if i == 1
                loglog(t(:,i),ratio,'LineWidth',4);
            elseif i == 2
                loglog(t(:,i),ratio, '-.','LineWidth',4);
            elseif i == 3
                loglog(t(:,i),ratio,'--','LineWidth',4);
            end
        else
            if i == 1
                plot(t(:,i),ratio,'LineWidth',4);
            elseif i == 2
                plot(t(:,i),ratio, '-.','LineWidth',4);
            elseif i == 3
                plot(t(:,i),ratio,'--','LineWidth',4);
            end
        end
        hold on
    end
    hold off
    xlabel('t (hr)'); ylabel(ratio_string);
    
    xlim([10^-2 10^3]);
    ylim([10^-10 10^1]);

    xticks([ 10^-2 10^0 10^3])
    xticklabels({'10^{-2}','10^{0}','10^{3}'})
    
    
    ttl = title('E)');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.39; 
    ttl.HorizontalAlignment = 'left';
    %title(ratio_title_string)

    subplot(2,3,6);
    for i = 1:length(dN_dt(1,:))
        if i == 1
            semilogx(t(:,i),dN_dt(:,i),'LineWidth',2);
        elseif i == 2
            semilogx(t(:,i),dN_dt(:,i), '-.','LineWidth',2);
        elseif i == 3
            semilogx(t(:,i),dN_dt(:,i),'--','LineWidth',2);
        end
%         if log_bool 
%             if i == 1
%                 loglog(t(:,i),dN_dt(:,i),'LineWidth',2);
%             elseif i == 2
%                 loglog(t(:,i),dN_dt(:,i), '-.','LineWidth',2);
%             elseif i == 3
%                 loglog(t(:,i),dN_dt(:,i),'--','LineWidth',2);
%             end
%         else
%             if i == 1
%                 plot(t(:,i),dN_dt(:,i),'LineWidth',2);
%             elseif i == 2
%                 plot(t(:,i),dN_dt(:,i), '-.','LineWidth',2);
%             elseif i == 3
%                 plot(t(:,i),dN_dt(:,i),'--','LineWidth',2);
%             end
%         end
        hold on
    end
    hold off
    xlabel('t (hr)'); ylabel('dT/dt (cells/mLhr)');
  
    xlim([10^-2 10^3]);
%     
% 
    xticks([ 10^-2 10^0 10^3])
    xticklabels({'10^{-2}','10^{0}','10^{3}'})

    ttl = title('F)');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.35; 
    ttl.HorizontalAlignment = 'left';
    %title('dT/dt')
    
    
    set(gcf, 'PaperUnits', 'inches');
    x_width=18 ;y_width=12;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]);  
    orient(fig,'landscape');
    if savefig_bool
        saveas(gcf, append(file_name, '.png'))
        %save(append(file_name, '.mat'),'t','S','N','G1','G2','dN_dt')
    end
end


function plot_population_times(kG_list, kN_list, population_time_list, population_type, scenario, savefig_bool, file_name)
    % Plot heatmaps using the kG, kN and t_est or t_fix results for a large
    % amount of simulations.
    % Input:
    % [list] lists of the kG, kN used and t_est or t_fix results, [str]
    % quantitative measure used ("Establishment" or "Fixation") for plot
    % titles (population_type), [int] scenario simulated for plot titles
    % and file names(scenario), [bool] to save figure or not
    % (savefig_bool),[str] name of file to save (with or without file type extension, default is .fig if there is none)
    
    % create figure
    fig = figure;
    
    % put input lists into a table to make the heatmaps with
    X = kG_list(:); Y = kN_list(:); Z = population_time_list(:);
    tbl = table(X,Y,Z);
    
    % create the heatmap 
    hHM = heatmap(tbl,'X', 'Y', 'ColorVariable','Z','ColorMethod', 'none', 'CellLabelColor','none',...
                  'GridVisible','off','MissingDataLabel', 'No Data', 'MissingDataColor','w', ...
                  'FontSize', 20);%,'ColorScaling','log');
    hHM.NodeChildren(3).YDir='normal'; % flips the y-axis to make it increasing (matlab default is decreasing)
    hHM.Colormap = jet;
    
    
    % set axis labels and title
    ylabel('k_N (/hr)')
    if scenario == 2
        xlabel('k_{G2} (/hr)')
        title('')%append(population_type,' Time as a Function of k_N and k_{G2} (k_N <= k_{G2})'))
    else
        xlabel('k_{G} (/hr)')
        title('')%append(population_type,' Time as a Function of k_N and k_{G} (k_N <= k_{G})'))        
    end
    colorbar; % adds a color bar for the heatmap

    set(gcf, 'PaperUnits', 'inches');
    x_width=10 ;y_width=7;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
    
    if population_type == "Establishment"
        annotation('textarrow',[0.98,0.98],[0.7,0.7],'string','{\it \tau_{est}} (hr)', 'FontSize', 22, ...
            'HeadStyle','none','LineStyle','none','HorizontalAlignment','right','TextRotation',90);
        annotation('textarrow',[0.025,0.025],[0.97,0.97],'string','A)', 'FontWeight', 'bold', ...
            'FontSize', 22, 'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');
    else
        annotation('textarrow',[0.98,0.98],[0.7,0.7],'string','{\it \tau_{fix}} (hr)', 'FontSize', 22, ...
            'HeadStyle','none','LineStyle','none','HorizontalAlignment','right','TextRotation',90);
        annotation('textarrow',[0.025,0.025],[0.97,0.97],'string','A)', 'FontWeight', 'bold', ...
            'FontSize', 22, 'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');
    end
    
    % rotate x-axis labels
    s = struct(hHM);
    s.XAxis.TickLabelRotation = 60;  % angled
    s.XDisplayLabels = compose('%.3f',str2double(hHM.XDisplayLabels));
    s.YDisplayLabels = compose('%.4f',str2double(hHM.YDisplayLabels));
    
    %s.FontSize = 60;
    % save the file 
    if savefig_bool
        saveas(gcf, append(file_name,'.png'))
        %save(append(file_name, '.mat'), 'kG_list', 'kN_list', 'population_time_list')
    end
    
    close(fig)
end

function plot_population_times_delta(dG_list, dN_list, population_time_list, population_type, kN_dsweep, kG_dsweep, savefig_bool, file_name)
    % create figure
    fig = figure;
    
    % put input lists into a table to make the heatmaps with
    X = dG_list(:); Y = dN_list(:); Z = population_time_list(:);
    tbl = table(X,Y,Z);
    
    % create the heatmap 
    hHM = heatmap(tbl,'X', 'Y', 'ColorVariable','Z','ColorMethod', 'none', 'CellLabelColor','none',...
                  'GridVisible','off','MissingDataLabel', 'No Data', 'MissingDataColor','w', ...
                  'FontSize', 20);
    hHM.NodeChildren(3).YDir='normal'; % flips the y-axis to make it increasing (matlab default is decreasing)
    hHM.Colormap = jet;
    
    
    % set axis labels and title
    ylabel('\delta_N (/hr)')
    xlabel('\delta_G (/hr)')
    title('')       

    colorbar; % adds a color bar for the heatmap
    
    set(gcf, 'PaperUnits', 'inches');
    x_width=10 ;y_width=7.2;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
    
    if population_type == "Establishment"
        annotation('textarrow',[0.98,0.98],[0.7,0.7],'string','{\it \tau_{est}} (hr)', 'FontSize', 22, ...
            'HeadStyle','none','LineStyle','none','HorizontalAlignment','right','TextRotation',90);
        annotation('textarrow',[0.025,0.025],[0.97,0.97],'string','A)', 'FontWeight', 'bold', ...
            'FontSize', 22, 'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');
    else
        annotation('textarrow',[0.983,0.983],[0.7,0.7],'string','{\it \tau_{fix}} (hr)', 'FontSize', 22, ...
            'HeadStyle','none','LineStyle','none','HorizontalAlignment','right','TextRotation',90);
        annotation('textarrow',[0.025,0.025],[0.97,0.97],'string','B)', 'FontWeight', 'bold', ...
            'FontSize', 22, 'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');
    end
    
    % rotate x-axis labels
    s = struct(hHM);
    s.XAxis.TickLabelRotation = 60;  % angled
    s.XDisplayLabels = compose('%.3f',str2double(hHM.XDisplayLabels));
    s.YDisplayLabels = compose('%.4f',str2double(hHM.YDisplayLabels));
    
 
    
    %s.FontSize = 60;
    % save the file 
    if savefig_bool
        saveas(gcf, append(file_name,'.png'))
        %save(append(file_name, '.mat'), 'kG_list', 'kN_list', 'population_time_list')
    end
    
    close(fig)
end
