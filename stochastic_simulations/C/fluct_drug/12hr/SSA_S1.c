// Joshua Guthrie, Charlebois Laboratory, University of Alberta

/* REACTIONS
   ********************************************************************
   0: S ----> 2S      (parameter = k_S*z(N_tot))
   1: N ----> 2N      (parameter = k_N*z(N_tot))
   2: G1 ----> 2G1    (parameter = k_G1*z(N_tot))
   3: S ----> N       (parameter = r_NS)
   4: N ----> S       (parameter = r_SN)
   5: S ----> G1      (parameter = r_G1S) 
   6: N ----> G1      (parameter = r_G1N) 
   7: S ----> 0       (parameter = delta_S)
   8: N ----> 0      (parameter = delta_N)
   9: G1 ----> 0     (parameter = delta_G1) 
   ********************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[]){
    /* Command line arguments: 
       argv[1] = k_N [float]
       argv[2] = k_G1 [float]
       argv[3] = delta_S [float]
       argv[4] = delta_N [float]
       argv[5] = N_i [int]
       argv[6] = t_end [float]
       argv[7] = num_runs [int]
       argv[8] = distributions outfile name [str]
       argv[9] = number for sample trajectory to save (0 for no save) [int]
       argv[10] = sample trajectory outfile name [str]
       argv[11] = stop at G1 appearance, establishment, or fixation (1, 2, or 3, use 0 otherwise) [int]
       argv[12] = measure relative to initial or total population (0 for initial (default), 1 for total) [int]
    */
    
    int sample_traj = atoi(argv[9]);

    int rel_to = atoi(argv[12]);

    FILE *outfile_est_fix = fopen(argv[8], "w");
    FILE *outfile = NULL;
    if (sample_traj > 0) {
        outfile = fopen(argv[10], "w");
    }
    

    /* INITIALIZATION */

    // Barayni-Hill coefficients
    double k = 1e+7;
    double n = 2;

    // reaction parameters c_i (monomolecular reactions, k_i = c_i)
    float k_S = 0.0;
    float k_N = atof(argv[1]);
    float k_G1 = atof(argv[2]); 
    float r_NS = 0.0625;
    float r_SN = 0.0035;
    float r_G1S = 0.0;
    float r_G1N = (1e-6)/3;
    float delta_S = atof(argv[3]);
    float delta_N = atof(argv[4]);
    float delta_G1 = 1/156; 

    // state changes for each reaction 
    int V[10][3] = {{+1, 0, 0},
                    {0, +1, 0},
                    {0, 0, +1},
                    {-1, +1, 0},
                    {+1, -1, 0},
                    {-1, 0, +1},
                    {0, -1, +1},
                    {-1, 0, 0},
                    {0, -1, 0},
                    {0, 0, -1}};
    
    
    float t_end = atof(argv[6]); // time to simulate, hours
    int num_runs = atoi(argv[7]); // Number of runs (trajectories)

    // initialize state (number of cells for each population) and time
    long int S = 550000;
    long int N = atoi(argv[5]);
    long int G1 = 0;
    double t = 0.0; // hours

    // conservation equations
    long int N_tot = S + N + G1;
    double hill = pow(k,n)/(pow(k,n) + pow(N_tot,n));
    double dNtot_dt = k_S*hill*S + k_N*hill*N + k_G1*hill*G1 - delta_S*S - delta_N*N - delta_G1*G1;
    
    // Print header for appearance, establishment, fixation time output file
    fprintf(outfile_est_fix, "#Stochastic Simulation Algorithm Results (Establishment/Fixation Times)\n#\n");

    fprintf(outfile_est_fix, "#Reaction Parameters (/hr)\n#-----------------------------------------------------------------------\n");
    fprintf(outfile_est_fix, "#Growth Rates: k_S = %.4f, k_N = %.4f, k_G1 = %.4f\n", k_S, k_N, k_G1);
    fprintf(outfile_est_fix, "#Switching Rates: r_NS = %.4f, r_SN = %.4f, r_G1S = %.4f, r_G1N = %.4e\n", r_NS, r_SN, r_G1S, r_G1N);
    fprintf(outfile_est_fix, "#Death Rates: d_S = %0.1f, d_N = %0.1f, d_G1 = %0.1f\n", delta_S, delta_N, delta_G1);
    fprintf(outfile_est_fix, "#-----------------------------------------------------------------------\n#\n");

    fprintf(outfile_est_fix, "#Simulation Parameters\n#-----------------------------------------------------------------------\n");
    fprintf(outfile_est_fix, "#Barayni-Hill: k = %.0e, n = %.0f\n", k, n);
    fprintf(outfile_est_fix, "#Time (hr): t_i = %0.1f, t_end = %0.1f\n", t, t_end);
    fprintf(outfile_est_fix, "#Initial State: S = %ld, N = %ld, G1 = %ld, N_tot = %ld, dN_tot/dt = %0.10f\n", S, N, G1, N_tot, dNtot_dt);
    fprintf(outfile_est_fix, "#-----------------------------------------------------------------------\n#\n");
    fprintf(outfile_est_fix, "#ESTABLISHMENT/FIXATION TIMES (if time is negative it wasn't found within t_end)\n");
    fprintf(outfile_est_fix, "#Measured relative to: %d\n", rel_to);
    fprintf(outfile_est_fix, "#Columns: Run (Trajectory) Number, G1 t_appearance, G1 t_est, G1 t_fix\n#-----------------------------------------------------------------------\n");
       
    if ((outfile != NULL)&&(sample_traj > 0)){
        // Print header for data file 
        fprintf(outfile, "#Stochastic Simulation Algorithm Results\n#\n");

        fprintf(outfile, "#Reaction Parameters (/hr)\n#-----------------------------------------------------------------------\n");
        fprintf(outfile, "#Growth Rates: k_S = %.4f, k_N = %.4f, k_G1 = %.4f\n", k_S, k_N, k_G1);
        fprintf(outfile, "#Switching Rates: r_NS = %.4f, r_SN = %.4f, r_G1S = %.4f, r_G1N = %.4e\n", r_NS, r_SN, r_G1S, r_G1N);
        fprintf(outfile, "#Death Rates: d_S = %0.1f, d_N = %0.1f, d_G1 = %0.1f\n", delta_S, delta_N, delta_G1);
        fprintf(outfile, "#-----------------------------------------------------------------------\n#\n");

        fprintf(outfile, "#Simulation Parameters\n#-----------------------------------------------------------------------\n");
        fprintf(outfile, "#Barayni-Hill: k = %.0e, n = %.0f\n", k, n);
        fprintf(outfile, "#Time (hr): t_i = %0.1f, t_end = %0.1f\n", t, t_end);
        fprintf(outfile, "#Initial State: S = %ld, N = %ld, G1 = %ld, N_tot = %ld, dN_tot/dt = %0.10f\n", S, N, G1, N_tot, dNtot_dt);
        fprintf(outfile, "#-----------------------------------------------------------------------\n#\n");

        fprintf(outfile, "#DATA (Columns: Row Number, Time, S, N, G1, N_tot, dN_tot/dt)\n#-----------------------------------------------------------------------\n");
        fprintf(outfile, "%d    %0.20f    %ld    %ld    %ld    %ld    %0.20f\n", 0, t, S, N, G1, N_tot, dNtot_dt);
    }

    
    srand(time(NULL)); // initialize random number generator

    int stop_at = atoi(argv[11]);
    for (int i = 1; i <= num_runs; i++){
        printf("Trajectory %d running...\n", i);

        // initialize state (number of cells for each population) and time
        long int S = 550000;
        long int N = atoi(argv[5]);
        long int G1 = 0;

        double t = 0.0; // hours
        long int N_0 = S + N + G1;
        
        // conservation equations
        long int N_tot = S + N + G1;

        double hill = pow(k,n)/(pow(k,n) + pow(N_tot,n));
        double dNtot_dt = k_S*hill*S + k_N*hill*N + k_G1*hill*G1 - delta_S*S - delta_N*N - delta_G1*G1;

        // loop control bools
        bool G1_app_found = false; 
        bool G1_est_found = false;
        bool G1_fix_found = false;


        // initialize establishment and fixation times (-1.0 is a dummy value)
        double G1_app = -1.0;
        double G1_est = -1.0;
        double G1_fix = -1.0;
        
        float G1_frac = -1.0;


        if (i == sample_traj){
            long int row_num = 1;
            while (t < t_end){
                /* Runs until t_end */
                // fluctuate the drug
                if (((t > 12.0)&&(t <= 24.0))||((t > 36.0)&&(t <= 48.0))||((t > 60.0)&&(t <= 72.0))||((t > 84.0)&&(t <= 96.0))
                    ||((t > 108.0)&&(t <= 120.0))||((t > 132.0)&&(t <= 144.0))||((t > 156.0)&&(t <= 168.0))||((t > 180.0)&&(t <= 192.0))
                    ||((t > 204.0)&&(t <= 216.0))||((t > 228.0)&&(t <= 240.0))){
                    // no drug condition
                    k_S = 0.3466;
                    k_N = 0.2600;
                    k_G1 = 0.1733; 
                    r_NS = 0.0625;
                    r_SN = 0.0035;
                    r_G1S = (1e-6)/3;
                    r_G1N = (1e-6)/3;
                    delta_S = 1/156;
                    delta_N = 1/156;
                    delta_G1 = 1/156;
                    }
                else {
                    // drug
                    k_S = 0.0;
                    k_N = atof(argv[1]);
                    k_G1 = atof(argv[2]); 
                    r_NS = 0.0625;
                    r_SN = 0.0035;
                    r_G1S = 0.0;
                    r_G1N = (1e-6)/3;
                    delta_S = atof(argv[3]);
                    delta_N = atof(argv[4]);
                    delta_G1 = 1/156;
                }
                
                // STEP 1: calculate propensity functions and their sum a_0
                double a[10] = {(k_S*hill)*S, (k_N*hill)*N, (k_G1*hill)*G1,
                                r_NS*S, r_SN*N, r_G1S*S, r_G1N*N, 
                                delta_S*S, delta_N*N, delta_G1*G1};
                
                double a_0 = 0.0;
                for (int i = 0; i < 10; i++) a_0 += a[i];

                // STEP 2: generate uniform random numbers in unit interval, calculate tau and j
                double r1 = (double)rand() / (double)RAND_MAX;
                while (r1 == 0.0) r1 = (double)rand() / (double)RAND_MAX;

                double tau = (1/a_0)*log(1/r1);

                double r2 = (double)rand() / (double)RAND_MAX;
                int j = 0;
                double a_sum = a[j]; 
                while (a_sum < r2*a_0){
                    j += 1;
                    a_sum += a[j];
                }

                // STEP 3: Update
                t += tau;

                S += V[j][0];
                N += V[j][1];
                G1 += V[j][2];

                N_tot = S + N + G1;
                hill = pow(k,n)/(pow(k,n) + pow(N_tot,n));
                dNtot_dt = k_S*hill*S + k_N*hill*N + k_G1*hill*G1 - delta_S*S - delta_N*N - delta_G1*G1;

                // print out the sample trajectory results
                fprintf(outfile, "%ld    %0.20f    %ld    %ld    %ld    %ld    %0.20f\n", row_num, t, S, N, G1, N_tot, dNtot_dt);

                // check for establishment and fixation for the G populations

                if (rel_to == 1){
                    G1_frac = (float) G1/N_tot;
                }
                else {
                    G1_frac = (float) G1/N_0;
                }

                // Check G1 and G2 first appearance
                if ((G1_app_found == false) && (G1 > 0)){
                    G1_app = t;
                    G1_app_found = true;
                }

                // Check G1 t_est and t_fix
                if ((G1_est_found == false) && (G1_frac > 0.05)){
                    G1_est = t;
                    G1_est_found = true;
                }
                else if ((G1_fix_found == false) && (G1_frac > 0.95)){
                    G1_fix = t;
                    G1_fix_found = true;
                }

                row_num++;
            }
        }
        else {
            while (t < t_end){
                /* Runs until t_end OR G2 fixates OR until argv[11] parameter is found */
            
                // fluctuate the drug
                if (((t > 12.0)&&(t <= 24.0))||((t > 36.0)&&(t <= 48.0))||((t > 60.0)&&(t <= 72.0))||((t > 84.0)&&(t <= 96.0))
                    ||((t > 108.0)&&(t <= 120.0))||((t > 132.0)&&(t <= 144.0))||((t > 156.0)&&(t <= 168.0))||((t > 180.0)&&(t <= 192.0))
                    ||((t > 204.0)&&(t <= 216.0))||((t > 228.0)&&(t <= 240.0))){
                    // no drug condition
                    k_S = 0.3466;
                    k_N = 0.2600;
                    k_G1 = 0.1733; 
                    r_NS = 0.0625;
                    r_SN = 0.0035;
                    r_G1S = (1e-6)/3;
                    r_G1N = (1e-6)/3;
                    delta_S = 1/156;
                    delta_N = 1/156;
                    delta_G1 = 1/156;
                    }
                else {
                    // drug
                    k_S = 0.0;
                    k_N = atof(argv[1]);
                    k_G1 = atof(argv[2]); 
                    r_NS = 0.0625;
                    r_SN = 0.0035;
                    r_G1S = 0.0;
                    r_G1N = (1e-6)/3;
                    delta_S = atof(argv[3]);
                    delta_N = atof(argv[4]);
                    delta_G1 = 1/156;
                }

                // STEP 1: calculate propensity functions and their sum a_0
                double a[10] = {(k_S*hill)*S, (k_N*hill)*N, (k_G1*hill)*G1, 
                                r_NS*S, r_SN*N, r_G1S*S, r_G1N*N, 
                                delta_S*S, delta_N*N, delta_G1*G1};
                
                double a_0 = 0.0;
                for (int i = 0; i < 10; i++) a_0 += a[i];

                // STEP 2: generate uniform random numbers in unit interval, calculate tau and j
                double r1 = (double)rand() / (double)RAND_MAX;
                while (r1 == 0.0) r1 = (double)rand() / (double)RAND_MAX;

                double tau = (1/a_0)*log(1/r1);

                double r2 = (double)rand() / (double)RAND_MAX;
                int j = 0;
                double a_sum = a[j]; 
                while (a_sum < r2*a_0){
                    j += 1;
                    a_sum += a[j];
                }

                // STEP 3: Update
                t += tau;

                S += V[j][0];
                N += V[j][1];
                G1 += V[j][2];

                N_tot = S + N + G1;
                hill = pow(k,n)/(pow(k,n) + pow(N_tot,n));

                // check for establishment and fixation for the G populations
                if (rel_to == 1){
                    G1_frac = (float) G1/N_tot;
                }
                else {
                    G1_frac = (float) G1/N_0;
                }

                // Check G1 and G2 first appearance
                if ((G1_app_found == false) && (G1 > 0)){
                    G1_app = t;
                    G1_app_found = true;
                    if (stop_at == 1) break;
                }

                // Check G1 t_est and t_fix
                if ((G1_est_found == false) && (G1_frac > 0.05)){
                    G1_est = t;
                    G1_est_found = true;
                    if (stop_at == 2) break;
                }
                else if ((G1_fix_found == false) && (G1_frac > 0.95)){
                    G1_fix = t;
                    G1_fix_found = true;
                    break;
                }
            }
        }
        fprintf(outfile_est_fix, "%d    %0.20f    %0.20f    %0.20f\n", i, G1_app, G1_est, G1_fix);
        printf("Results: G1_app = %0.20f, G1_est = %0.20f, G1_fix = %0.20f\n\n", G1_app, G1_est, G1_fix);
    }

    fclose(outfile_est_fix);
    if (outfile != NULL) fclose(outfile);
    return 0;
}
