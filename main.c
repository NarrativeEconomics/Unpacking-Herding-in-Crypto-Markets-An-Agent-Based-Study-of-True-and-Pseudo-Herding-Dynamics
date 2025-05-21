#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//Note: These are macros: every time the compiler sees them in the cnetworke, it replaces them with the given value.

#define INITIAL_CAPACITY 100
#define LINE_BUFFER_SIZE 1024

#define NUM_ITERATIONS 10000 // Number of SIMULATED steps
#define MAX_NUMBER_OF_WINDOWS  10 // This relates to computing the market trend function

#define NUMBER_OF_AGENTS 1000
#define NUMBR_OF_SIMULATED_DAYS 655
#define NUMBER_OF_LOOPS (int) ((NUMBER_OF_AGENTS * NUMBR_OF_SIMULATED_DAYS)*0.5)
#define D_BETA 2
#define D_SIGMA 0.233
#define TICK_SIZE 0.01 // This is the size of a market tick

//Note: These are not standard C macros, but are used in the provided code snippet.
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define CLAMP(x, min_val, max_val) (MIN(MAX((x), (min_val)), (max_val)))

//Note: random number generators ".

// return a random number between 0 and 1
double rand_uniform() {
    return (double)rand() / RAND_MAX;
}

// Generate a normally distributed number with given mean and stddev
double rand_normal(double mean, double stddev) {
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    double z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
    return z0 * stddev + mean;
}
//Note: this is a general function compares two values, and returns -1 if val_a > val_b, 1 if val_a < val_b,
//and 0 if they are equal.

// `const void*` means: “a pointer to an object of unknown type (generic) that we promise not to change.”
// this is used in the qsort function to sort the array of doubles in descending order.
// `qsort` is a standard library function that sorts an array using the quicksort algorithm.
int compare_desc(const void* a, const void* b) {
    double val_a = *(double*)a;
    double val_b = *(double*)b;

    if (val_a > val_b) return -1;  // larger value comes first
    if (val_a < val_b) return 1;   // smaller value comes later
    return 0;                      // equal values stay unchanged
}

/*====================================================================================================================*/
/*====================================================================================================================*/
/*                                            General reusable Functions                                              */
/*=================================================================================================================== */
/*====================================================================================================================*/

double* load_csv_data_by_column(const char* filename, int column_index) // functions only has one return so we keep the counetr as a pointer
{

    FILE* csvfileptr = fopen(filename, "r");
        if (csvfileptr == NULL) {
            printf("Error opening file\n%s\n", filename);
            exit(EXIT_FAILURE);
        }
    char line_buffer[LINE_BUFFER_SIZE]; // Creates a character array (a string buffer) where one line from the file will be stored temporarily. Each fgets() call will fill this buffer with a line
    int count = 0; // Initialize counter to 0
    int line_num = 0; // Initialize the line number counter to 0
    int capacity = INITIAL_CAPACITY; // Initialize the capacity of the array

    double* data_array = malloc(capacity * sizeof(double)); // Dynamically allocates memory for an array of double values. (stores the pointer (address) of that memory block.)


    while (fgets(line_buffer, sizeof(line_buffer), csvfileptr)!= NULL) //Read up to 1024 characters into line (which is 1024 bytes big).”
    {
        int current_colum = 0; // Initialize the current column counter to 0

        if (line_num++ == 0) continue; // skip the first line (header)

        char* token = strtok(line_buffer, ","); // This is the first colomn strtok() splits a string into a series of tokens. It takes two arguments: the string to split and the delimiter (in this case, a comma).

        while  (token && current_colum < column_index) {
            token = strtok(NULL, ","); // This moves the pointer to the next token. 	Continue splitting the same string
            current_colum++;
        }
        if (token) {
            if (count >= capacity) {
                double* temp_array = realloc(data_array, 2 * capacity * sizeof(double)); // realloc() changes the size of the memory block pointed to by ptr to size bytes.
                if (temp_array == NULL) {
                    printf("Error reallocating memory\n");
                    exit(EXIT_FAILURE);
                }
                data_array = temp_array;
                capacity *= 2; // Double the capacity
            }

            data_array[count++] = atof(token); // atof converts a string to a double.
        }
    }
    fclose(csvfileptr); // Close the file

    return data_array; // Return the pointer to the array
}
/*====================================================================================================================*/
/*====================================================================================================================*/
/*                                            Structs and functions                                                    */
/*=================================================================================================================== *
/*====================================================================================================================*/

/*============================================Trader and Trader Group ================================================*/

//Note: Define Trader struct
typedef struct {
    double opinion;
    double normalize_opinion;
    double private_info;
    double beta_0, beta_1, beta_2;
    double beta_3, beta_4, beta_5;
} Trader;


//Note: Define TradersGroup struct
typedef struct {
    char name[50]; // Declares a fixed-length character array that holds the group name.
    Trader* traders[NUMBER_OF_AGENTS]; //  an array of pointers to Trader objects, Each entry is a pointer to an Trader, not a copy — so changes to Traders reflect inside the group.
    int count; // Keeps track of how many traders are currently in the group.
} TradersGroup;

//Note: functions used by the Trader and TradersGroup structs



void add_trader(TradersGroup* group_pointer, Trader* trader_pointer) {
    if (group_pointer->count < NUMBER_OF_AGENTS) {
        group_pointer->traders[group_pointer->count] = trader_pointer;
        group_pointer->count++;
    }
}

void remove_trader(TradersGroup* group, Trader* trader) {
    int found_index = -1;

    // Step 1: Find the trader in the list
    for (int i = 0; i < group->count; ++i) {
        if (group->traders[i] == trader) {
            found_index = i;
            break;
        }
    }

    // Step 2: If found, shift remaining traders to the left
    if (found_index != -1) {
        for (int i = found_index; i < group->count - 1; ++i) {
            group->traders[i] = group->traders[i + 1];
        }
        group->count--;  // Decrease group size
    }
}

double average_opinion(TradersGroup* group) {
    if (group->count == 0) return 0.0;  // Avoid division by zero

    double sum = 0.0;

    for (int i = 0; i < group->count; ++i) {
        sum += group->traders[i]->opinion;
    }

    return sum / group->count;
};


double owa_opinion(TradersGroup* group) {

    if (group->count == 0) return 0.0;  // Avoid division by zero

    int n = group->count;
    double opinion_array[n];
    double weights[n];
    double weighted_average = 0.0;

    // Step 1: Copy opinions to array
    for (int i = 0; i < n; ++i) {
        if (group->traders[i] == NULL) {
            printf("owa_opinion(): NULL trader at index %d in group '%s'\n", i, group->name);
            continue;
        }
        opinion_array[i] = group->traders[i]->opinion;
    }

    // Step 2: Sort in descending order
    qsort(opinion_array, n, sizeof(double), compare_desc); // Sort in descending order

    // Step 3: Compute triangular weight denominator to normalize weights giving more weight to stronger opinions
    double denominator = (n * (n + 1)) / 2.0;

    // Step 4: Compute weights based on position favor stronger opinions
    for (int i = 0; i < n; ++i) {
        weights[i] = (n - i) / denominator;
    }

    // Step 5: Compute weighted average
    for (int i = 0; i < n; ++i) {
        weighted_average += opinion_array[i] * weights[i];
    }

    return weighted_average;
};

/*======================================Network and opinions update ==================================================*/

//Note: Define Network struct
typedef struct {
    int N; // total number of traders
    int T; // total number of days
    Trader* traders[NUMBER_OF_AGENTS]; // array of pointers to Trader objects
    TradersGroup NF; // group of nonfundamental traders
    TradersGroup F; // qroup of fundamental traders
    int ag_ids[NUMBER_OF_AGENTS]; // array of ids of participating agent
    int current_day;
    int  ag_ids_len; // length of ag_ids array
    int id1; // selectd agent id
    int id2; // selected agent id
} Opinion_Diffusion;

//Note: functions used by Opinion_Diffusion struct

// select a pair of nonfundamental and fundamental traders randomly
void select_pair(Opinion_Diffusion* Network) {

    int idx_f = rand() % Network->F.count;
    int idx_nf = rand() % Network->NF.count;
    Network->id1=  idx_f;
    Network->id2 = idx_nf;

}
//Note: opinion update function/s -------------------------------------------------------------------------------------

// update opinion of a single trader in the network using Euler's method'
void euler_update(Trader* trader, double beta_0_or_3, double beta_1_or_4, double beta_2_or_5, double I_group, double shared_tanh)
{
    double dt = 1.0; // time step
    double tau = 1.0; // time scale parameter

    // Update opinion of the trader using Euler's method.'
    // Bizyaeva, A., Franci, A., & Leonard, N. E. (2022). "Nonlinear opinion dynamics with tunable sensitivity."

    double dx_dt = (-trader->opinion + beta_0_or_3 * shared_tanh +
        beta_1_or_4 * trader->private_info +
        beta_2_or_5 * I_group) / tau;
    trader->opinion += dx_dt * dt;
}

// update opinion of a single trader in the network
void update_single_trader(Opinion_Diffusion* net, Trader* trader, double I_nf, double I_f, double shared_tanh, double trend) {
    trader->private_info = trend + rand_uniform();

    // find the trader in which group

    for (int i = 0; i < net->NF.count; ++i) {
        if (net->NF.traders[i] == trader) {
            euler_update(trader, trader->beta_0, trader->beta_1, trader->beta_2, I_nf, shared_tanh);
            return; // exit the function as soon as we find the trader
        }
    }
    for (int j = 0; j < net->F.count; ++j) {
        if (net->F.traders[j] == trader) {
            euler_update(trader, trader->beta_3, trader->beta_4, trader->beta_5, I_f, shared_tanh);
            return; // exit the function as soon as we find the trader
        }
    }
}

// the main function to be called to update opinions in the network
void update_opinions(Opinion_Diffusion* network, int id1, int id2, double I_nf, double I_f, double trend) {

    double X_nf = owa_opinion(&network->NF);
    double X_f = owa_opinion(&network->F);
    double shared_tanh = tanh(X_nf + X_f);

    //Hint: maintenance
    if (id1 < 0 || id1 >= NUMBER_OF_AGENTS || network->traders[id1] == NULL) {
        printf("Invalid id1 = %d or null trader pointer!\n", id1);
        return;
    }

    if (id2 < 0 || id2 >= NUMBER_OF_AGENTS || network->traders[id2] == NULL) {
        printf("Invalid id2 = %d or null trader pointer!\n", id2);
        return;
    }
    //Hint: --------
    update_single_trader(network, network->traders[id1], I_nf, I_f, shared_tanh, trend);
    update_single_trader(network, network->traders[id2], I_nf, I_f, shared_tanh, trend);
}

//Note: End opinion update function/s ----------------------------------------------------------------------------------


// Note: 2 version of the pick_agents function. One version uses pow() and the other uses noisy_approximated_volume

//Hint: these functions return the selected agents ids (of the participating agaent)
void pick_agents_1(Opinion_Diffusion* network) {
    int n = (int)network->N;
    double probs[NUMBER_OF_AGENTS];

    for (int i = 0; i < n; ++i) {
        probs[i] = rand_uniform(); // return a float between 0 and 1
    }
    int counter = 0;

    //Hint: this probability has a power law with exponent 2.44
    double p = pow(network->T-network->current_day +1, -2.44) + 0.01; // from prediction markets paper

    for (int j = 0; j < n; ++j) {
            if (probs[j] < p) {
                network->ag_ids[counter++] = j; //select the agent whose id is j and add it to ag_ids array
            }
        }
    network->ag_ids_len = counter;
    }

void pick_agents_2(Opinion_Diffusion* network, double* noisy_approximated_volume) {
    int n = (int)network->N; // make sure the size is intger
    double probs[NUMBER_OF_AGENTS];

    for (int i = 0; i < n; ++i) {
        probs[i] = rand_uniform(); //return a float between 0 and 1
                        }
    int counter = 0;

    for (int i = 0; i < n; ++i) {
        if (probs[i] < (noisy_approximated_volume)[i]) {
            network->ag_ids[counter++] = i;
                    }
             }
    network->ag_ids_len = counter;
    }
// Note: End of 2 version of the pick_agents function-------------------------------------------------------------------


void update_op_seies(Opinion_Diffusion* network, double*  noisy_approximated_volume)
        {
    pick_agents_2(network, noisy_approximated_volume);

}
void lunch(Opinion_Diffusion* network, double I_nf, double I_n, double trend) {
    select_pair(network);
    update_opinions(network, network->id1, network->id2, I_nf, I_n, trend);
}


/*======================================Market struct and functions==================================================*/

typedef struct {
    double pt[NUM_ITERATIONS]; // market price
    double last_pt; // last market price
    double current_trend[MAX_NUMBER_OF_WINDOWS];

    //design parameters
    float d_beta;
    float d_sigma;

    // Other parameters...
    double ED; // excess demand
    double temp_dem[NUM_ITERATIONS];
    double ag_part[NUM_ITERATIONS]; // the participating traders
    int step_tracker;
    int current_day;
}Financial_Market;

//Note: function to update demand and price in the market
void update_demand_and_price(Financial_Market* market, Opinion_Diffusion* network, int* step_tracker) {
    int t = *step_tracker;

    if (t >= NUM_ITERATIONS) return;

    double D = 0.0;
    int count_part = 0;

    for (int i = 0; i < network->ag_ids_len; ++i) {
        int id = network->ag_ids[i];
        if (id < 0 || id >= NUMBER_OF_AGENTS) {
            printf("Warning: Skipping invalid agent ID %d at index %d\n", id, i);
            continue;
        }

        D += network->traders[id]->opinion;
        count_part++;
    }

    if (D != 0.0) {
        double noisy_D = D + rand_normal(0, market->d_sigma);
        double beta_by_noisy_D = market->d_beta * noisy_D;


        // Hint: we normalize to keep it within the range of [0.01, 1.0]
        if (fabs(beta_by_noisy_D) > 1.0) {
            beta_by_noisy_D /= fabs(beta_by_noisy_D);
        }



        double pi_up = fmax(0.0, beta_by_noisy_D);
        double pi_down = fmax(0.0, -beta_by_noisy_D);

        double r = rand_uniform();

        if (r < pi_up) {
            market->last_pt += TICK_SIZE;
        } else if (r < pi_up + pi_down) {
            market->last_pt -= TICK_SIZE;
        }
    }

    // Always update these regardless of D

    market->ED = D;
    market->temp_dem[t] = D;
    market->ag_part[t] = count_part;
    market->last_pt = MAX(0.01, market->last_pt);
    market->pt[t] = market->last_pt;



}


//Note: function to update trends in the market

int* all_possible_window_sizes(int available_length, int* count) {
    *count = available_length -1;
    int* window_sizes = malloc(*count * sizeof(int));
    for (int i = 0; i < *count; i++) {
        window_sizes[i] = i+2;
    }
    return window_sizes;
}
double MSE_compute_slope(double* x, double* y, int n) {
    double sum_x = 0;
    double sum_y = 0;
    double sum_xy = 0;
    double sum_x_squared = 0;

    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x_squared += x[i] * x[i];
    }
    double numerator = n * sum_xy - sum_x * sum_y;
    double denominator = n * sum_x_squared - sum_x * sum_x;
    double slope = numerator / denominator;
    return slope;
}



void update_trends(Financial_Market* market, int* step_tracker) {
    // Initialize trends to 0.0
    for (int i = 0; i < MAX_NUMBER_OF_WINDOWS; i++) {
        market->current_trend[i] = 0.0;
    }

    // Calculate trends
    int available_length = (*step_tracker);
    //find all the available windows sizes from 2 to available_length
    int available_window_sizes_count;

    int* available_window_sizes =  all_possible_window_sizes(available_length, &available_window_sizes_count);

    // SELECT FIVE WINDOWS OF SIZE AT RANDOM

    int selected_window_sizes[5];
    for (int i = 0; i < 5; i++) {
        selected_window_sizes[i] = available_window_sizes[rand() % available_window_sizes_count];
    }


    // FIND THE PRICES (VALUES) OF THE SELECTED WINDOWS
    // NOTE: For each selected window size, extract a window of that many most-recent prices from market->pt[],
    for (int w = 0; w < 5; w++) // loop over all possible window sizes
    {
        int size = selected_window_sizes[w];
        if (size > (*step_tracker)) continue; // Skip if the window size is larger than the available length
        double x[size];  // time steps
        double y[size]; // most recent prices
        for (int i = 0; i < size; i++) {
            x[i] = i;
            y[i] = market->pt[(*step_tracker) - size + i];
        }
        // HINT: using the least squares method to find the slope (m) of the trend line
        double slope = MSE_compute_slope(x, y, size);
        market->current_trend[w] = slope; // store the slope in the market->current_trend array
    }

    // normalize the trends to the range [-1, 1]
    // find the maximum slope
    double max_slope = 0;
    // normalize the trends to the range [-1, 1]
    for (int i = 0; i < 5; i++) {
        double abs_slope = fabs(market->current_trend[i]);
        if (abs_slope > max_slope)
            max_slope = abs_slope;
    }

    if (max_slope > 0.0) {
        for (int i = 0; i < 5; i++) {
            market->current_trend[i] /= max_slope;
            market->current_trend[i] = CLAMP(market->current_trend[i], -1, 1);
        }
    }
    //HINT: free the dynamically allocated memory
    free(available_window_sizes);
    }


double get_random_trend(Financial_Market* market) {
    int index = rand() % 5;
    return market->current_trend[index];  // safe: value may be 0.0, but never NULL
}

//Note: function to lunch market
void lunch_market(Financial_Market* market, Opinion_Diffusion* network) {

    update_demand_and_price(market, network, &market->step_tracker);

}

/* ==================================================================================================================*/
/* ==================================================================================================================*/
/* =================================THIS IS THE MAIN FUNCTION TO RUN SIMULATION===================================== */
/* ==================================================================================================================*/
/* ==================================================================================================================*/

void FC_run_sim(double* params) {
    double beta_0 = params[0];
    double beta_2 = params[1];
    double beta_3 = params[2];
    double beta_5 = params[3];

    double* noisy_approximated_volume_DA = load_csv_data_by_column("BTC_estimated_noisy_volumes.csv", 1);
    double* I_nf_ptr = load_csv_data_by_column("merged_bear_regime_data.csv", 4);
    double* I_f_ptr = load_csv_data_by_column("merged_bear_regime_data.csv", 5);




    double trend = 0.0; // initialize before loop


    FILE* price_ed_file = fopen("prices_and_ed.csv", "w");
    FILE* opinion_file = fopen("opinions.csv", "w");

    if (!price_ed_file || !opinion_file) {
        printf("Error opening output files.\n");
        exit(EXIT_FAILURE);
    }

    //  Write CSV headers
    fprintf(price_ed_file, "time,price,excess_demand\n");


    Opinion_Diffusion network;


    // Opinion network settings
    network.N = NUMBER_OF_AGENTS;
    network.T = NUMBR_OF_SIMULATED_DAYS;
    network.ag_ids_len = 0;
    network.current_day = 0;
    network.id1 = -1;
    network.id2 = -1;





    for (int i = 0; i < NUMBER_OF_AGENTS; i++) {
        network.traders[i] = malloc(sizeof(Trader));
    }
    network.ag_ids_len = NUMBER_OF_AGENTS;

    TradersGroup NF = {.count = 0};
    TradersGroup F  = {.count = 0};




    // Step 5: Set the group's name
    strcpy(NF.name, "Non_Fundamentals");    // strcpy is used to copy strings into char[] in C
    strcpy(F.name, "Fundamentals");


    for (int i = 0; i < NUMBER_OF_AGENTS; i++) {
        if (i % 2 == 0) {
            add_trader(&NF, network.traders[i]);
            network.traders[i]->opinion = 0.0;
            network.traders[i]->beta_0 = beta_0;
            network.traders[i]->beta_1 = (1-beta_0);
            network.traders[i]->beta_2 = beta_2;
            network.traders[i]->beta_3 = 0.0;
            network.traders[i]->beta_4 = 0.0;
            network.traders[i]->beta_5 = 0.0;
            network.traders[i]->private_info = 0.0;
        }
        else {
            add_trader(&F, network.traders[i]);
            network.traders[i]->opinion = 0.0;
            network.traders[i]->beta_0 = 0.0;
            network.traders[i]->beta_1 = 0.0;
            network.traders[i]->beta_2 = 0.0;
            network.traders[i]->beta_3 = beta_3;
            network.traders[i]->beta_4 = (1-beta_3);
            network.traders[i]->beta_5 = beta_5;
            network.traders[i]->private_info = 0.0;
        }
    }

    network.F = F;
    network.NF = NF;

    Financial_Market market;
    // Market settings
    market.step_tracker = 0;
    market.ED = 0.0;
    market.last_pt = TICK_SIZE; // or any starting price
    market.d_beta = D_BETA;
    market.d_sigma = D_SIGMA;



    update_op_seies(&network,noisy_approximated_volume_DA);


    //Note: loop through the simulation
    for (int loop_index=0; loop_index < NUMBER_OF_LOOPS; loop_index++) {

        if (loop_index % (NUMBER_OF_LOOPS/30) == 0)
            {
            trend = get_random_trend(&market);
            }

        market.current_day = loop_index/(int)(network.N *0.5);



        if (market.current_day < NUMBR_OF_SIMULATED_DAYS) {

            lunch(&network, I_nf_ptr[market.current_day], I_f_ptr[market.current_day], trend);
        } else {
            printf("Skipping lunch(): current_day = %d out of bounds\n", market.current_day);
            break;
        }

        if (loop_index % (int)(network.N *0.5)== 0) {
            update_op_seies(&network, noisy_approximated_volume_DA);
            for (int i = 0; i < NUMBER_OF_AGENTS; i++) {
                fprintf(opinion_file, "%f", network.traders[i]->opinion);

                if (i < NUMBER_OF_AGENTS - 1) {
                    fprintf(opinion_file, ",");
                }
            }
            fprintf(opinion_file, "\n");


            lunch_market(&market, &network);
            int t = market.step_tracker;
            fprintf(price_ed_file, "%d,%.6f,%.6f\n", t, market.pt[t], market.ED);

        market.step_tracker++; // Advance time
        }
//Note: loop ends



    }
    for (int i = 0; i < NUMBER_OF_AGENTS; i++) {
        free(network.traders[i]);
    }
    fclose(price_ed_file);
    fclose(opinion_file);

    free(noisy_approximated_volume_DA); // free the dynamically allocated memory after use
    free(I_nf_ptr);
    free(I_f_ptr);
}











/*====================================================================================================================*/


int main() {
    double params[4] = {0.3, 0.2, 0.1, 0.4}; // example parameters
    FC_run_sim(params);
    printf("Simulation complete.\n"); // helpful sanity check
    return 0;
}
