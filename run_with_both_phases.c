//
// Created by oc21380 on 5/30/2025.
//
//
// Created by oc21380 on 5/29/2025.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//Note: These are macros: every time the compiler sees them, it replaces them with the given value.

#define LINE_BUFFER_SIZE 1024 // Maximum size of the line buffer (reading files)
#define MAX_NUMBER_OF_WINDOWS  5 // This relates to computing the market trend function (trend function)
#define NUMBER_OF_AGENTS 200 // Number of agents in the market


#define NUMBR_OF_SIMULATED_DAYS 2549   //Bull phase:1,604 days  Bear phase: 945 days

#define NUMBER_OF_LOOPS ((NUMBER_OF_AGENTS * NUMBR_OF_SIMULATED_DAYS) / 2) // half of the total number of agent-days half a chance to  participate every day


#define D_BETA 6 //Design paremeter for price update based on Lux 2000
#define D_SIGMA 0.233 //Design paremeter for price update based on Lux 2000
#define D_LAMBDA 20 //Design paremeter for price update based on farmer and Joshi 2002
#define D_SIGMA_ZETA 0.04 //Design paremeter for price update based on farmer and Joshi 2002





#define TICK_SIZE 0.01 //Size of each price tick

#define PRIVATE_INFO_NOISE_STDDEV 0.05
#define SLOPE_EPSILON 1e-8 // this is used to avoid division by zero in the OWA_slope function


// Note: These macros (CLAMP, MIN, MAX) are not part of the C standard library.
// Note:They are user-defined utility macros commonly used for range-limiting values.

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define CLAMP(x, min_val, max_val) (MIN(MAX((x), (min_val)), (max_val)))

#define DATE_STRING_LENGTH 11  // "YYYY-MM-DD" + null terminator

static char DateStr[NUMBR_OF_SIMULATED_DAYS][DATE_STRING_LENGTH];

const char* OUTPUT_DIR = "C:\\Users\\oc21380\\PycharmProjects\\J_1\\EMSS25\\C_codes\\Claibrate_Black_it_C_code\\Test_Winner\\dump\\";

const char*  INPUT_PATH = "C:\\Users\\oc21380\\CLionProjects\\EMSS25_complete_code\\cmake-build-debug\\BTC_estimated_noisy_combined_volumes_with_betas.csv";
//Note: random number generators.


// Note:Returns a pseudo-random double in the range [0.0, 1.0].
// Note: rand() has limited precision and poor randomness quality.
// TODO: Implement a better random number generator.
double rand_uniform() {
    return (double) rand() / RAND_MAX;
}

// Returns a normally distributed random number with given mean and standard deviation,
// using the Box-Muller transform.
// Clamps u1 away from 0 to avoid log(0) and NaN.
double rand_normal(double mean, double stddev) {
    double u1 = rand_uniform();
    double u2 = rand_uniform();

    // Prevent log(0) and sqrt of negative infinity
    if (u1 < 1e-8) u1 = 1e-8;

    double z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
    return z0 * stddev + mean;
}

//Note: this is a general function compares two values, and returns -1 if val_a > val_b, 1 if val_a < val_b,
//and 0 if they are equal.

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

void load_dates_from_csv(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[LINE_BUFFER_SIZE];
    int count = 0;

    fgets(line, sizeof(line), file); // skip header

    while (fgets(line, sizeof(line), file) && count < NUMBR_OF_SIMULATED_DAYS) {
        char* token = strtok(line, ",");
        if (token) {
            strncpy(DateStr[count], token, DATE_STRING_LENGTH - 1);
            DateStr[count][DATE_STRING_LENGTH - 1] = '\0'; // ensure null-termination
            count++;
        }
    }

    fclose(file);
}




//Note: This function reads a CSV file and stores the values in a dynamically allocated array.
typedef struct {
    double* data;
    int length;
} CsvColumn;

CsvColumn load_csv_data_by_column(const char* filename, int column_index) {


    static double noisy_data[NUMBR_OF_SIMULATED_DAYS];     // index 9  → Noisy Volume
    static double inf_data[NUMBR_OF_SIMULATED_DAYS];       // index 5  → fgi_Index (sentiment info)
    static double iff_data[NUMBR_OF_SIMULATED_DAYS];       // index 6  → MKT_RF_scaled
    static double icf_data[NUMBR_OF_SIMULATED_DAYS];       // index 7  → Fundamental_Index
    static double beta_0_3[NUMBR_OF_SIMULATED_DAYS];       // index 12 → beta_0_3
    static double beta_2[NUMBR_OF_SIMULATED_DAYS];         // index 13 → beta_2
    static double beta_5[NUMBR_OF_SIMULATED_DAYS];         // index 14 → beta_5


    double* data_array = NULL;

    // Choose buffer based on column index
    if (column_index == 5) {
        data_array = inf_data;        // fgi_Index (Sentiment)
    } else if (column_index == 6) {
        data_array = iff_data;        // MKT_RF_scaled (Market factor)
    } else if (column_index == 7) {
        data_array = icf_data;        // Fundamental_Index commposite
    } else if (column_index == 9) {
        data_array = noisy_data;      // Noisy Volume
    } else if (column_index == 12) {
        data_array = beta_0_3;        // Beta for regime-based reaction
    } else if (column_index == 13) {
        data_array = beta_2;
    } else if (column_index == 14) {
        data_array = beta_5;
    } else {
        printf("Unsupported column index: %d\n", column_index);
        exit(EXIT_FAILURE);
    }


    FILE* csvfileptr = fopen(filename, "r");
    if (csvfileptr == NULL) {
        printf("Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line_buffer[LINE_BUFFER_SIZE];
    int count = 0;
    int line_num = 0;
    double last_valid = 0.0;

    while (fgets(line_buffer, sizeof(line_buffer), csvfileptr) != NULL) {
        int current_column = 0;
        if (line_num++ == 0) continue;

        char* token = strtok(line_buffer, ",");
        while (token && current_column < column_index) {
            token = strtok(NULL, ",");
            current_column++;
        }

        if (token == NULL) continue;

        while (*token == ' ' || *token == '\t') token++;

        if (*token == '\0' || *token == '\n') {
            data_array[count++] = last_valid;
        } else {
            double val = atof(token); // Convert string to double
            last_valid = val;
            data_array[count++] = val;
        }

        if (count >= NUMBR_OF_SIMULATED_DAYS) break;
    }

    fclose(csvfileptr);

    CsvColumn result = { data_array, count };
    return result;
}





//Note: this is to shuffle indexes
typedef struct {
    int even_indices[NUMBER_OF_AGENTS / 2];
    int odd_indices[NUMBER_OF_AGENTS / 2];
    int even_pos;
    int odd_pos;
} AgentSampler;

void init_sampler(AgentSampler* s, int N) {
    int even_count = 0, odd_count = 0;

    for (int i = 0; i < N; ++i) {
        if (i % 2 == 0) s->even_indices[even_count++] = i;
        else            s->odd_indices[odd_count++] = i;
    }

    // Fisher-Yates shuffle
    for (int i = even_count - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        int tmp = s->even_indices[i];
        s->even_indices[i] = s->even_indices[j];
        s->even_indices[j] = tmp;
    }

    for (int i = odd_count - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        int tmp = s->odd_indices[i];
        s->odd_indices[i] = s->odd_indices[j];
        s->odd_indices[j] = tmp;
    }

    s->even_pos = 0;
    s->odd_pos = 0;
}
int next_even(AgentSampler* s, int total_even) {
    if (s->even_pos >= total_even) {
        init_sampler(s, total_even * 2); // re-shuffle
    }
    return s->even_indices[s->even_pos++];
}

int next_odd(AgentSampler* s, int total_odd) {
    if (s->odd_pos >= total_odd) {
        init_sampler(s, total_odd * 2); // re-shuffle
    }
    return s->odd_indices[s->odd_pos++];
}






/*====================================================================================================================*/
/*====================================================================================================================*/
/*                                            Structs and functions                                                    */
/*=================================================================================================================== *
/*====================================================================================================================*/

/*============================================Trader and Trader Group ================================================*/



//Note: Define Trader struct
typedef struct {
    int index;
    char group[50];
    double opinion;
    double last_opinion;
    double normalize_opinion;
    double last_normalized_opinion;
    double private_info;
    double beta_0, beta_1, beta_2;
    double beta_3, beta_4, beta_5;
    double position;
    double last_position;
} Trader;

typedef struct {
    double position;
} MarketMaker;

//Note: Define TradersGroup struct
typedef struct {
    char name[50];
    Trader* traders[NUMBER_OF_AGENTS]; //  an array of pointers to Trader objects
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

    // Step 1: Find the trader in the array
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

    // Step 1: Copy opinions to an array
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
}

/*======================================Network and opinions update ==================================================*/

//Note: Define Network struct
typedef struct Opinion_Diffusion {
    int N; // total number of traders
    int T; // total number of days
    Trader* traders[NUMBER_OF_AGENTS]; // array of pointers to Trader objects (keep track of all traders in the network)
    TradersGroup NF; // group of nonfundamental traders
    TradersGroup F; // qroup of fundamental traders
    int ag_ids[NUMBER_OF_AGENTS]; // array of ids of participating agents
    int current_day;
    int  ag_ids_len; // length of ag_ids array
    int id1; // selectd agent id
    int id2; // selected agent id
} Opinion_Diffusion;

//Note: functions used by Opinion_Diffusion struct

// select a pair of nonfundamental and fundamental traders randomly
void select_pair_old(Opinion_Diffusion* Network) {
    int max_even = (Network->N % 2 == 0) ? Network->N - 2 : Network->N - 1;
    int max_odd = (Network->N % 2 == 0) ? Network->N - 1 : Network->N - 2;

    // Random even index for F group
    Network->id1 = 2 * (rand() % ((max_even / 2) + 1));

    // Random odd index for NF group
    Network->id2 = 2 * (rand() % ((max_odd / 2) + 1)) + 1;


}
void select_pair(Opinion_Diffusion* Network, AgentSampler* sampler) {

    Network->id1 = next_even(sampler, Network->F.count);
    Network->id2 = next_odd(sampler, Network->NF.count);

}


//Note: opinion update function/s -------------------------------------------------------------------------------------

// update opinion of a single trader in the network using Euler's method'
void euler_update_nf(Trader* trader, double beta_0_or_3, double beta_1_or_4, double beta_2_or_5, double I_group, double shared_tanh)
{
    double dt = 1.0; // time step
    double tau = 1.0; // time scale parameter

    // Update opinion of the trader using Bizyaeva, A., Franci, A., & Leonard, N. E. (2022).
    //"Nonlinear opinion dynamics with tunable sensitivity." + my modifications

    double dx_dt = (-trader->opinion + beta_0_or_3 * shared_tanh +
        beta_1_or_4 * trader->private_info +
        beta_2_or_5 * I_group) / tau;
    trader->last_opinion = trader->opinion; // save the previous opinion before updating it
    trader->opinion += dx_dt * dt;
    trader->last_normalized_opinion = trader->normalize_opinion;
    trader->normalize_opinion =  ( trader->opinion  -(-2)) / (2-(-2)); // normalize opinion

}

void euler_update_f(Trader* trader, double beta_0_or_3, double beta_1_or_4, double beta_2_or_5, double I_group, double I_cf, double shared_tanh)
{
    double dt = 1.0; // time step
    double tau = 1.0; // time scale parameter

    // Update opinion of the trader using Bizyaeva, A., Franci, A., & Leonard, N. E. (2022).
    //"Nonlinear opinion dynamics with tunable sensitivity." + my modifications

    double dx_dt = (-trader->opinion + beta_0_or_3 * shared_tanh +
    beta_1_or_4 * trader->private_info +
        beta_2_or_5 * I_group+
        (1 - beta_2_or_5) * I_cf) / tau;
    trader->last_opinion = trader->opinion; // save the previous opinion before updating it
    trader->opinion += dx_dt * dt;
    trader->last_normalized_opinion = trader->normalize_opinion;
    trader->normalize_opinion =  ( trader->opinion  -(-2)) / (2-(-2)); // normalize opinion

}

// update opinion of a single trader in the network
void update_single_trader(Opinion_Diffusion* net, Trader* trader, double I_nf, double I_f, double I_cf,  double X_nf,double X_f, double trend,
    double beta_0_3, double beta_2, double beta_5) {

    trader->private_info = CLAMP(trend + rand_normal(0, PRIVATE_INFO_NOISE_STDDEV),-1,1);
    //printf("trend: %.2f, Private info: %.2f\n", trend, trader->private_info);
    double shared_tanh = tanh(X_nf+ X_f);

    // find the trader in which group

    for (int i = 0; i < net->NF.count; ++i) {
        if (net->NF.traders[i] == trader) {
            trader->beta_0 = beta_0_3;
            trader->beta_1 = 1-beta_0_3;
            trader->beta_2 = beta_2;
            euler_update_nf(trader, trader->beta_0, trader->beta_1, trader->beta_2, I_nf, shared_tanh);
            return; // exit the function as soon as we find the trader
        }
    }
    for (int j = 0; j < net->F.count; ++j) {
        if (net->F.traders[j] == trader) {
            trader->beta_3 = beta_0_3;
            trader->beta_4 = 1-beta_0_3;
            trader->beta_5 = beta_5;  // beta_5 is now beta_2 in the paper (see equation 3 in the paper)
            euler_update_f(trader, trader->beta_3, trader->beta_4, trader->beta_5, I_f, I_cf, shared_tanh);
            return; // exit the function as soon as we find the trader
        }
    }
}

// the main function to be called to update opinions in the network
void update_opinions(Opinion_Diffusion* network, int id1, int id2, double I_nf, double I_f, double I_cf, double trend, double beta_0_3,
    double beta_2, double beta_5) {

    double X_nf = owa_opinion(&network->NF);
    double X_f = owa_opinion(&network->F);


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

    update_single_trader(network, network->traders[id1], I_nf, I_f, I_cf, X_nf, X_f, trend, beta_0_3, beta_2, beta_5);
    update_single_trader(network, network->traders[id2], I_nf, I_f, I_cf,  X_nf, X_f, trend, beta_0_3, beta_2, beta_5);
}

//Note: End opinion update function/s ----------------------------------------------------------------------------------


// Note: 2 versions of the pick_agents function. One version uses pow_law and the other uses noisy_approximated_volume

//Hint: these functions return the selected agents ids (of the participating agaent)


void pick_agents_2(Opinion_Diffusion* network, const double* noisy_approximated_volume) {

    int n = (int)network->N;
    double probs[NUMBER_OF_AGENTS];

    for (int i = 0; i < n; ++i) {
        probs[i] = rand_uniform(); //return a float between 0 and 1
                        }
    int counter = 0;

    for (int j = 0; j < n; ++j) {
        if (probs[j] < (noisy_approximated_volume)[j]) {

            if (counter < NUMBER_OF_AGENTS) {

                network->ag_ids[counter++] = j;//select the agent whose id is j and add it to ag_ids array
            }
        }
    }
    network->ag_ids_len = counter;
    }
// Note: End of 2 version of the pick_agents function-------------------------------------------------------------------


void update_op_series(Opinion_Diffusion* network, double*  noisy_approximated_volume)
        {
    pick_agents_2(network, noisy_approximated_volume);

}


/*void execute_trade(Opinion_Diffusion* net, int id1, int id2, double trade_volume) {
    Trader* a = net->traders[id1];
    Trader* b = net->traders[id2];

    // Decide roles: the one with higher opinion buys
    if (a->opinion > b->opinion) {
        a->position += trade_volume;
        b->position -= trade_volume;
    } else {
        a->position -= trade_volume;
        b->position += trade_volume;
    }
}*/


void lunch_opinion(Opinion_Diffusion* network, double I_nf, double I_n, double I_cf, double trend, AgentSampler* sampler, double beta_0_3,
    double beta_2, double beta_5, MarketMaker* mm) {
    //select_pair_old(network);

    select_pair(network, sampler);
    update_opinions(network, network->id1, network->id2, I_nf, I_n, I_cf, trend, beta_0_3,beta_2, beta_5);
    //double trade_volume = 1.0;  // or some function of opinion difference
    //execute_trade(network, network->id1, network->id1, trade_volume);

}



/*======================================Market struct and functions==================================================*/

typedef struct Financial_Market{
    double* pt; // market price array
    double last_pt; // last market price
    double current_trend[MAX_NUMBER_OF_WINDOWS];
    double prev_avg_opinion;  // store the previous average opinion

    //design parameters
    float d_beta;
    float d_sigma;

    // Other parameters...
    double ED; // excess demand
    double* temp_dem; //TODO: do I need this?
    double* ag_part; // the participating traders
    int step_tracker;
    int current_day;
}Financial_Market;

//Note: function to update demand and price in the market
void update_demand_and_price_lux(Financial_Market* market, Opinion_Diffusion* network, const int* const step_tracker) {

    int t = *step_tracker;
    if (t >= NUMBER_OF_LOOPS) return;

    double opinion_sum = 0.0;
    int count_part = 0;

    for (int i = 0; i < network->ag_ids_len; ++i) {
        int id = network->ag_ids[i];
        opinion_sum += network->traders[id]->opinion;
        count_part++;
    }

    if (count_part > 0) {
        double avg_opinion = opinion_sum / count_part;
        //double delta_opinion = avg_opinion - market->prev_avg_opinion;
       // market->ED = count_part * delta_opinion; // demand
        market->ED = avg_opinion;
        //printf("avg_opinion=%f, delta_opinion=%f, ED=%f\n", avg_opinion, delta_opinion, market->ED);
       // printf("avg_opinion=%f, ED=%f\n", avg_opinion, market->ED);
        double mu = rand_normal(0.0, market->d_sigma); // Noise term ~ N(0, sigma^2)
        double dt = 1; // Time increment (simulating "infinitesimal")


        // Compute price up/down probabilities based on Lux-Marchesi equation
        double prob_up = MAX(0.0,  (D_BETA * market->ED + mu)* dt) ;
        double prob_down =  -MIN(0.0, (D_BETA * market->ED + mu) * dt);

        // Draw random number to determine direction
        double u = rand_uniform();
        //printf("prob_up=%f, prob_down=%f, u=%f\n", prob_up, prob_down, u);


        if (u < prob_up) {
            //printf("Price goes up\n");
            market->last_pt += TICK_SIZE;
        } else if (u < prob_up + prob_down) {
            //printf("Price goes down\n");
            market->last_pt -= TICK_SIZE;
        } else {
            market->last_pt = market->last_pt;
            //printf("Price stays the same\n");
            // No price change (optional, Lux paper assumes this happens too)
        }

        market->prev_avg_opinion = avg_opinion;  // Save for next step


    }

    market->ag_part[t] = count_part;
    market->last_pt = MAX(0.01, market->last_pt);
    market->pt[t] = market->last_pt;


}




void update_demand_and_price(Financial_Market* market, Opinion_Diffusion* network, const int* const step_tracker) {

    int t = *step_tracker;
    if (t >= NUMBER_OF_LOOPS) return;

    double opinion_sum = 0.0;
    int count_part = 0;

    for (int i = 0; i < network->ag_ids_len; ++i) {
        int id = network->ag_ids[i];
        opinion_sum += network->traders[id]->opinion;
        count_part++;
    }

    if (count_part > 0) {
        double avg_opinion = opinion_sum / count_part;
        double delta_opinion = avg_opinion - market->prev_avg_opinion;

        market->last_pt = market->last_pt + delta_opinion;
        market->ED = delta_opinion;
        market->prev_avg_opinion = avg_opinion;  // Save for next step

    }

    market->ag_part[t] = count_part;
    market->last_pt = MAX(0.01, market->last_pt);
    market->pt[t] = market->last_pt;


}

double compute_recent_mean(Financial_Market* market, int window_size) {
    int t = market->step_tracker;
    if (t < window_size) return market->last_pt;  // fallback if too early

    double sum = 0.0;
    for (int i = t - window_size; i < t; i++) {
        sum += market->pt[i];
    }
    return sum / window_size;
}

void update_demand_and_price_farmer2002(Financial_Market* market, Opinion_Diffusion* network, const int* const step_tracker) {

    int t = *step_tracker;
    if (t >= NUMBER_OF_LOOPS) return;

    double net_demand= 0.0;
    int count_part = 0;
    double D_i;
    double log_p_t;
    double zeta;


    for (int i = 0; i < network->ag_ids_len; ++i) {
        int id = network->ag_ids[i];
        D_i = network->traders[id]->opinion - network->traders[id]->last_opinion;
        net_demand += D_i;
        count_part++;
    }
    market->ag_part[t] = count_part;
    if (count_part > 0) {
        market->ED = net_demand/count_part;
    }
    else {market->ED = net_demand;}


    zeta = rand_normal(0.0,D_SIGMA_ZETA); // Noise term ~ N(0, sigma^2)

    double mu_t = compute_recent_mean(market, 500);  // Same window as in data
    double theta =  0.0158
;
    log_p_t = market->last_pt
            + (market->ED / D_LAMBDA)
            + zeta;
          //  + theta * (mu_t - market->last_pt);  // mean-reverting pull
    // random noise



    market->last_pt = log_p_t;
    market->pt[t] = market->last_pt;
    //printf("t=%d | Net demand: %f | log_p_t: %f | exp(log_p_t): %f\n", t, market->ED, log_p_t, market->last_pt);

}

//Note: function to update trends in the market
void all_possible_window_sizes(int available_length, int* out_array, int* count) {
    if (available_length < 2) {
        *count = 0;
        return;
    }

    *count = available_length - 1;

    for (int i = 0; i < *count; i++) {
        out_array[i] = i + 2;
    }
}


//Note: computing the slope as a closed-form soluiton for the Ordinary Least Squares (OLS) formula
//https://www.geeksforgeeks.org/least-square-method/

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


void update_trends(Financial_Market* market, const int* const step_tracker) {
    // Initialize trends to 0.0
    for (int i = 0; i < MAX_NUMBER_OF_WINDOWS; i++) {
        market->current_trend[i] = 0.0;
    }

    int available_length = *step_tracker;
    int available_window_sizes_count;

    // Preallocate max possible size (safe upper bound)
    int available_window_sizes[NUMBR_OF_SIMULATED_DAYS];  // assuming MAX_WINDOW_SIZE ≥ expected max
    all_possible_window_sizes(available_length, available_window_sizes, &available_window_sizes_count);

    if (available_window_sizes_count < 5) {
        printf("Not enough window sizes available. Skipping trend update.\n");
        return;
    }

    // Pick 5 random window sizes
    int selected_window_sizes[5];
    for (int i = 0; i < 5; i++) {
        selected_window_sizes[i] = available_window_sizes[rand() % available_window_sizes_count];
    }

    // Predeclare largest needed arrays for x and y
    double x[NUMBR_OF_SIMULATED_DAYS];
    double y[NUMBR_OF_SIMULATED_DAYS];

    // Compute slope for each selected window
    for (int w = 0; w < 5; w++) {
        int size = selected_window_sizes[w];
        if (size > *step_tracker || size > NUMBR_OF_SIMULATED_DAYS) continue;

        for (int i = 0; i < size; i++) {
            x[i] = i;
            y[i] = market->pt[*step_tracker - size + i];
        }

        double slope = MSE_compute_slope(x, y, size);
        market->current_trend[w] = slope;
    }

    // Normalize only if meaningful variation exists
    double max_slope = fabs(market->current_trend[0]);
    double min_slope = fabs(market->current_trend[0]);

    for (int i = 1; i < 5; i++) {
        double abs_slope = fabs(market->current_trend[i]);
        if (abs_slope > max_slope) max_slope = abs_slope;
        if (abs_slope < min_slope) min_slope = abs_slope;
    }

    double range = max_slope - min_slope;
    if (range > 1e-4 && max_slope > SLOPE_EPSILON) {
        for (int i = 0; i < 5; i++) {
            market->current_trend[i] /= max_slope;
            market->current_trend[i] = CLAMP(market->current_trend[i], -1, 1);
        }
    }
}


double get_random_trend(Financial_Market* market) {
    int index = rand() % 5; //Hint: we have 5 trends in the market->current_trend array
    return market->current_trend[index];
}

//Note: function to lunch market
void lunch_market(Financial_Market* market, Opinion_Diffusion* network) {

  update_demand_and_price_farmer2002(market, network, &market->step_tracker);
    //update_demand_and_price(market, network, &market->step_tracker);

}


/* ==================================================================================================================*/
/* ==================================================================================================================*/
/* =================================THIS IS THE MAIN FUNCTION TO RUN SIMULATION===================================== */
/* ==================================================================================================================*/
/* ==================================================================================================================*/
//void FC_run_sim(const char* price_filename, const char* opinion_filename)
void FC_run_sim(const char* price_filename)
    {

  //  srand(1);


    static double result[NUMBER_OF_LOOPS];
    static Trader trader_pool[NUMBER_OF_AGENTS]; // array of traders
    static Opinion_Diffusion network_storage;
    static Financial_Market market_storage;
    static double pt_array[NUMBER_OF_LOOPS];
    static double temp_dem_array[NUMBER_OF_LOOPS];
    static double ag_part_array[NUMBER_OF_LOOPS];

    memset(result, 0, sizeof(result));
    //memset(trader_pool, 0, sizeof(trader_pool));
    memset(pt_array, 0, sizeof(pt_array));
    memset(temp_dem_array, 0, sizeof(temp_dem_array));
    memset(ag_part_array, 0, sizeof(ag_part_array));


    Opinion_Diffusion* network = &network_storage;
    Financial_Market* market = &market_storage;



    CsvColumn noisy     = load_csv_data_by_column(INPUT_PATH, 9);   // Noisy Volume
    CsvColumn inf       = load_csv_data_by_column(INPUT_PATH, 5);   // fgi_Index
    CsvColumn iff       = load_csv_data_by_column(INPUT_PATH, 6);   // MKT_RF_scaled
    CsvColumn icf       = load_csv_data_by_column(INPUT_PATH, 7);   // Fundamental_Index
    CsvColumn beta_0_3  = load_csv_data_by_column(INPUT_PATH, 12);  // beta_0_3
    CsvColumn beta_2    = load_csv_data_by_column(INPUT_PATH, 13);  // beta_2
    CsvColumn beta_5    = load_csv_data_by_column(INPUT_PATH, 14);  // beta_5




    load_dates_from_csv(INPUT_PATH);

    //------------------------------------OUTPUT FILES PERPATION -----------------------------------------------------//

    FILE* price_ed_file = fopen(price_filename, "w");
    if (!price_ed_file) {
        fprintf(stderr, "Error: could not open %s for writing\n", price_filename);
        exit(EXIT_FAILURE);
    }

   //FILE* opinion_file = fopen(opinion_filename, "w");
   //if (!opinion_file) {
   //     fprintf(stderr, "Error: could not open %s for writing\n", opinion_filename);
   //     exit(EXIT_FAILURE);
   // }


    //  Write CSV headers
    fprintf(price_ed_file, "Date, step_traker,day,price,excess_demand\n");
    //----------------------------------------------------------------------------------------------------------------//



    network->N = NUMBER_OF_AGENTS;
    network->T = NUMBR_OF_SIMULATED_DAYS;
    network->ag_ids_len = 0;
    network->current_day = 0;


    TradersGroup NF = {.count = 0};
    TradersGroup F  = {.count = 0};

    strcpy(NF.name, "Non_Fundamentals");
    strcpy(F.name, "Fundamentals");


    static MarketMaker mm = { .position = 0.0 };  // initialized with zero position

    // initialize traders
    for (int i = 0; i < NUMBER_OF_AGENTS; i++) {
        network->traders[i] = &trader_pool[i];
        network->traders[i]->index = i;
        network->traders[i]->opinion = 0.0;
        network->traders[i]->last_opinion = network->traders[i]->opinion;
        network->traders[i]->private_info = 0.0;
        network->traders[i]->position = 0.0;
        network->traders[i]->last_position= 0.0;

        if (i % 2 == 0) {
            add_trader(&NF, network->traders[i]);
            strcpy( network->traders[i]->group, NF.name);
            network->traders[i]->beta_0 = 0.0;
            network->traders[i]->beta_1 = 0.0;
            network->traders[i]->beta_2 = 0.0;
            network->traders[i]->beta_3 = 0.0;
            network->traders[i]->beta_4 = 0.0;
            network->traders[i]->beta_5 = 0.0;
        } else {
            add_trader(&F, network->traders[i]);
            strcpy(network->traders[i]->group, F.name);
            network->traders[i]->beta_0 = 0.0;
            network->traders[i]->beta_1 = 0.0;
            network->traders[i]->beta_2 = 0.0;
            network->traders[i]->beta_3 = 0.0;
            network->traders[i]->beta_4 = 0.0;
            network->traders[i]->beta_5 = 0.0;
        }
    }

    network->F = F;
    network->NF = NF;

    market->pt = pt_array;
    market->temp_dem = temp_dem_array;
    market->ag_part = ag_part_array;
    market->step_tracker = 0;
    market->current_day = 0;
    market->prev_avg_opinion = 0.0;
    market->ED = 0.0;
    market->last_pt = 9;
    market->d_beta = D_BETA;
    market->d_sigma = D_SIGMA;

    AgentSampler sampler;
    init_sampler(&sampler, NUMBER_OF_AGENTS);

    double trend = 0.0;

    for (int loop_index = 0; loop_index < NUMBER_OF_LOOPS; loop_index++) {
        if (market->step_tracker >= 6 && market->current_day % 30 == 0) {
            update_trends(market, &market->step_tracker);
            trend = get_random_trend(market);


        }

        market->current_day = loop_index / (NUMBER_OF_AGENTS / 2);


        lunch_opinion(network, inf.data[market->current_day], iff.data[market->current_day],icf.data[market->current_day], trend, &sampler,
            beta_0_3.data[market->current_day], beta_2.data[market->current_day], beta_5.data[market->current_day], &mm);





        if (loop_index % (NUMBER_OF_AGENTS / 2) == 0)
            {
            update_op_series(network, noisy.data);
          //  fprintf(opinion_file, "%s",DateStr[market->current_day]);
          //  fprintf(opinion_file, ",");
          //  for (int i = 0; i < NUMBER_OF_AGENTS; i++) {
          //      fprintf(opinion_file, "%f", network->traders[i]->opinion);
           //     if (i < NUMBER_OF_AGENTS - 1) {
           //         fprintf(opinion_file, ",");
            //    }
          //  }
         //   fprintf(opinion_file, "\n");
            }



        lunch_market(market, network);

        //result[loop_index] = market->pt[market->step_tracker];
        fprintf(price_ed_file, "%s,%d,%d,%f,%f\n",
         DateStr[market->current_day],
         loop_index,
         market->current_day,
         market->pt[loop_index],
         market->ED);


        market->step_tracker++;
    }


    fclose(price_ed_file);
   // fclose(opinion_file);
}








// Forward declaration of FC_run_sim—with its new signature.



//void FC_run_sim(const char* price_filename, const char* opinion_filename);
void FC_run_sim(const char* price_filename);

int main(void) {
    // 1) Define two sub‐directories (must already exist on disk):
    const char* prices_dir   =
        "C:\\Users\\oc21380\\PycharmProjects\\J_1\\EMSS25\\C_codes"
        "\\Claibrate_Black_it_C_code_Farmer_Joshi2002\\Test_Winner\\dump\\prices\\";
   // const char* opinions_dir =
    //    "C:\\Users\\oc21380\\PycharmProjects\\J_1\\EMSS25\\C_codes"
   //     "\\Claibrate_Black_it_C_code_Farmer_Joshi2002\\Test_Winner\\dump\\opinions\\";

    // 2) Buffers to hold the two filenames
    char price_fn[512];
    char opinion_fn[512];

    // 3) Loop over 1000 runs (0..999)
    for (int run =0; run < 1000; ++run) {
        srand(run);  // seed randomness uniquely per run



        // 3a) Build the price filename inside dump\prices\
        //     e.g. "…\dump\prices\run_0000_prices_and_ed.csv"
        sprintf(price_fn,
                "%srun_%04d_prices_and_ed_Farmer_Joshi2002_full.csv",
                prices_dir,
                run);

        // 3b) Build the opinion filename inside dump\opinions\
        //     e.g. "…\dump\opinions\run_0000_opinions.csv"
      // sprintf(opinion_fn,
      //          "%srun_%04d_opinions_Farmer_Joshi2002_full.csv",
       //         opinions_dir,
      //          run);

        // 3c) Call the simulation, passing those two paths
       // FC_run_sim(price_fn, opinion_fn);
         FC_run_sim(price_fn);

        // 3d) Optional progress message
        printf("Completed run %d / 1000\n", run + 1);
    }

    return 0;
}