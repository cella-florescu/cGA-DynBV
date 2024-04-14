//
//  main.cpp
//  cGA
//
//  Created by Cella on 26/11/2023.
//

#include <iostream>
#include <algorithm>
#include <random>
#include <assert.h>

using namespace std;

const int MAXN = 1e5;
const int MAXITER = 2e6;
const double EPS = 1e-8;

// INPUT
int N, K;
bool DYNBV;

double UPPERBOUND, LOWERBOUND, CHANGE;

double distribution[MAXN];
int permutation[MAXN];
bool offspring1[MAXN], offspring2[MAXN], lb[MAXN], ub[MAXN], lbaub[MAXN];

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

inline void sample_individual(bool individual[]) {
    for (int i = 0; i < N; ++i) {
        individual[i] = (dis(gen) < distribution[i]);
    }
}

random_device rdd;
mt19937 genn(rdd());

inline int compare_offspring(bool offspring1[], bool offspring2[]) { // -1 if o1 is better, 0 if tie, 1 if o2 is better
    if (DYNBV) {
        shuffle(permutation, permutation + N, genn);
        int index = 0;
        
        while (index < N && offspring1[permutation[index]] == offspring2[permutation[index]]) {
            ++index;
        }

        if (index == N) {
            return 0;
        } else if (offspring1[permutation[index]] > offspring2[permutation[index]]) {
            return -1;
        }
        return 1;
    } else {
        int s1 = 0, s2 = 0;
        for (int i = 0; i < N; ++i) {
            s1 += offspring1[i];
            s2 += offspring2[i];
        }

        if (s1 > s2) {
            return -1;
        } else if (s2 > s1) {
            return 1;
        }
        return 0;
    }
}

inline void update_distribution() {
    
    bool *winner, *loser;
    int result = compare_offspring(offspring1, offspring2);

    if (result == -1) {
        winner = offspring1;
        loser = offspring2;
    } else if (result == 1) {
        winner = offspring2;
        loser = offspring1;
    } else {
        return;
    }
    
    for (int i = 0; i < N; ++i) {
        if (winner[i] > loser[i]) {
            distribution[i] += CHANGE;
        } else if (winner[i] < loser[i]) {
            distribution[i] -= CHANGE;
        }
        distribution[i] = min(max(LOWERBOUND, distribution[i]), UPPERBOUND);
    }
}

inline bool stopping_criterion_achieved() {
    for (int i = 0; i < N; ++i) {
        if (distribution[i] < UPPERBOUND - EPS) {
            return false;
        }
    }
    return true;
}

vector<int> num_probabilities_high, num_probabilities_low, num_probabilities_medium;
vector<double> track_potential;

inline void collect_stats() {
    int h = 0, l = 0, m = 0;
    double potential = N - 1.0;

    for (int i = 0; i < N; ++i) {
        potential -= distribution[i];
        if (distribution[i] > UPPERBOUND - EPS) {
            ++h;
            ub[i] = true;
        } else if (distribution[i] < LOWERBOUND + EPS) {
            ++l;
            lb[i] = true;
            lbaub[i] |= ub[i];
        } else if (distribution[i] > 1.0 / 6.0 && distribution[i] < 5.0 / 6.0) {
            ++m;
        }
    }
    num_probabilities_high.push_back(h);
    num_probabilities_low.push_back(l);
    num_probabilities_medium.push_back(m);
    track_potential.push_back(potential);
}

inline bool is_maximum(bool offspring[]) {
    for (int i = 0; i < N; ++i) {
        if (!offspring[i]) {
            return false;
        }
    }
    return true;
}

int main() {
    int collect_data = 0;
    cin >> N >> K >> DYNBV >> collect_data;
     //N = 1e3; K = 10; DYNBV = false;

    assert(N <= MAXN);

    //UPPERBOUND = 1.0 - 1.0 / (N * log(N));
    //LOWERBOUND = 1.0 / (N * log(N));

    UPPERBOUND = 1.0 - 1.0 / N;
    LOWERBOUND = 1.0 / N;
    CHANGE = 1.0 / K;

    for (int i = 0; i < N; ++i) {
        distribution[i] = 0.5;
        permutation[i] = i;
    }
    
    int iteration = 0;

    bool done = false;

    while (!done && !stopping_criterion_achieved() && iteration < MAXITER) {
        sample_individual(offspring1);
        sample_individual(offspring2);

        done = is_maximum(offspring1) || is_maximum(offspring2);

        update_distribution();
        if (collect_data) {
            collect_stats();
        }
        ++iteration;
    }

    if (collect_data == 1) {

        for (int it : num_probabilities_low) {
            cout << it << " ";
        }
        cout << '\n';

        for (int it : num_probabilities_medium) {
            cout << it << " ";
        }
        cout << '\n';

        for (int it : num_probabilities_high) {
            cout << it << " ";
        }
        cout << '\n';

        for (double it : track_potential) {
            cout << it << " ";
        }
        cout << '\n';
    } else if (collect_data == 2) {
        int ans = 0;
        for (int i = 0; i < N; ++i) {
            ans += lb[i];
        }
        cout << ans;
    } else if (collect_data == 3) {
        int ans = 0;
        for (int i = 0; i < N; ++i) {
            ans += lbaub[i];
        }
        cout << ans;
    } else {
        cout << iteration << '\n';
    }    

    return 0;
}
