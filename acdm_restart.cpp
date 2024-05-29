#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
// #pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <random>
#include <chrono>
#include <queue>
#include <unistd.h>

#if PARALLEL
#include <omp.h>
#endif

using namespace std;
#include "utilities.h"

vector<double> b_star;

void acdm(vector<vector<pair<int, int>>> &Adj,
          vector<int> &edgeMap,
          long long iterations,
          long long n,
          long long m,
          int recalc,
          bool provided_best_loads = false)
{

    vector<double> z(2 * m), b(n), sw(n), sy(n), y(2 * m), w(2 * m);
    double theta;
    vector<int> perm;
    for (int i = 0; i < m; ++i)
        perm.push_back(i);
    z = get_initial(Adj, n, m);

    for (int i = 0; i < n; ++i)
    {
        b[i] = 0;
        sw[i] = 0;
        sy[i] = 0;
        for (auto &j : Adj[i])
        {
            b[i] += z[j.second];
        }
    }

    for (int i = 0; i < 2 * m; ++i)
    {
        w[i] = 0;
        y[i] = z[i];
    }
    double density = current_density(Adj, b, z, n, m);
    double density_sorting = 0.0;
    double fxk = 0;
    for (int i = 0; i < n; ++i)
        fxk += b[i] * b[i];
    double last_fxk = 0;
    for (int t = 1; t <= iterations; ++t)
    {
        auto start1 = std::chrono::high_resolution_clock::now();
        std::random_shuffle(perm.begin(), perm.end());
        for (int s = 0; s < m; ++s)
        {
            if (s == 0 && last_fxk < fxk)
            {
                for (int i = 0; i < 2 * m; ++i)
                {
                    z[i] = theta * theta * w[i] + y[i];
                    y[i] = z[i];
                    w[i] = 0;
                }

                for (int i = 0; i < n; ++i)
                {
                    b[i] = 0;
                    sw[i] = 0;
                    sy[i] = 0;
                    for (auto &j : Adj[i])
                    {
                        b[i] += z[j.second];
                        sw[i] += w[j.second];
                        sy[i] += y[j.second];
                    }
                }

                theta = 1.0 / m;
            }
            else
            {
                theta = (sqrt(pow(theta, 4) + 4 * pow(theta, 2)) - pow(theta, 2)) / 2;
            }
            // pick a random edge
            int e = perm[s];
            int eu = 2 * e;
            int u = edgeMap[eu];
            int ev = 2 * e + 1;
            int v = edgeMap[ev];

            double y_eu = y[eu] - 1 / (theta * 2 * m) * (theta * theta * sw[u] + sy[u]);
            double y_ev = y[ev] - 1 / (theta * 2 * m) * (theta * theta * sw[v] + sy[v]);

            if (abs(y_eu - y_ev) <= 1)
            {
                y_eu = (y_eu - y_ev + 1) / 2;
            }
            else if (y_eu - y_ev < -1)
            {
                y_eu = 0;
            }
            else
            {
                y_eu = 1;
            }
            y_ev = 1 - y_eu;

            double w_eu = w[eu] - (1 - m * theta) / (theta * theta) * (y_eu - y[eu]);
            double w_ev = w[ev] - (1 - m * theta) / (theta * theta) * (y_ev - y[ev]);

            last_fxk = fxk;

            fxk = fxk - b[u] * b[u] - b[v] * b[v];

            sw[u] = sw[u] - w[eu] + w_eu;
            sw[v] = sw[v] - w[ev] + w_ev;
            sy[u] = sy[u] - y[eu] + y_eu;
            sy[v] = sy[v] - y[ev] + y_ev;

            b[u] = theta * theta * sw[u] + sy[u];
            b[v] = theta * theta * sw[v] + sy[v];

            fxk = fxk + b[u] * b[u] + b[v] * b[v];

            y[eu] = y_eu;
            y[ev] = y_ev;
            w[eu] = w_eu;
            w[ev] = w_ev;
        }

        for (int i = 0; i < 2 * m; ++i)
        {
            z[i] = theta * theta * w[i] + y[i];
        }
        auto finish1 = std::chrono::high_resolution_clock::now();

        /*
        Print diagnostics for experimental section, can safely remove this for non experiments.
        */
        cout << "Iteration=" << t << endl;
        cout << "Time=" << std::chrono::duration_cast<chrono::milliseconds>(finish1 - start1).count() << " miliseconds" << endl;
        if (provided_best_loads)
        {
            double E = 0;
            for (int i = 0; i < n; ++i)
            {
                E += (b[i] - b_star[i]) * (b[i] - b_star[i]);
            }
            E = sqrt(E);
            cout << "Error=" << E << endl;
        }
        for (int i = 0; i < n; ++i)
            b[i] *= t;
        if (t % recalc == 0)
            density = max(density, current_density(Adj, b, z, n, m, t));

        cout << "Best density=" << density << endl;
        density_sorting = max(density_sorting, current_density_sorting(Adj, b, n, m));
        cout << "Sorting density=" << density_sorting << endl;

        double sum_of_squares = 0.0;
        for (int i = 0; i < n; ++i)
            sum_of_squares += b[i] * b[i] / (t * t);
        sum_of_squares = sqrt(sum_of_squares);
        cout << "Sum of squares=" << sum_of_squares << endl;

        cout << endl;

        if (t >= iterations)
        {
            cout << "b vector=[";
            for (int i = 0; i < n; ++i)
                cout << b[i] / t << ",";
            cout << "]" << endl;
            break;
        }
    }
}

int main(int argc, char **argv)
{
    ios_base::sync_with_stdio(0);
    cin.tie(0);

    cout << "Starting Graph reading" << endl;

    auto start = std::chrono::high_resolution_clock::now();
    long long n, m;
    cin >> n >> m;

    vector<vector<pair<int, int>>> Adj(n); // neighbor and edge index
    vector<int> edgeMap(2 * m);

    int e_idx = 0;
    for (int e = 0; e < m; ++e)
    {
        int i, j;
        cin >> i >> j;
        Adj[i].push_back({j, e_idx});
        edgeMap[e_idx] = i;
        Adj[j].push_back({i, e_idx + 1});
        edgeMap[e_idx + 1] = j;
        e_idx += 2;
    }
    string line;
    getline(cin, line); // read \n

    bool provided_best_loads = false;
    b_star = vector<double>(n, 0); // optimal load vector
    int i = 0;
    while (getline(cin, line))
    {
        provided_best_loads = true;
        b_star[i++] = stod(line);
    }

    auto end = std::chrono::high_resolution_clock::now();
    cout << "Time for reading input is " << std::chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    cout << endl;

    int iters = atoi(argv[1]);
    int recalc = 1;
    // bool print_diagnostics = true;

    // if (argc > 2)
    //     recalc = atoi(argv[2]);
    acdm(Adj, edgeMap, iters, n, m, recalc, provided_best_loads);
}
