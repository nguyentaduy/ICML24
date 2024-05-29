#pragma GCC optimize("Ofast")
// #pragma GCC optimize("unroll-loops")
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

using namespace std;
#include "utilities.h"

vector<double> b_star;

void rcdm_permutation(vector<vector<pair<int, int>>> &Adj,
          vector<int> &edgeMap,
          long long iterations,
          long long n,
          long long m,
          int recalc,
          bool provided_best_loads = false)
{
  vector<double> z(2 * m), b(n);

  vector<int> perm;
  for (int i = 0; i < m; ++ i)
    perm.push_back(i);

  z = get_initial(Adj, n, m);
  for (int i = 0; i < n; ++i)
  {
    for (auto &j : Adj[i])
    {
      b[i] += z[j.second];
    }
  }

  double density = current_density(Adj, b, z, n, m);
  double density_sorting = 0.0;
  double tk = 1, tknew = 1;
  double fxk = 0;
  for (int i = 0; i < n; ++i)
    fxk += b[i] * b[i];

  for (int t = 1; t <= iterations; ++t)
  {
    auto start1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
    {
      b[i] = 0;
      for (auto &j : Adj[i])
      {
        b[i] += z[j.second];
      }
    }
    std::random_shuffle(perm.begin(), perm.end());
    for (int idx = 0; idx < m; ++idx)
    {
      int s = perm[idx];

      // pick a random edge
      int e = rand() % m;
      int eu = 2 * e;
      int u = edgeMap[eu];
      int ev = 2 * e + 1;
      int v = edgeMap[ev];

      double z_eu = z[eu] - b[u];
      double z_ev = z[ev] - b[v];

      if (abs(z_eu - z_ev) <= 1)
      {
        z_eu = (z_eu - z_ev + 1) / 2;
      }
      else if (z_eu - z_ev < -1)
      {
        z_eu = 0;
      }
      else
      {
        z_eu = 1;
      }
      z_ev = 1 - z_eu;
      b[u] = b[u] - z[eu] + z_eu;
      z[eu] = z_eu;
      b[v] = b[v] - z[ev] + z_ev;
      z[ev] = z_ev;
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
  bool print_diagnostics = true;
  int recalc = 1;
  if (argc > 2)
    recalc = atoi(argv[2]);
  rcdm_permutation(Adj, edgeMap, iters, n, m, recalc, provided_best_loads);
}
