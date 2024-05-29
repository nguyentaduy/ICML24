

#pragma GCC optimize("Ofast")
// #pragma GCC optimize ("unroll-loops")
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
// #include <omp.h>

using namespace std;
#include "utilities.h"

vector<double> b_star;

void frankwolfe(
    vector<vector<pair<int, int>>> &Adj,
    long long iterations,
    long long n,
    long long m,
    bool provided_best_loads = false)
{

  long long Delta = max_element(Adj.begin(), Adj.end(), [](vector<pair<int, int>> &v1, vector<pair<int, int>> &v2) -> bool
                                { return v1.size() < v2.size(); })
                        ->size();

  vector<double> x(2 * m), b(n);
  // x = get_initial(Adj, n, m);
  for (int i = 0; i < 2 * m; ++i)
    x[i] = 0.5;
  vector<double> new_x(2 * m);
  for (int i = 0; i < n; ++i)
  {
    for (auto &j : Adj[i])
    {
      b[i] += x[j.second];
    }
  }

  double density = 1.0 * m / n;
  double density_sorting = 0.0;
  for (int t = 1; t <= iterations; ++t)
  {
    auto start = std::chrono::high_resolution_clock::now();
    double lr = 2.0 / (t + 2);
    for (int i = 0; i < n; ++i)
    {
      for (auto &[j, e_idx] : Adj[i])
      {
        int sister_idx = e_idx % 2 ? e_idx - 1 : e_idx + 1;
        if (b[i] < b[j])
        {
          new_x[e_idx] = 1;
          new_x[sister_idx] = 0;
        }
        else
        {
          new_x[e_idx] = 0;
          new_x[sister_idx] = 1;
        }
      }
    }
    for (int i = 0; i < 2 * m; ++i)
    {
      x[i] = (1 - lr) * x[i] + lr * new_x[i];
    }
    for (int i = 0; i < n; ++i)
    {
      for (auto &j : Adj[i])
      {
        b[i] += x[j.second];
      }
    }
    auto end = std::chrono::high_resolution_clock::now();
    density = max(density, current_density(Adj, b, x, n, m, t));
    cout << "Iteration=" << t << endl;
    cout << "Time=" << std::chrono::duration_cast<chrono::milliseconds>(end - start).count() << " miliseconds" << endl;
    cout << "Best density=" << density << endl;

    density_sorting = max(density_sorting, current_density_sorting(Adj, b, n, m));
    cout << "Sorting density=" << density_sorting << endl;
    if (provided_best_loads)
    {
      double E = 0;
      for (int i = 0; i < n; ++i)
      {
        E += (b[i] / t - b_star[i]) * (b[i] / t - b_star[i]);
      }
      E = sqrt(E);
      cout << "Error=" << E << endl;
    }

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
  // long long n = getInt(), m=getInt();
  long long n, m;
  cin >> n >> m;
  vector<vector<pair<int, int>>> Adj(n); // neighbor and edge index

  int e_idx = 0;
  for (int e = 0; e < m; ++e)
  {
    int i, j;
    cin >> i >> j;
    Adj[i].push_back({j, e_idx});
    Adj[j].push_back({i, e_idx + 1});
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
  if (argc > 2)
    recalc = atoi(argv[2]);
  frankwolfe(Adj, iters, n, m, provided_best_loads);
}
