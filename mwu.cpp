

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

using namespace std;
#include "utilities.h"

vector<double> b_star;

void mwu(
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
  x = get_initial(Adj, n, m);
  for (int i = 0; i < n; ++i)
  {
    for (auto &j : Adj[i])
    {
      b[i] += x[j.second];
    }
  }
  // vector<int> indices(n);
  auto indices = new int[n];
  auto deleted = new int[n];
  for (int i = 0; i < n; ++i)
    deleted[i] = false;
  // vector<bool> deleted(n, false);
  iota(indices, indices + n, 0);
  double density = 1.0 * m / n;
  double density_sorting = 0.0;
  for (int t = 1; t <= iterations; ++t)
  {
    auto start = std::chrono::high_resolution_clock::now();
    sort(indices, indices + n, [&b](int &i, int &j) -> bool
         { return b[i] < b[j]; });
    // auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; ++i)
      deleted[i] = false;
    int M = m;
    int N = n;
    for (int k = 0; k < n; ++k)
    {
      int i = indices[k];
      int score_i = 0;
      for (auto &[j, e_idx] : Adj[i])
      {
        M -= !deleted[j];
        score_i += !deleted[j];
      }
      N--;
      b[i] += score_i;
      if (N)
        density = max(density, 1.0 * M / N);
      deleted[i] = true;
    }
    auto end = std::chrono::high_resolution_clock::now();

    // Manage xs, not part of algorithm, just calculated to see how well fractional peeling does
    for (int i = 0; i < n; ++i)
      deleted[i] = false;
    for (int k = 0; k < n; ++k)
    {
      int i = indices[k];
      for (auto &[j, e_idx] : Adj[i])
      {
        x[e_idx] += !deleted[j] * 1.0 / (t + 1);
      }
      deleted[i] = true;
    }

    cout << "Iteration=" << t << endl;
    cout << "Time=" << std::chrono::duration_cast<chrono::milliseconds>(end - start).count() << " miliseconds" << endl;
    cout << "Best density=" << density << endl;
    // double density_sorting = current_density_sorting(Adj, b, n, m);
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
  delete indices;
  delete deleted;
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
  mwu(Adj, iters, n, m, provided_best_loads);
}
