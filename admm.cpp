

#pragma GCC optimize("Ofast")
#pragma GCC optimize ("unroll-loops")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
 
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

////////////////////////////////////////////////////////////////////////////////////////
// Helper for fast input, From Flowless paper

inline char GET_CHAR(){
	const int maxn = 131072;
	static char buf[maxn],*p1=buf,*p2=buf;
	return p1==p2&&(p2=(p1=buf)+fread(buf,1,maxn,stdin),p1==p2)?EOF:*p1++;
}
inline int getInt() {
	int res(0);
	char c = GET_CHAR();
	while(c < '0') c = GET_CHAR();
	while(c >= '0') {
		res = res * 10 + (c - '0');
		c = GET_CHAR();
	}
	return res;
}

////////////////////////////////////////////////////////////////////////////////////////

vector<double> get_initial(vector<vector<pair<int, int>>>& Adj, long long n,  long long m){
    vector<double> x(2*m, 0.0);

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > pq;
    for(int i=0; i<n; ++i)pq.push({Adj[i].size(), i});
    vector<bool> deleted(n, false);
    vector<int> degree(n, 0);
    for(int i=0; i<n; ++i)degree[i]=Adj[i].size();
    long long N = n;
    long long M = m;
    double density = 1.0*M/N;
    while(!pq.empty()){
        auto [d, i] = pq.top();
        pq.pop();
        if (deleted[i])continue;
        
        for(auto & [j, e_idx] : Adj[i]){
            if ( !deleted[j] ){
                if (i!=j){
                    pq.push({--degree[j], j});
                    M--;
                }
                x[e_idx]=1;
            }
        }
        N--;
        if (N>0){
            density = max(density, 1.0*M/N);
        }
        deleted[i] = true;
    }
    // assert (abs(accumulate(x.begin(), x.end(), 0)-m)<=1e-5 );
    return x;
}

double current_density(vector<vector<pair<int, int>>>& Adj, vector<double>& b, vector<double>& x,  long long n, long long m, long long t=1){
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>> > pq;

    for(int i=0; i<n; ++i)pq.push({b[i]+Adj[i].size(), i});
    vector<bool> deleted(n, false);
    vector<double> degree(n, 0);
    for(int i=0; i<n; ++i)degree[i]=b[i]+Adj[i].size();
    long long N = n;
    long long M = m;
    double density = (double)(1.0*M/N);
    while(!pq.empty()){
        auto [d, i] = pq.top();
        pq.pop();
        if (deleted[i])continue;
        // cout << i << endl;
        for(auto & [j, e_idx] : Adj[i]){
            int sister_idx = e_idx % 2 ? e_idx-1 : e_idx + 1;
            if ( !deleted[j] ){
              pq.push({degree[j] - t * x[sister_idx]-1, j});
              degree[j] -= (t*x[sister_idx]+1);
              M--;
            }
        }
        N--;
        if (N>0){
            density = max(density, 1.0*M/N);
        }
        deleted[i] = true;
    }
    // assert(abs(M)<1e-5);
    return density;
}

vector<double> b_star;
void admm(vector<vector<pair<int, int>>>& Adj, 
            vector<tuple<int, int, bool, int>>& edges,  
            long long iterations, 
            long long n, 
            long long m, 
            bool provided_best_loads=false){
  long long Delta = max_element(Adj.begin(), Adj.end(), [](vector<pair<int, int>>& v1, vector<pair<int, int>>& v2)->bool {return v1.size()<v2.size();})->size();
  double density = 1.0*m/n;

  map<pair<int, int>, int> L_idx, z_idx, x_idx;
  int idx_one = 0;
  int idx_two = 0;
  for(int k=0; k<edges.size(); ++k){
      auto & [e_idx,sister_idx, bb, i] = edges[k];
      int j = get<3>(edges[sister_idx]);
      if (k%2==0){
          L_idx[{i,j}] = e_idx/2;
          L_idx[{j,i}] = e_idx/2;
      }
      z_idx[{i,j}] = e_idx;
      x_idx[{i,j}] = e_idx;
  }
  vector<double> b(n, 0), x(2*m, 0), L(m, 0), z(2*m, 0);
  x = get_initial(Adj, n, m);
  double c = 1000;

  for(int t = 1; t<=iterations; ++t){
    auto start = std::chrono::high_resolution_clock::now(); 
    for(int i=0; i<n; ++i){
      vector<double> best_x(Adj[i].size(), 0);
      vector<pair<double, int>> sorted_edges;
      for(int idx=0; idx<Adj[i].size(); ++idx){
        int j = Adj[i][idx].first;
        sorted_edges.push_back({-c*z[z_idx[{i,j}]]+L[L_idx[{i,j}]], idx});
      }
      sort(sorted_edges.begin(), sorted_edges.end());
      vector<double> ps;
      double ans = 0.0;
      for (auto & [p, j] : sorted_edges){
        ans += -p;
        ps.push_back(ans);
      }
      int best_j = -1;
      for(int j=0; j<sorted_edges.size(); ++j){
        double p = sorted_edges[j].first;
        if (-p - 2.0/(2.0*(j+1)+c)*ps[j] >=0){
          best_j = j;
        }
      }
      for(int j=0; j<=best_j; ++j){
        auto [p, real_j] = sorted_edges[j];
        best_x[real_j] = (-p - 2/(2*(best_j+1)+c)*ps[best_j])/c;
      }
      for(int ttt=0; ttt<Adj[i].size(); ++ttt){
        int j = Adj[i][ttt].first;
        x[x_idx[{i, j}]] = best_x[ttt];
      }
    }

    auto new_L = L;
    for(int i=0; i<n; ++i){
      for(auto & [j, e_idx] : Adj[i]){
        if (i<j){
          auto idx = L_idx[{i,j}];
          new_L[idx] = L[idx] + c/n * (x[x_idx[{i,j}]]+x[x_idx[{j,i}]]-1 );
        }
      }
    }

    double res_s = 0.0; //residual of dual 
    for(int i=0; i<n; ++i){
      double this_res = 0.0;
      for(auto & [j, e_idx] : Adj[i]){
        double new_z = x[x_idx[{i,j}]] + (L[L_idx[{i,j}]] - new_L[L_idx[{i,j}]])/c;
        this_res += c*(new_z - z[z_idx[{i,j}]]);
        z[z_idx[{i,j}]] = new_z;
      }
      res_s += this_res*this_res; 
    }
    res_s = sqrt(res_s);


    double res_x = 0;
    for(int i=0; i<n; ++i){
      for(auto & [j, e_idx] : Adj[i]){
        double temp = x[x_idx[{i,j}]] + x[x_idx[{j,i}]] - 1;
        res_x += (temp*temp); 
      }
    }
    res_x = sqrt(res_x);

    L = new_L;

    if (res_x > 10*res_s){
      cout << "Doubling penalty parameter" << endl;
      c*= 2;
    }
    else if (res_x < res_s/10){
      cout << "Halfing penalty parameter" << endl;
      c/=2;
    }

    auto end = std::chrono::high_resolution_clock::now();
    cout << "Iteration="<<t << endl;
    cout << "Time=" << std::chrono::duration_cast<chrono::milliseconds	>(end - start).count() << endl;
    for(int i=0; i<n; ++i){
      b[i]=0;
      for(auto & [j, e_idx] : Adj[i]){
        b[i] += x[x_idx[{i,j}]];
      }
      b[i] *= t;
    }
    density = max(density, current_density(Adj, b, x,  n, m, t));
    cout << "Best density=" << density << endl;
    if (provided_best_loads){
        double E = 0;
        for(int i=0; i<n; ++i){
          E += (b[i]/t - b_star[i])*(b[i]/t - b_star[i]);
        }
        E = sqrt(E);
        cout << "Error=" << E << endl;
    }
  }
}


int main(int argc, char** argv)
{
  ios_base::sync_with_stdio(0);
  cin.tie(0);
    
  cout << "Starting Graph reading" << endl;

  auto start = std::chrono::high_resolution_clock::now();
  // long long n = getInt(), m=getInt();
  long long n, m;
  cin >> n >> m;
  vector<vector<pair<int, int>>> Adj(n); //neighbor and edge index
  vector<tuple<int, int, bool, int>> edges; //index of edge, sister edge, and whether edge is primary

  int e_idx = 0;
  for(int e=0; e<m; ++e){
    int i, j;
    // i = getInt();
    // j = getInt();
    cin >> i >> j;
    Adj[i].push_back({j, e_idx});
    Adj[j].push_back({i, e_idx+1});
    edges.push_back({e_idx, e_idx+1, i<j, i});
    edges.push_back({e_idx+1, e_idx, j<i, j});
    e_idx+=2;
  }
  string line;
  getline(cin, line); // read \n 

  bool provided_best_loads = false;
  b_star = vector<double>(n, 0); //optimal load vector 
  int i = 0;
  while(getline(cin, line)){
    provided_best_loads = true;
    b_star[i++] = stod(line); 
  }
  auto end = std::chrono::high_resolution_clock::now();
  cout << "Time for reading input is " << std::chrono::duration_cast<chrono::milliseconds	>(end - start).count() << " ms" << endl;
  
  int iters = atoi(argv[1]);
  int recalc = 1;
  if (argc>2)recalc=atoi(argv[2]);
  admm(Adj, edges, iters, n, m, provided_best_loads);
}




