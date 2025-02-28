#include "declarations.h"

// DSU (Disjoint Set Union) data structure
struct DSU {
    vector<int> parent, rank;
    DSU(int n) {
        parent.resize(n);
        rank.resize(n, 1);
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int x) {
        if (x != parent[x]) parent[x] = find(parent[x]);
        return parent[x];
    }
    void unite(int x, int y) {
        int xr = find(x), yr = find(y);
        if (xr != yr) {
            if (rank[xr] < rank[yr]) swap(xr, yr);
            parent[yr] = xr;
            if (rank[xr] == rank[yr]) rank[xr]++;
        }
    }
};

template <typename Func, typename... Args>
void measure_time(Func func, Args&&... args) {
    // Get the start time
    auto start = high_resolution_clock::now();

    // Call the function
    func(forward<Args>(args)...);

    // Get the end time
    auto end = high_resolution_clock::now();

    // Calculate the duration
    auto duration = duration_cast<microseconds>(end - start);

    // Output the time taken
    cout << "Time taken: " << duration.count() << " microseconds" << endl;
}

template<typename T>
void debug(T value) {
    cout << value << endl;
}

template<typename T, typename... Args>
void debug(T first, Args... args) {
    cout << first << " ";
    debug(args...);  // Recursive call with the remaining arguments
}

// Extended GCD
int extended_gcd(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int x1, y1;
    int d = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

// Modular inverse
int mod_inverse(int a, int m) {
    int x, y;
    int g = extended_gcd(a, m, x, y);
    return (g != 1) ? -1 : (x % m + m) % m;
}

// Matrix multiplication
matrix multiply(matrix& a, matrix& b) {
    int n = a.size();
    matrix result(n, vector<long long>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] = (result[i][j] + a[i][k] * b[k][j] % mod) % mod;
            }
        }
    }
    return result;
}

// Operator overload for matrix multiplication
matrix operator*(matrix& a, matrix& b) {
    return multiply(a, b);
}

// GCD and LCM
int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}
int lcm(int a, int b) {
    return a / gcd(a, b) * b;
}

// Reading functions
inline void read(vi& v, int n) {
    for (int i = 0; i < n; i++) {
        int h;
        cin >> h;
        v.push_back(h);
    }
}
inline void read(umi& m, int n) {
    for (int i = 0; i < n; i++) {
        int h;
        cin >> h;
        m[h]++;
    }
}
inline void read(mi& m, int n) {
    for (int i = 0; i < n; i++) {
        int h;
        cin >> h;
        m[h]++;
    }
}

// Prefix sum
inline void to_prefix_sum(vi& v) {
    v.insert(v.begin(), 0);
    int n = v.size();
    for (int i = 1; i < n; i++) {
        v[i] += v[i - 1];
    }
}

// Prime generator
vector<int> generate_primes(int n) {
    vector<bool> is_prime(n + 1, true);
    vector<int> primes;
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= n; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j <= n; j += i) {
                is_prime[j] = false;
            }
        }
    }
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

// Binomial coefficient
long long binomial_coefficient(int n, int k) {
    if (k > n) return 0;
    if (k > n - k) k = n - k;
    long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - i + 1) / i;
    }
    return result;
}

// Modular exponentiation
template<class T>
T mypow(T m1, T n) {
    T a = m1, result = 1;
    for (T i = n; i > 0; i >>= 1) {
        if (i & 1) {
            result = result * a % mod;
        }
        a = a * a % mod;
    }
    return result;
}


// Overload << operator for vector
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& vec) {
    os << "{ ";
    for (auto it = vec.begin(); it != vec.end(); ++it) {
        os << *it;
        if (next(it) != vec.end()) os << ", ";
    }
    os << " }";
    return os;
}

// Overload << operator for set
template<typename T>
ostream& operator<<(ostream& os, const set<T>& s) {
    os << "{ ";
    for (auto it = s.begin(); it != s.end(); ++it) {
        os << *it;
        if (next(it) != s.end()) os << ", ";
    }
    os << " }";
    return os;
}

// Overload << operator for map
template<typename K, typename V>
ostream& operator<<(ostream& os, const map<K, V>& m) {
    os << "{ ";
    for (auto it = m.begin(); it != m.end(); ++it) {
        os << "(" << it->first << " : " << it->second << ")";
        if (next(it) != m.end()) os << ", ";
    }
    os << " }";
    return os;
}

// Overload << operator for matrix (2D vector)
template<typename T>
ostream& operator<<(ostream& os, const vector<vector<T>>& matrix) {
    os << "{\n";
    for (const auto& row : matrix) {
        os << "  { ";
        for (auto it = row.begin(); it != row.end(); ++it) {
            os << *it;
            if (next(it) != row.end()) os << ", ";
        }
        os << " }\n";
    }
    os << "}";
    return os;
}

// Function to compute the LPS (Longest Prefix Suffix) array for the pattern
vector<int> compute_kmp_prefix(const string& pattern) {
    int m = pattern.size();
    vector<int> lps(m, 0);
    int len = 0;
    for (int i = 1; i < m; ++i) {
        while (len > 0 && pattern[i] != pattern[len]) {
            len = lps[len - 1];
        }
        if (pattern[i] == pattern[len]) {
            ++len;
        }
        lps[i] = len;
    }
    return lps;
}

// Function to perform KMP search and find all occurrences of the pattern in the text
vector<int> KMP_search(const string& text, const string& pattern) {
    vector<int> lps = compute_kmp_prefix(pattern);  // Step 1: Compute the LPS array
    vector<int> result;  // To store the indices of matches
    int n = text.size();
    int m = pattern.size();

    int i = 0, j = 0;  // i for text, j for pattern
    while (i < n) {
        if (text[i] == pattern[j]) {  // Match found
            i++;
            j++;
        }

        // If a full match of the pattern is found
        if (j == m) {
            result.push_back(i - j);  // Store the start index of the match
            j = lps[j - 1];  // Use the LPS array to skip unnecessary checks
        }
        // Mismatch after j matches
        else if (i < n && text[i] != pattern[j]) {
            // If j > 0, use the LPS array to shift the pattern
            if (j != 0) {
                j = lps[j - 1];
            }
            else {
                i++;
            }
        }
    }
    return result;
}

// Function to compute the LPS (Longest Prefix Suffix) array for the pattern
template<typename T>
vector<int> compute_kmp_prefix(const vector<T>& pattern) {
    int m = pattern.size();
    vector<int> lps(m, 0);
    int len = 0;
    for (int i = 1; i < m; ++i) {
        while (len > 0 && pattern[i] != pattern[len]) {
            len = lps[len - 1];
        }
        if (pattern[i] == pattern[len]) {
            ++len;
        }
        lps[i] = len;
    }
    return lps;
}

// Function to perform KMP search and find all occurrences of the pattern in the text
template<typename T>
vector<int> KMP_search(const vector<T>& text, const vector<T>& pattern) {
    vector<int> lps = compute_kmp_prefix(pattern);  // Step 1: Compute the LPS array
    vector<int> result;  // To store the indices of matches
    int n = text.size();
    int m = pattern.size();

    int i = 0, j = 0;  // i for text, j for pattern
    while (i < n) {
        if (text[i] == pattern[j]) {  // Match found
            i++;
            j++;
        }

        // If a full match of the pattern is found
        if (j == m) {
            result.push_back(i - j);  // Store the start index of the match
            j = lps[j - 1];  // Use the LPS array to skip unnecessary checks
        }
        // Mismatch after j matches
        else if (i < n && text[i] != pattern[j]) {
            // If j > 0, use the LPS array to shift the pattern
            if (j != 0) {
                j = lps[j - 1];
            }
            else {
                i++;
            }
        }
    }
    return result;
}

// Function to compute the LPS (Longest Prefix Suffix) array for the pattern
long long mod_exp(long long base, long long exp, long long mod) {
    long long result = 1;
    while (exp > 0) {
        if (exp % 2 == 1)
            result = (result * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}

// Miller-Rabin primality test
bool miller_rabin(long long n, int k = 40) {
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    long long d = n - 1;
    int r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<long long> dis(2, n - 2);

    for (int i = 0; i < k; i++) {
        long long a = dis(gen);
        long long x = mod_exp(a, d, n);
        if (x == 1 || x == n - 1) continue;

        for (int j = 0; j < r - 1; j++) {
            x = (x * x) % n;
            if (x == n - 1) break;
        }
        if (x != n - 1) return false;
    }

    return true;
}

long long next_prime(long long n) {
    while (!miller_rabin(n)) {
        n++;
    }
    return n;
}


