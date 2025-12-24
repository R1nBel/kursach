#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <string>

#include "RingMatrix.cpp"

RingMatrix random_matrix(int n, int mod, std::mt19937& rng)
{
    RingMatrix M(n, n, mod);
    std::uniform_int_distribution<int> dist(0, mod - 1);
    size_t total = (size_t)n * (size_t)n;
    for (size_t i = 0; i < total; ++i) 
    {
        M.data[i] = dist(rng);
    }
    return M;
}

bool verify_characteristic_polynomial(const RingMatrix& A, const int64_t* p)
{
    uint64_t n = A.rows;
    uint64_t mod = A.mod;

    std::vector<RingMatrix> powers;
    powers.reserve(n + 1);

    RingMatrix I(n, n, mod);
    for (int i = 0; i < n; ++i) I.data[(size_t)i * n + i] = 1 % mod;
    powers.emplace_back(I);

    for (int t = 1; t <= n; ++t) 
    {
        RingMatrix next = powers.back() * A;
        powers.emplace_back(next);
    }

    RingMatrix S = powers[n] * p[0];

    for (int k = 1; k <= n; ++k) 
    {
        RingMatrix term = powers[n - k] * p[k];
        S = S + term;
    }

    size_t total = (size_t)n * (size_t)n;
    for (size_t idx = 0; idx < total; ++idx) 
    {
        if (S.data[idx] % mod != 0) return false;
    }
    return true;
}

std::string timestamp_string()
{
    std::time_t now = std::time(nullptr);
    std::tm local_tm;
    localtime_s(&local_tm, &now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d_%H-%M-%S", &local_tm);
    return std::string(buf);
}

void run_experiments(
    int min_n,
    int max_n,
    int step_n,
    int mod,
    unsigned int seed,
    int quantity_per_n,
    const std::string& output_dir = ".")
{
    if (min_n <= 0 || max_n < min_n || step_n <= 0 || quantity_per_n <= 0)
    {
        std::cerr << "Invalid parameters for experiments." << std::endl;
        return;
    }

    std::string ts = timestamp_string();
    std::ostringstream fname;
    fname << output_dir << "experiments_danilevsky_" << ts << ".csv";

    std::ofstream csv(fname.str());
    if (!csv.is_open())
    {
        std::cerr << "Cannot open CSV file for writing: " << fname.str() << std::endl;
        return;
    }

    csv << "experiment_index;n;mod;avg_time_microseconds" << std::endl;

    std::mt19937 rng(seed);
    int experiment_index = 0;

    std::cout << "Danilevsky timing experiments (mod = " << mod << ")\nVerification is using only as experemential creteria. It is not a strict proof of the correctness of the work." << std::endl;
    std::cout << "n from " << min_n << " to " << max_n << " step " << step_n
        << ", each n: " << quantity_per_n << " matrices" << std::endl;
    std::cout << "Results will be written to " << fname.str() << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    for (int n = min_n; n <= max_n; n += step_n)
    {
        std::cout << "Processing n = " << n << " ..." << std::endl;
        long long sum_time_us = 0;
        int successful_tests = 0;

        for (int q = 0; q < quantity_per_n; ++q)
        {
            try 
            {
                RingMatrix A = random_matrix(n, mod, rng);

                auto t0 = std::chrono::high_resolution_clock::now();
                int64_t* coeffs = A.characteristicPolynomialDanilevsky();
                auto t1 = std::chrono::high_resolution_clock::now();

                bool ok = verify_characteristic_polynomial(A, coeffs);
                if (!ok)
                {
                    std::cerr << "Warning: verification failed for n=" << n << ", experiment " << q << std::endl;
                }

                delete[] coeffs;

                auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
                sum_time_us += (long long)dur.count();
                successful_tests++;

                std::cout << "  Experiment #" << experiment_index
                    << " n=" << n
                    << ", time=" << dur.count() << " micros" << std::endl;
            }
            catch (const RingMatrixException& e) 
            {
                std::cerr << "  Experiment #" << experiment_index
                    << " n=" << n << " skipped: " << e.what() << std::endl;
            }
            catch (const std::exception& e) 
            {
                std::cerr << "  Experiment #" << experiment_index
                    << " n=" << n << " skipped due to unexpected error: " << e.what() << std::endl;
            }

            ++experiment_index;
        }

        if (successful_tests > 0) 
        {
            long long avg_us = sum_time_us / successful_tests;

            int first_idx_for_n = experiment_index - quantity_per_n;
            csv << first_idx_for_n << ";" << n << ";" << mod << ";" << avg_us << std::endl;

            std::cout << "AVERAGE for n=" << n << " (based on " << successful_tests
                << " successful tests): " << avg_us << " us" << std::endl;
        }
        else 
        {
            std::cout << "No successful tests for n=" << n << std::endl;
        }
        std::cout << "------------------------------------------------" << std::endl;
    }

    csv.close();
    std::cout << "All experiments finished. CSV: " << fname.str() << std::endl;
}

void manual_test()
{
    try 
    {
     int n, mod;
        std::cout << "======MANUAL TEST=====\nEnter n : ";
        std::cin >> n;

        std::string matrix_string;

        std::cout << "Enter modulus: ";
        std::cin >> mod;

        matrix_string = std::to_string(mod) + " " +
            std::to_string(n) + " " +
            std::to_string(n) + " ";

        std::cout << "\nEnter " << n * n << " matrix elements (separated by spaces):\n";
        int element;
        for (int i = 0; i < n * n; i++)
        {
            std::cin >> element;
            matrix_string += std::to_string(element);
            if (i < n * n - 1)
            {
                matrix_string += " ";
            }
        }

        RingMatrix A = RingMatrix::fromString(matrix_string.c_str());

        std::cout << "Matrix A:" << std::endl;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                std::cout << A.at(i, j) << (j + 1 == n ? "" : " ");
            }
            std::cout << std::endl;
        }

        int64_t* coeffs = A.characteristicPolynomialDanilevsky();
        std::cout << "Characteristic polynomial coefficients (p[0]=1):" << std::endl;
        for (int k = 0; k <= n; ++k)
        {
            std::cout << coeffs[k] << (k == n ? "" : " ");
        }
        std::cout << std::endl;

        bool ok = verify_characteristic_polynomial(A, coeffs);
        std::cout << "Verification p(A) == 0 : " << (ok ? "OK" : "FAIL") << "\n(Verification is using only as experemential creteria. It is not a strict proof of the correctness of the work.)" << std::endl;

        delete[] coeffs;
    }

    catch (const RingMatrixException& e)
    {
        std::cerr << "Manual test skipped: " << e.what() << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Manual test skipped due to unexpected error: " << e.what() << std::endl;
    }
   
}

int main(int argc, char** argv)
{
    int min_n = 2;
    int max_n = 100;
    int step_n = 1;
    int mod = 1000003;
    unsigned int seed = 12345;
    int quantity_per_n = 5;
    std::string output_dir = "../experiments/";

    //manual_test();

    run_experiments(min_n, max_n, step_n, mod, seed, quantity_per_n, output_dir);

    return 0;
}