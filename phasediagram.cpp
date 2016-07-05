/* 
 * File:   phasediagram.cpp
 * Author: Abuenameh
 *
 * Created on 15 November 2013, 16:18
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <complex>
#include <deque>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <utility>
#include <map>
#include <cmath>
#include <algorithm>
#include <queue>

//using namespace std;

using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::deque;
using std::ostream;
using std::numeric_limits;
using std::setprecision;
using std::string;
using std::ostringstream;
using std::allocator;
using std::flush;
using std::pair;
using std::make_pair;
using std::max;
using std::map;
using std::ptr_fun;
using std::max_element;
using std::queue;

#define BOOST_THREAD_USE_LIB

#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/progress.hpp>
#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/date_time.hpp>
#include <boost/foreach.hpp>

using boost::lexical_cast;
using boost::multi_array;
using boost::extents;
using boost::progress_display;
using boost::random::mt19937;
using boost::random::uniform_real_distribution;
using boost::filesystem::path;
using boost::filesystem::exists;
using boost::filesystem::ofstream;
using boost::thread;
using boost::ref;
using boost::reference_wrapper;
using boost::bind;
using boost::mutex;
using boost::lock_guard;
using boost::thread_group;
using boost::posix_time::microsec_clock;
using boost::posix_time::ptime;
using boost::posix_time::time_period;

using namespace boost::algorithm;

#define USE_CXXLAPACK

//#include <flens/flens.cxx>
//
//using namespace flens;

#ifdef AMAZON
#include <lapacke.h>
extern "C" void openblas_set_num_threads(int num_threads);
#else
#include <clapack.h>
#endif

#include <armadillo>

using arma::vec;
using arma::mat;

//typedef IndexOptions<int, 0> ZeroBased;
//typedef GeMatrix<FullStorage<double, ColMajor> > DMatrix;
//typedef GeMatrix<FullStorage<complex<double>, ColMajor> > ZMatrix;
//typedef DenseVector<Array<double> > DVector;
//typedef DenseVector<Array<complex<double> > > ZVector;


#define ij(i, j) L * i + j

#include "mathematica.hpp"

template<typename T> void printMath(ostream& out, string name, T& t) {
    out << name << "=" << ::math(t) << ";" << endl;
}

template<typename T> void printMath(ostream& out, string name, int i, T& t) {
    out << name << "[" << i << "]" << "=" << ::math(t) << ";" << endl;
}

#define L 50
#define nmax 7
#define dim (nmax+1)
//#define I complex<double>(0, 1)

typedef mat::fixed<dim, dim> Operator;
typedef Operator Hamiltonian;
typedef vec::fixed<L> MeanField;
typedef vec::fixed<L> Observable;
typedef vec::fixed<L> Parameter;
typedef mat::fixed<L, L> SiteMatrix;
typedef vec::fixed<L> SiteVector;
typedef mat::fixed<L, L> HoppingMatrix;
typedef vec::fixed<dim> Eigenvalues;
typedef mat::fixed<dim, dim> Eigenvectors;
typedef vec::fixed<dim> SiteState;

typedef boost::array<double, L> ArrayL;
typedef boost::array<double, L*L> JMatrix;

typedef boost::array<double, dim> Diagonal;
typedef boost::array<double, dim*dim> LapackEigenvectors;
typedef boost::array<double, 2 * dim - 2 > Workspace;

inline int mod(int i) {
    return (i + L) % L;
    //    return ((i - 1 + L) % L) + 1;
}

mutex progress_mutex;

double chopeps = 0;

double chop(double x) {
    return (fabs(x) < chopeps) ? 0 : x;
}

double M = 1000;
double g13 = 2.5e9;
double g24 = 2.5e9;
double delta = 1.0e12;
double Delta = -2.0e10;
double alpha = 1.1e7;

double g = sqrt(M) * g13;

double JW(double W) {
    return alpha * (W * W) / (g * g + W * W);
}

double JWij(double Wi, double Wj) {
    return alpha * (Wi * Wj) / (sqrt(g * g + Wi * Wi) * sqrt(g * g + Wj * Wj));
}

ArrayL JW(ArrayL W) {
    ArrayL v;
    for (int i = 0; i < L; i++) {
        v[i] = W[i] / sqrt(g * g + W[i] * W[i]);
    }
    ArrayL J;
    for (int i = 0; i < L - 1; i++) {
        J[i] = alpha * v[i] * v[i + 1];
    }
    J[L - 1] = alpha * v[L - 1] * v[0];
    return J;
}

Parameter JW(Parameter W) {
    Parameter v = W / sqrt(g * g + W % W);
    Parameter J;
    for (int i = 0; i < L - 1; i++) {
        J[i] = alpha * v[i] * v[i + 1];
    }
    J[L - 1] = alpha * v[L - 1] * v[0];
    return J;
}

double UW(double W) {
    return -(g24 * g24) / Delta * (g * g * W * W) / ((g * g + W * W) * (g * g + W * W));
}

ArrayL UW(ArrayL W) {
    ArrayL U;
    for (int i = 0; i < L; i++) {
        U[i] = -(g24 * g24) / Delta * (g * g * W[i] * W[i]) / ((g * g + W[i] * W[i]) * (g * g + W[i] * W[i]));
    }
    return U;
}

Parameter UW(Parameter W) {
    return -(g24 * g24) / Delta * (g * g * W % W) / ((g * g + W % W) % (g * g + W % W));
}

mutex points_mutex;

struct Point {
    int i;
    int j;
    double x;
    double mu;
};

bool precalculate = false;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void phasepoints(Parameter& xi, double eps, queue<Point>& points, multi_array<double, 2 >& fc, multi_array<double, 2 >& fs, multi_array<MeanField, 2>& ares, progress_display& progress) {

    Operator Id, ni, bai, bci, bi;

    Id.eye();
    ni.zeros();
    bai.zeros();
    for (int n = 0; n <= nmax; n++) {
        ni(n, n) = n;
        if (n > 0) {
            bai(n - 1, n) = sqrt(n);
        }
    }
    bci = bai.t();
    bi = bai + bci;

    Operator n2i = ni * ni;
    Operator n2ni = n2i - ni;

    MeanField a, anew, da;
    Observable n;
    double Da;
    deque<MeanField> as, das;

    Hamiltonian H;
    Eigenvalues eigvals;
    Eigenvectors eigvecs;

    SiteMatrix A;
    SiteVector b;
    SiteVector phi;

    SiteMatrix rho;
    SiteVector evals;
    
	char job = 'V';
	int dimvar = dim;
	Diagonal D, E;
	LapackEigenvectors Z;
	int ldz = dim;
	Workspace work;
	int info;

        for (;;) {
        Point point;
        {
            mutex::scoped_lock lock(points_mutex);
            if (points.empty()) {
                break;
            }
            point = points.front();
            points.pop();
        }

        a.fill(0.1);

            double x = point.x;
            Parameter Wrand = x * xi;
            double Uav = 2 * UW(x); //1; //UW(Wrand).sum()/L;
            double Jav = mean(JW(Wrand)); //0;// = JW(Wrand).mean(); //sum()/L;
            Parameter U = UW(Wrand); // / Uav;
            U /= Uav;
            HoppingMatrix J;
            J(0, 0) = Jav / Uav; //J[i];
            for (int i = 0; i < L; i++) {
                int j1 = mod(i - 1);
                int j2 = mod(i + 1);
                J(i, j1) = JWij(Wrand[i], Wrand[j1]) / Uav; //J[i];
                J(i, j2) = JWij(Wrand[i], Wrand[j2]) / Uav; //J[i];
            }
            double mu = point.mu;

                    a.fill(0.1);
            anew.zeros();
            Da = numeric_limits<double>::infinity();
            //            das.clear();
            
            //        double Eij = 0;
            while (Da > eps) {
                //            Eij = 0;
                for (int i = 0; i < L; i++) {
                    int j1 = mod(i - 1);
                    int j2 = mod(i + 1);

                    double aj = J(i, j1) * a[j1] + J(i, j2) * a[j2];

				for (int n = 0; n <= nmax; n++) {
					D[n] = U[i] * n * (n - 1) - mu * n;
					if (n < nmax) {
						E[n] = aj * sqrt(n + 1);
					}
				}

#ifdef AMAZON
                    LAPACKE_dstev(LAPACK_COL_MAJOR, job, dim, D.data(), E.data(), Z.data(), ldz);
#else
				dstev_(&job, &dimvar, D.data(), E.data(), Z.data(), &ldz, work.data(), &info);
#endif
                                
				anew[i] = 0;
				for (int n = 0; n < nmax; n++) {
					anew[i] += sqrt(n + 1) * Z[n] * Z[n + 1];
				}

				n[i] = 0;
				for (int m = 0; m <= nmax; m++) {
					n[i] += m * Z[m] * Z[m];
				}

                }

                Da = max(abs(a) - abs(anew));
                a = anew;
            }
//            ares[point.i][point.j] = a;
            //        Eres[point.i][point.j] = Eij;

            double N = 0;
            for (int i = 0; i < L; i++) {
                N += n[i];
            }
            
            A.zeros();
            b.zeros();

            for (int j = 0; j < L; j++) {
                int j1 = mod(j - 1);
                int j2 = mod(j + 1);

                for (int i = 0; i < L; i++) {
                    if (i == j) {
                        A(j, i) = J(j, j1) * a[j1] + J(j, j2) * a[j2];
                    }
                    if ((i == j1) || (i == j2)) {
                        A(j, i) = -J(j, i) * a[i];
                    }
                }

                b[j] = J(j, j2) * a[j2] - J(j, j1) * a[j1];
            }
            phi = pinv(A) * b;

            double fsij = 0;
            for (int i = 0; i < L; i++) {
                int j1 = mod(i - 1);
                int j2 = mod(i + 1);

                fsij += J(i, j1) * a[j1] * a[i] * (phi(i) - phi(j1) + 1) * (phi(i) - phi(j1) + 1);
                fsij += J(i, j2) * a[j2] * a[i] * (phi(i) - phi(j2) - 1) * (phi(i) - phi(j2) - 1);
            }
            fsij /= 2 * J(0, 0) * N;

            fs[point.i][point.j] = fsij;

            rho.zeros();
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    if (i == j) {
                        rho(i, i) = n[i];
                    } else {
                        rho(i, j) = a[i] * a[j];
                    }
                }
            }
            rho /= N;

            evals = eig_sym(rho);
            fc[point.i][point.j] = evals[L - 1];

            //		cout << fc(point.i, point.j) << "\t" << fs(point.i, point.j) << endl;

            {
                mutex::scoped_lock lock(progress_mutex);
                ++progress;
            }

    }
}

/*
 *
 */
int main(int argc, char** argv) {
    
#ifdef AMAZON
    openblas_set_num_threads(1);
#endif

    mt19937 rng;
    uniform_real_distribution<> uni(-1, 1);

    int seed = lexical_cast<int>(argv[1]);
    int nseed = lexical_cast<int>(argv[2]);

    double xmin = lexical_cast<double>(argv[3]);
    double xmax = lexical_cast<double>(argv[4]);
    int nx = lexical_cast<int>(argv[5]);

    deque<double> x(nx);
    if (nx == 1) {
        x[0] = xmin;
    } else {
        double dx = (xmax - xmin) / (nx - 1);
        for (int ix = 0; ix < nx; ix++) {
            x[ix] = xmin + ix * dx;
        }
    }

    double mumin = lexical_cast<double>(argv[6]);
    double mumax = lexical_cast<double>(argv[7]);
    int nmu = lexical_cast<int>(argv[8]);

    deque<double> mu(nmu);
    if (nmu == 1) {
        mu[0] = mumin;
    } else {
        double dmu = (mumax - mumin) / (nmu - 1);
        for (int imu = 0; imu < nmu; imu++) {
            mu[imu] = mumin + imu * dmu;
        }
    }

    double D = lexical_cast<double>(argv[9]);
    double eps = lexical_cast<double>(argv[10]);

    int numthreads = lexical_cast<int>(argv[11]);

    int resi = lexical_cast<int>(argv[12]);

    precalculate = lexical_cast<bool>(argv[13]);

#ifdef AMAZON
    path resdir("/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/Gutzwiller Phase Diagram Additional");
#else
//    path resdir("/Users/Abuenameh/Dropbox/Amazon EC2/Simulation Results/Gutzwiller Phase Diagram");
    path resdir("/Users/Abuenameh/Documents/Simulation Results/Gutzwiller Phase Diagram Additional");
#endif
    if (!exists(resdir)) {
        cerr << "Results directory " << resdir << " does not exist!" << endl;
        exit(1);
    }
    for (int iseed = 0; iseed < nseed; iseed++, seed++) {
        ptime begin = microsec_clock::local_time();


        ostringstream oss;
        oss << "res." << resi << ".txt";
        path resfile = resdir / oss.str();
        while (exists(resfile)) {
            resi++;
            oss.str("");
            oss << "res." << resi << ".txt";
            resfile = resdir / oss.str();
        }
        if (seed < 0) {
            resi = seed;
            oss.str("");
            oss << "res." << resi << ".txt";
            resfile = resdir / oss.str();
        }

        Parameter xi;
        xi.fill(1);
        //        xi.assign(1);
        rng.seed(seed);
        if (seed > -1) {
            for (int j = 0; j < L; j++) {
                xi[j] = (1 + D * uni(rng));
            }
        }

        int Lres = L;

        ofstream os(resfile);
        printMath(os, "Lres", resi, Lres);
        printMath(os, "seed", resi, seed);
        printMath(os, "eps", resi, eps);
        printMath(os, "Delta", resi, D);
        printMath(os, "xres", resi, x);
        printMath(os, "mures", resi, mu);
        printMath(os, "xires", resi, xi);
        printMath(os, "precalc", resi, precalculate);
        os << flush;

        cout << "Res: " << resi << endl;

        multi_array<double, 2 > fcres(extents[nx][nmu]);
        multi_array<double, 2 > fsres(extents[nx][nmu]);
//        multi_array<double, 2> dur(extents[nx][nmu]);
//        multi_array<int, 2> iterres(extents[nx][nmu]);
//        multi_array<deque<MeanField>, 2> avecres(extents[1][1]);
//        multi_array<deque<MeanField>, 2> a2vecres(extents[1][1]);
        multi_array<MeanField, 2> ares(extents[1][1]);
//        multi_array<double, 2> Eres(extents[nx][nmu]);

        progress_display progress(nx * nmu);

        queue<Point> points;
        for (int imu = 0; imu < nmu; imu++) {
            queue<Point> rowpoints;
                    for (int ix = 0; ix < nx; ix++) {
                Point point;
                point.i = ix;
                point.j = imu;
                point.x = x[ix];
                point.mu = mu[imu];
                points.push(point);
            }
        }

        thread_group threads;
        for (int i = 0; i < numthreads; i++) {
            threads.create_thread(bind(&phasepoints, boost::ref(xi), eps, boost::ref(points), boost::ref(fcres), boost::ref(fsres), boost::ref(ares), boost::ref(progress)));
        }
        threads.join_all();


        printMath(os, "fcres", resi, fcres);
        printMath(os, "fsres", resi, fsres);
//                printMath(os, "dur", resi, dur);
//        printMath(os, "iters", resi, iterres);
//                printMath(os, "avecres", resi, avecres);
//                printMath(os, "a2vecres", resi, a2vecres);
//                printMath(os, "ares", resi, ares);
        //        printMath(os, "Eres", resi, Eres);

        ptime end = microsec_clock::local_time();
        time_period period(begin, end);
        cout << endl << period.length() << endl << endl;

        os << "runtime[" << resi << "]=\"" << period.length() << "\";" << endl;
    }

    //    cout << Eigen::SimdInstructionSetsInUse() << endl;

    return 0;
}

