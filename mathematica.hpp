/* 
 * File:   mathematica.hpp
 * Author: Abuenameh
 *
 * Created on 15 November 2013, 16:20
 */

#ifndef MATHEMATICA_HPP
#define	MATHEMATICA_HPP


#include <iostream>
#include <iomanip>
#include <deque>
#include <vector>
#include <map>
#include <string>
#include <complex>

using std::cout;
using std::string;
using std::pair;
using std::deque;
using std::ostream;
using std::ostringstream;
using std::numeric_limits;
using std::setprecision;
using std::complex;
using std::vector;

#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

using boost::multi_array;

using namespace boost::algorithm;

#include <armadillo>

using arma::vec;
using arma::mat;
using arma::cx_vec;
using arma::cx_mat;
using arma::cx_double;
using arma::uword;

template<typename T>
class mathematica {
public:

    mathematica(T& v_) : v(v_) {
    }

    T& v;
};

template<>
class mathematica<int> {
public:

    mathematica(int i_) : i(i_) {
    }

    int i;
};

template<>
class mathematica<double> {
public:

    mathematica(double d_) : d(d_) {
    }

    double d;
};

template<>
class mathematica<bool> {
public:

    mathematica(bool b_) : b(b_) {
    }

    bool b;
};

template<>
class mathematica<complex<double> > {
public:

    mathematica(complex<double> c_) : c(c_) {
    }

    complex<double> c;
};

ostream& operator<<(ostream& out, const mathematica<int> m) {
    out << m.i;
    return out;
}

ostream& operator<<(ostream& out, const mathematica<double> m) {
    double d = m.d;
    ostringstream oss;
    oss << setprecision(numeric_limits<double>::digits10) << d;
    out << replace_all_copy(oss.str(), "e", "*^");
    return out;
}

ostream& operator<<(ostream& out, const mathematica<bool> m) {
    out << (m.b ? "True" : "False");
    return out;
}

ostream& operator<<(ostream& out, const mathematica<complex<double> > m) {
    complex<double> c = m.c;
    out << "(" << mathematica<double>(c.real()) << ")+I(" << mathematica<double>(c.imag()) << ")";
    return out;
}

template<uword rows>
ostream& operator<<(ostream& out, const mathematica<vec::fixed<rows> >& m) {
    vec::fixed<rows>& v = m.v;
    out << "{" << mathematica<double > (v[0]);
    for (int i = 1; i < rows; i++) {
        out << "," << mathematica<double > (v[i]);
    }
    out << "}";
    return out;
}

template<uword rows>
ostream& operator<<(ostream& out, const mathematica<cx_vec::fixed<rows> >& m) {
    cx_vec::fixed<rows>& v = m.v;
    out << "{" << mathematica<cx_double > (v[0]);
    for (int i = 1; i < rows; i++) {
        out << "," << mathematica<cx_double > (v[i]);
    }
    out << "}";
    return out;
}

template<uword rows, uword cols>
ostream& operator<<(ostream& out, const mathematica<cx_mat::fixed<rows, cols> >& m) {
    cx_mat::fixed<rows, cols>& mat = m.v;
    out << "{{" << mathematica<cx_double > (mat(0, 0));
    for (int j = 1; j < cols; j++) {
        out << "," << mathematica<cx_double > (mat(0, j));
    }
    out << "}";
    for (int i = 1; i < rows; i++) {
        out << ",{" << mathematica<cx_double > (mat(i, 0));
        for (int j = 1; j < cols; j++) {
            out << "," << mathematica<cx_double > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}

template<uword rows, uword cols>
ostream& operator<<(ostream& out, const mathematica<mat::fixed<rows, cols> >& m) {
    mat::fixed<rows, cols>& mat = m.v;
    out << "{{" << mathematica<double > (mat(0, 0));
    for (int j = 1; j < cols; j++) {
        out << "," << mathematica<double > (mat(0, j));
    }
    out << "}";
    for (int i = 1; i < rows; i++) {
        out << ",{" << mathematica<double > (mat(i, 0));
        for (int j = 1; j < cols; j++) {
            out << "," << mathematica<double > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}

ostream& operator<<(ostream& out, const mathematica<vec>& m) {
    vec& v = m.v;
    out << "{" << mathematica<double> (v[0]);
    for (int i = 1; i < v.n_elem; i++) {
        out << "," << mathematica<double> (v[i]);
    }
    out << "}";
    return out;
}

ostream& operator<<(ostream& out, const mathematica<cx_vec>& m) {
    cx_vec& v = m.v;
    out << "{" << mathematica<cx_double > (v[0]);
    for (int i = 1; i < v.n_elem; i++) {
        out << "," << mathematica<cx_double > (v[i]);
    }
    out << "}";
    return out;
}

ostream& operator<<(ostream& out, const mathematica<mat>& m) {
    mat& mat = m.v;
    int rows = mat.n_rows;
    int cols = mat.n_cols;
    out << "{{" << mathematica<double > (mat(0, 0));
    for (int j = 1; j < cols; j++) {
        out << "," << mathematica<double > (mat(0, j));
    }
    out << "}";
    for (int i = 1; i < rows; i++) {
        out << ",{" << mathematica<double > (mat(i, 0));
        for (int j = 1; j < cols; j++) {
            out << "," << mathematica<double > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}

ostream& operator<<(ostream& out, const mathematica<cx_mat>& m) {
    cx_mat& mat = m.v;
    int rows = mat.n_rows;
    int cols = mat.n_cols;
    out << "{{" << mathematica<cx_double > (mat(0, 0));
    for (int j = 1; j < cols; j++) {
        out << "," << mathematica<cx_double > (mat(0, j));
    }
    out << "}";
    for (int i = 1; i < rows; i++) {
        out << ",{" << mathematica<cx_double > (mat(i, 0));
        for (int j = 1; j < cols; j++) {
            out << "," << mathematica<cx_double > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}

/*template<typename T>
ostream& operator<<(ostream& out, const mathematica<DenseVector<Array<T> > >& m) {
    typedef typename DenseVector<Array<T> >::IndexType IndexType;
    DenseVector<Array<T> > & v = m.v;
    IndexType first = v.firstIndex();
    IndexType last = v.lastIndex();
//    int length = v.length();
    out << "{" << mathematica<T > (v(first));
    for (int i = first + 1; i <= last; i++) {
        out << "," << mathematica<T > (v(i));
    }
    out << "}";
    return out;
}

template<typename T, StorageOrder Order>
ostream& operator<<(ostream& out, const mathematica<GeMatrix<FullStorage<T, Order> > > m) {
    typedef typename GeMatrix<FullStorage<T, Order> >::IndexType IndexType;
    GeMatrix<FullStorage<T, Order> >& mat = m.v;
    IndexType firstRow = mat.firstRow();
    IndexType firstCol = mat.firstCol();
    IndexType lastRow = mat.lastRow();
    IndexType lastCol = mat.lastCol();
//    IndexType r = mat.numRows();
//    IndexType c = mat.numCols();
    out << "{{" << mathematica<T > (mat(firstRow, firstCol));
    for (int j = firstCol + 1; j <= lastCol; j++) {
        out << "," << mathematica<T > (mat(firstRow, j));
    }
    out << "}";
    for (int i = firstRow + 1; i <= lastRow; i++) {
        out << ",{" << mathematica<T > (mat(i, firstCol));
        for (int j = firstCol + 1; j <= lastCol; j++) {
            out << "," << mathematica<T > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}*/

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<deque<T, Alloc> >& m) {
    deque<T, Alloc>& d = m.v;
    out << "{" << mathematica<T > (d[0]);
    for (int i = 1; i < d.size(); i++) {
        out << "," << mathematica<T > (d[i]);
    }
    out << "}";
    return out;
}

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<vector<T, Alloc> >& m) {
    vector<T, Alloc>& d = m.v;
    out << "{" << mathematica<T > (d[0]);
    for (int i = 1; i < d.size(); i++) {
        out << "," << mathematica<T > (d[i]);
    }
    out << "}";
    return out;
}

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<multi_array<T, 2, Alloc > >& m) {
    multi_array<T, 2, Alloc> & ma = m.v;
    int r = ma.shape()[0];
    int c = ma.shape()[1];
    out << "{{" << mathematica<T > (ma[0][0]);
    for (int j = 1; j < c; j++) {
        out << "," << mathematica<T > (ma[0][j]);
    }
    out << "}";
    for (int i = 1; i < r; i++) {
        out << ",{" << mathematica<T > (ma[i][0]);
        for (int j = 1; j < c; j++) {
            out << "," << mathematica<T > (ma[i][j]);
        }
        out << "}";
    }
    out << "}";
    return out;
}

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<multi_array<T, 3, Alloc > >& m) {
    multi_array<T, 3, Alloc> & ma = m.v;
    int r = ma.shape()[0];
    int c = ma.shape()[1];
    int d = ma.shape()[2];
    out << "{{{" << mathematica<T > (ma[0][0][0]);
    for (int k = 1; k < d; k++) {
        out << "," << mathematica<T > (ma[0][0][k]);
    }
    out << "}";
    for (int j = 1; j < c; j++) {
        out << ",{" << mathematica<T > (ma[0][j][0]);
        for (int k = 1; k < d; k++) {
            out << "," << mathematica<T > (ma[0][j][k]);
        }
        out << "}";
    }
    out << "}";
    for (int i = 1; i < r; i++) {
        out << ",{{" << mathematica<T > (ma[i][0][0]);
        for (int k = 1; k < d; k++) {
            out << "," << mathematica<T > (ma[i][0][k]);
        }
        out << "}";
        for (int j = 1; j < c; j++) {
            out << ",{" << mathematica<T > (ma[i][j][0]);
            for (int k = 1; k < d; k++) {
                out << "," << mathematica<T > (ma[i][j][k]);
            }
            out << "}";
        }
        out << "}";
    }
    out << "}";
    return out;

}

template<typename T, typename U>
ostream& operator<<(ostream& out, const mathematica<pair<T, U> >& p) {
    out << "{" << mathematica<T>(p.v.first) << "," << mathematica<U>(p.v.second) << "}";
    return out;
}

template<typename T, size_t N>
ostream& operator<<(ostream& out, const mathematica<boost::array<T, N> >& arr) {
    boost::array<T, N>& a = arr.v;
    out << "{" << mathematica<T>(a[0]);
    for (int i = 1; i < N; i++) {
        out << "," << mathematica<T>(a[i]);
    }
    out << "}";
    return out;
}

template<typename T>
mathematica<T> math(T& t) {
    return mathematica<T > (t);
}

mathematica<double> math(double d) {
    return mathematica<double>(d);
}

mathematica<complex<double> > math(complex<double> c) {
    return mathematica<complex<double> >(c);
}




#endif	/* MATHEMATICA_HPP */

