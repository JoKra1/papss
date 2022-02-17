// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <vector>
#include <memory>
#include <cmath>
#include <limits>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]


// ##################################### Classes #####################################


// Penalty abstract class (interface)
// Is implemented by a specific penalty (probably first identity).
// Needs a public method to get the raw penalty matrix and the same matrix
// parameterized with a lambda term.
// These penalty terms are discussed extensively in Wood (2017)
// 'Generalized Additive Models : An introduction with R, Second Edition'
class Penalty
{
protected:
    int dim;

public:
    virtual Eigen::MatrixXd getPenalty() = 0;
    virtual Eigen::MatrixXd parameterizePenalty(double l) = 0;
    int getDim();
};

// Dimension getter
int Penalty::getDim()
{
    return dim;
}

// Identity penalty class implements Penalty interface.
class IdentityPenalty : public Penalty
{
public:
    // Constructors
    IdentityPenalty(int dim);
    // Implement Penalty interface
    Eigen::MatrixXd getPenalty();
    Eigen::MatrixXd parameterizePenalty(double l);
};

// Constructor for IdentityPenalty.
// Parameters:
// int dim: dimensionalty (rows & cols) of the penalty
IdentityPenalty::IdentityPenalty(int dim)
{
    this->dim = dim;
}

// Returns identity penalty of dimensions (dim, dim)
Eigen::MatrixXd IdentityPenalty::getPenalty()
{
    Eigen::MatrixXd S = Eigen::MatrixXd::Identity(dim, dim);
    return S;
}

// Returns a lambda parameterized penalty of dimensions (dim, dim)
Eigen::MatrixXd IdentityPenalty::parameterizePenalty(double l)
{
    Eigen::MatrixXd Sl = getPenalty().array() * l;
    return Sl;
}

// LambdaTerm class. As discussed in Wood (2017), multiple penalties (for individual
// terms) can share the same lambda value (e.g., this is achieved with the 'by' keyword
// in mgcv). Thus, the LambdaTerm class stores smart pointers to all the penalties
// associated with this lambda value. It comes with a method to extract these penalties
// (Not sure whether that is needed, just thought it might be handy..), a method to embedd these
// penalties in a zero-padded matrix (S in Wood 2017) to represent them in quadratic form cf' * S * cf
// where cf contains a weight for each column in the model matrix (Wood, 2017 s. 4.3.1).
// This same method can also be used to embed multiple lambda terms in the same zero-padded matrix, which is
// required for the generalized Fellner Schall update described in Wood & Fasiolo (2017):
// 'A generalized Fellner-Schall method for smoothing parameter optimization with application
// to Tweedie location, scale and shape models'
// This update is also implemented as a method here so that the lambda values of each term can be updated.
class LambdaTerm
{
private:
    std::vector<std::unique_ptr<Penalty>> penalties; // Penalties associated with lambdaterm
    int nPenalties;                                  // Number of elements in penalties
    double lambda;

public:
    LambdaTerm(int n, int dim);
    const std::vector<std::unique_ptr<Penalty>> &getPenalties() const; // Will never change corresponding LambdaTerm instance.
    void embeddInS(Eigen::MatrixXd &embS, int &cIndex);
    void stepFellnerSchall(const Eigen::MatrixXd &embS, const Eigen::MatrixXd &cf, const Eigen::MatrixXd &inv,
                           const Eigen::MatrixXd &gInv, int &cIndex, double sigma);
};

// Associate n penalties of dimension (dim) with this LambdaTerm.
LambdaTerm::LambdaTerm(int n, int dim)
{
    nPenalties = n;
    lambda = 0.1;

    // Create associated penalty objects
    for (int i = 0; i < n; ++i)
    {
        /*
        Penalty is an abstract class, so if we later want to loop over the penalties
        in a vector, we cannot do for (Penalty S: penalties) {...} (1). Thus we store pointers
        to Penalties in the penalties vector (2,3). We need to move the unique pointer to the penalties,
        since it cannot be copied (2,3). Later we can then call the virtual functions required by Penalty
        implementations through the pointer (4).

        (1) See: https://en.cppreference.com/w/cpp/language/abstract_class
        (2) See: https://en.cppreference.com/w/cpp/memory/unique_ptr
        (3) See: https://stackoverflow.com/questions/31410858/adding-elements-to-stdvector-of-an-abstract-class
        (4) See: https://en.cppreference.com/w/cpp/language/virtual
        */
        std::unique_ptr<Penalty> newPenalty = std::make_unique<IdentityPenalty>(dim);
        penalties.push_back(std::move(newPenalty));
    }
}

// Return a const reference to the penalties. Method is itself declared as const
// since it should never modfiy an instance which calls this.
// See: https://isocpp.org/wiki/faq/const-correctness
const std::vector<std::unique_ptr<Penalty>> &LambdaTerm::getPenalties() const
{
    return penalties;
}

// Embed all penalties belonging to this term in a zero-padded matrix embS as discussed by Wood (2017, s. 4.3.1).
// cIndex corresponds to the starting index (row & column) at which we want to start embedding the penalties.
// This is handy in case the model matrix contains un-penalized terms. Thus, any model matrix used by this implementation
// should be ordered - starting with unpenalized columns/terms following by penalized terms.
void LambdaTerm::embeddInS(Eigen::MatrixXd &embS, int &cIndex)
{
    /*
    As described in the constructor of the LambdaTerm class, we here loop over the
    ptrs in the penaltyList (i.e., again over there references) and access the
    virtual functions implemented by the specific penalties to fill the full matrix block.
    */
    for (const std::unique_ptr<Penalty> &S : penalties)
    {
        int dimS = S->getDim();
        embS.block(cIndex, cIndex, dimS, dimS) = S->parameterizePenalty(lambda);
        cIndex += dimS;
    }
}

// Perform a generalized Fellner Schall update step for a lambda term. This update rule is
// discussed in Wood & Fasiolo (2017). In the paper, the authors provide the update in terms
// of X and y (as well as X and z for generalized models). However, as discussed in Wood (2017)
// and Wood (2011): 'Fast stable restricted maximum likelihood and marginal likelihood estimation
// of semiparametric generalized linear models: Estimation of Semiparametric Generalized Linear Models'
// the explicit calculation of (X' * X + embS)^-1 is undesirable. Thus we here invoke the update on
// 'updated terms' (see below for a quick overview and Wood, 2011; Wood, 2017 for more details)
// obtained after repeated QR factorization and Cholesky decomposition, as described extensively
// in Wood (2011, 2017) to improve on the ill-conditioned nature of the former term. Wood (2017) shows
// that X can first be decomposed into X = Q * R. By obtaining f = t(Q) * y, the normal
// solution to a least squares problem cf = (X' * X)^-1 * X' * y can be re-expressed in terms of
// R and f as cf = (R)^-1 * f, allowing to consider the least squares problem purely in terms of R and f.
// In case of a penalized least squares problem (as the one considered here), involving a penalty
// matrix like embS, Wood (2011, 2017) shows that after first obtaining the Cholesky factor cholS of
// embS = cholS' * cholS and then forming a second QR decomposition of the row-wise concatenation [R,cholS] = Q2 * R2,
// R2 corresponds to the root of X' * X + embS - thus yielding (X' * X + embS)^-1 = (R2)^-1 * (R2)^-1'
// (Note that additional updates to f are necessary, involving Q and Q2, see: Wood, 2017).
// Thus, the update rule as discussed in Wood & Fasiolo (2017):
// lambda_i+1 = sigma^2 * ((tr(ginv(embS) * embJ) - tr((X' * X + embS)^-1) * embJ) / (cf' * embJ * cf)) * lambda_i
// is here calculated as:
// lambda_i+1 = sigma^2 * ((tr(ginv(embS) * embJ) - tr((R2)^-1 * (R2)^-1' * embJ) / (cf' * embJ * cf)) * lambda_i
// Where tr() denotes the trace calculation and ginv() corresponds to forming the generalized MP. inverse.
// Similarly, instead of calculating sigma^2 as:
// sigma^2 = ((y - X * cf).(y - X * cf)) / (n - tr((X' * X + embS)^-1 * X` * X))
// We here calculate sigma^2 as:
// sigma^2 = ((f - R * cf).(f - R * cf) + ((y).(y) - (f).(f))) / (n - tr((R2)^-1 * (R2)^-1' * R' * R))
// Here n is the number of rows in X and ().() refers to the dot product between the vectors in () and ().
void LambdaTerm::stepFellnerSchall(const Eigen::MatrixXd &embS, const Eigen::MatrixXd &cf, const Eigen::MatrixXd &inv,
                                   const Eigen::MatrixXd &pInv, int &cIndex, double sigma)
{

    // Embed penalties belonging to this lambda term in a zero-padded matrix embJ of size embS.
    int dimS = embS.rows();
    Eigen::MatrixXd embJ = Eigen::MatrixXd::Zero(dimS, dimS);
    this->embeddInS(embJ, cIndex);


    // Calculate Numerator and Demoninator of the FS update, as described above.
    double num = (pInv * embJ).trace() - (inv * embJ).trace();
    Eigen::MatrixXd denom = cf.transpose() * embJ * cf;

    // Now calculate the lambda update for this term.
    lambda = sigma * num / denom(0, 0) * lambda;
}


// ##################################### Functions #####################################


// Enforces positivity constraints on a vector. Based on the work by
// Hoeks & Levelt (1993): 'Pupillary dilation as a measure of attention: A quantitative system analysis'
// we require all 'attention spikes', i.e., the weights in our cf vector to be positive.
// Thus, we unfortunately cannot rely on a closed solution for our optimization problem but have to
// resort to perform projected gradient optimization, as discussed here: https://angms.science/doc/NMF/nnls_pgd.pdf.
void enforceConstraints(Eigen::VectorXd &cf, const std::vector<char> &constraints)
{
    int idx = 0;
    for (char c : constraints)
    {
        if (c == 'c')
        {
            if (cf(idx) < 0)
            {
                cf(idx) = 0;
            }
        }
        ++idx;
    }
}
