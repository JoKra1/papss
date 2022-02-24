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
// 'updated terms' (see code below for a quick overview and Wood, 2011; Wood, 2017 for more details)
// obtained after repeated QR factorization and Cholesky decomposition, as described extensively
// in Wood (2011, 2017) to improve on the ill-conditioned nature of the former term.
void LambdaTerm::stepFellnerSchall(const Eigen::MatrixXd &embS, const Eigen::MatrixXd &cf, const Eigen::MatrixXd &inv,
                                   const Eigen::MatrixXd &gInv, int &cIndex, double sigma)
{

    // Embed penalties belonging to this lambda term in a zero-padded matrix embJ of size embS.
    int dimS = embS.rows();
    Eigen::MatrixXd embJ = Eigen::MatrixXd::Zero(dimS, dimS);
    this->embeddInS(embJ, cIndex);

    // Calculate Numerator and Demoninator of the FS update, as described above.
    double num = (gInv * embJ).trace() - (inv * embJ).trace();
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
void enforceConstraints(Eigen::VectorXd &cf, const Rcpp::StringVector &constraints)
{

    for (int idx = 0; i < constraints.size(); ++i)
    {
        std::string c = Rcpp::as<std::string>(constraints(i));
        // Equal comparison returns 0
        // See: https://www.cplusplus.com/reference/string/string/compare/
        if (c.compare("c") == 0)
        {
            if (cf(idx) < 0)
            {
                cf(idx) = 0;
            }
        }
    }
}

// Gradient descent optimizer with momentum and restarts and Paul Tseng's momentum update
// rule that allows to optimize a penalized least squares problem.
// Further permits for optimizing a penalized Non-negative least squares problem
// (penalized NNLS) by means of projection (see link below for details).
// For discussion of algorithm see: https://angms.science/doc/NMF/nnls_pgd.pdf
// code is based on the algorithm/pseudo-code discussed in aforementioned
// source and has been adapted to solve a penalized NNLS instead of a simple
// NNLS. Gradient of the penalized least squares loss function is from
// Wood (2017).
void agdTOptimize(Eigen::VectorXd &cf, const Eigen::MatrixXd &R, const Eigen::VectorXd &f,
                  const Eigen::MatrixXd &embS, const std::vector<char> &constraints, double r,
                  int maxiter, double tol)
{
    // Notation follows https://angms.science/doc/NMF/nnls_pgd.pdf.
    Eigen::MatrixXd Rt = R.transpose();
    Eigen::MatrixXd Q = Rt * R;
    Eigen::MatrixXd p = Rt * f;

    // Pre-calculate learning rate.
    Eigen::MatrixXd QQ = Q.array() * Q.array(); // Element wise Q_i**2
    double lr = 1 / sqrt(QQ.sum());

    // Add Penalty to R term.
    Q += embS;

    // Calculate momentum related coefficients.
    Eigen::VectorXd ycf0 = Eigen::VectorXd(cf);
    Eigen::VectorXd ycf = ycf0;
    Eigen::VectorXd prevCf = ycf0;

    // Error increase check.
    double prevErr = std::numeric_limits<double>::max();

    for (int i = 0; i < maxiter; ++i)
    {
        // Take an accelerated gradient step.
        cf = ycf - (lr * ((Q * ycf) - p));

        // Enforce constraints.
        enforceConstraints(cf, constraints);

        // Calculate momentum.
        double bk = (i) / (i + 3);
        ycf = cf + bk * (cf - prevCf);

        // Error calculation.
        Eigen::VectorXd err = f - R * cf;
        double errDot = err.dot(err);

        // Restart handling.
        if (errDot > prevErr)
        {
            cf = prevCf - (lr * ((Q * prevCf) - p));
            enforceConstraints(cf, constraints);
            ycf = cf;
        }

        // Crude convergence check.
        double absErrDiff = errDot > prevErr ? errDot - prevErr : prevErr - errDot;
        if (absErrDiff < tol)
        {
            break;
        }

        // Prepare next iter.
        prevCf = cf;
        prevErr = errDot;
    }
}

// Fits an additive model based on the stable LS solutions discussed in Wood (2011,2017). 
Eigen::MatrixXd solveAM(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::VectorXd> &initCf,
                        const Rcpp::StringVector &constraints, const std::vector<int> &lambdaTermFreq, int startIndex,
                        int maxIter, int maxIterOptim, double tol = 0.001)
{
    // Get dimension of X for re-use later.
    int rowsX = X.rows();
    int colsX = X.cols();

    // First compute QR decomposition of X.
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
    Eigen::MatrixXd R = qr.matrixR().template triangularView<Eigen::Upper>();

    // We do not need the zero rows in R (i.e., we care only about the "reduced/thin" QR decomposition).
    // Below, Eigen is forced to place the block in a temporary variable before assigning it to R.
    // See: http://eigen.tuxfamily.org/dox/TopicLazyEvaluation.html
    // See: https://stackoverflow.com/questions/30145771/shrink-matrix-with-eigen-using-block-in-assignment
    R = R.block(0, 0, colsX, colsX).eval();

    // We also do not need the full Q.
    // See: https://forum.kde.org/viewtopic.php?f=74&t=91271
    Eigen::MatrixXd Q = qr.householderQ().setLength(qr.nonzeroPivots()) * Eigen::MatrixXd::Identity(rowsX, colsX);

    // Extract permutation matrix applied by Householder QR decomp and apply it to R.
    // See: https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR.html
    // See: https://en.wikipedia.org/wiki/QR_decomposition#Column_pivoting
    Eigen::MatrixXd P = qr.colsPermutation();
    R *= P.transpose();

    // Now calculate f and r terms from Wood (2017)
    Eigen::VectorXd f = Q.transpose() * y;
    double r = y.dot(y) - f.dot(f);

    // Create a set of coefficients
    Eigen::VectorXd cf = Eigen::VectorXd(initCf);

    // Prepare lambda terms.
    std::vector<std::unique_ptr<LambdaTerm>> lambdaContainer;
    for (int freq : lambdaTermFreq)
    {
        std::unique_ptr<LambdaTerm> LDT = std::make_unique<LambdaTerm>(freq, colsX);
        lambdaContainer.push_back(std::move(LDT));
    }

    /*
    for(const std::unique_ptr<LambdaTerm> &LDT: lambdaContainer){
        const std::vector<std::unique_ptr<Penalty>> &penalties = (LDT)->getPenalties();
    }
    */

    // Error increase check.
    double prevErr = std::numeric_limits<double>::max();

    // Now iteratively optimize cf and then the current smoothness penalties.
    for (int i = 0; i < maxIter; ++i)
    {
        // Create embedded S term.
        int cInd = startIndex;
        Eigen::MatrixXd embS = Eigen::MatrixXd::Zero(colsX, colsX);

        // Now embed all lambda terms into S.
        for (const std::unique_ptr<LambdaTerm> &LDT : lambdaContainer)
        {
            (LDT)->embeddInS(embS, cInd);
        }

        // Optimize for cf given current lambda values.
        agdTOptimize(cf, R, f, embS, constraints, r, maxIterOptim, tol);

        // Now calculate the next step in the stable LS approach:
        // QR decomposition based on R + Cholesky factor of embS.

        // First the Cholesky factor using LDL' decomposition:
        // See: https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html
        Eigen::LDLT<Eigen::MatrixXd> ldlt(embS);

        // Now get the un-pivoted "Cholesky" factor of embS
        // Notation follows: https://services.math.duke.edu/~jdr/2021f-218/materials/week11.pdf
        // Start with un-pivoting the lower triangular matrix L.
        Eigen::MatrixXd L = ldlt.matrixL();
        Eigen::MatrixXd PL = ldlt.transpositionsP().transpose() * L;

        // The .abs() is a cheat here so that the root does not fail.
        // This is necessary since if embS is not positive semidefinite
        // D will contain very small negative elements.
        Eigen::VectorXd D = D.array().abs().sqrt();

        // Final factor
        Eigen::MatrixXd L1 = PL * D.asDiagonal();

        // Concatenate R and L1
        // See: https://stackoverflow.com/questions/21496157/eigen-how-to-concatenate-matrix-along-a-specific-dimension
        Eigen::MatrixXd RL1(colsX + L1.rows(), colsX);
        RL1 << R, L1;

        // Now form the next QR decomposition (See above).
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr2(RL1);
        Eigen::MatrixXd R2 = qr2.matrixR().template triangularView<Eigen::Upper>();
        R2 = R2.block(0, 0, colsX, colsX).eval();
        Eigen::MatrixXd P2 = qr2.colsPermutation();
        R2 *= P2.transpose();

        // Q here needs an additional update that was not necessary in the first step.
        Eigen::MatrixXd Q2 = qr2.householderQ().setLength(qr2.nonzeroPivots()) * Eigen::MatrixXd::Identity(RL1.rows(), colsX);
        Q2 = Q * Q2.block(0, 0, colsX, colsX).eval();

        // Now form the inverse of R2 (which is really the root of the inverse of X' %*% X + embS - see Wood (2017))
        Eigen::MatrixXd rInv = R2.inverse();

        // Now we can calculate the actual inverse, not the root, of the term above:
        Eigen::MatrixXd Inv = rInv * rInv.transpose();

        // We also need the pseudo inverse of embS for the EFS update.
        // See: https://eigen.tuxfamily.org/dox/classEigen_1_1CompleteOrthogonalDecomposition.html
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> ginvDecomp(embS);
        Eigen::MatrixXd gInv = ginvDecomp.pseudoInverse();

        // f and r also need to be recomputed.
        Eigen::VectorXd f2 = Q2.transpose() * y;
        double r2 = y.dot(y) - f2.dot(f2);

        // Finally we need to calculate the current estimate of sigma^2.
        Eigen::VectorXd res = f2 - R2 * cf;
        double errDot = res.dot(res) + r2;
        double sigma = errDot / (rowsX - (Inv * R.transpose() * R).trace());

        // Crude convergence control
        double absErrDiff = errDot > prevErr ? errDot - prevErr : prevErr - errDot;
        if (absErrDiff < tol)
        {
            break;
        }

        prevErr = errDot;

        // Now we can update all lamda terms.
        cInd = startIndex;
        for (const std::unique_ptr<LambdaTerm> &LDT : lambdaContainer)
        {
            // Update individual term.
            (LDT)->stepFellnerSchall(embS, cf, Inv, gInv, cInd, sigma);
        }
        
    }
    return cf;
}