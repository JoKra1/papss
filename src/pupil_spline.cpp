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

// Penalty abstract class (interface)
// Is implemented by a specific penalty (probably first identity).
// Needs a public method to get the raw penalty matrix and the same matrix
// parameterized with a lambda term.
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