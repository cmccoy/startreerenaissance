#include <vector>
#include <iostream>

#include <Eigen/Core>


typedef Eigen::Matrix<double, 4, 4> DNAMatrix;
typedef Eigen::Array<double, 4, 1> DNAFreqArray;

struct Sequence
{
    DNAMatrix transitions;
    double distance;
};

struct GTRModel
{
    DNAMatrix q;
    DNAFreqArray pi;
};

int main()
{
    return 0;
}
