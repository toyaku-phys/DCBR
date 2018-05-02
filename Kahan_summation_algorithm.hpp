#ifndef KSA_H
#define KSA_H
#include <boost/multiprecision/cpp_int.hpp>
namespace mp = boost::multiprecision;

class Kahan
{
   public:
      double sum;
      double c;
      mp::cpp_int counter;
      Kahan();
      Kahan(const double& val);
      Kahan& operator+=(const double& val);
      double get_av()const;
};

Kahan::Kahan()
{
   sum = 0.0;
   c   = 0.0;
   counter=0;
}

Kahan::Kahan(const double& val)
{
   sum = val;
   c   = 0.0;
   ++counter;
}

double Kahan::get_av()const
{
   const auto res = sum/(counter.convert_to<double>());
   return res;
}

#pragma clang optimize off
inline Kahan& Kahan::operator+=(const double& val)
{
   const auto y = val - c;
   const auto t = (this->sum) + y;
   this->c      = (t-(this->sum))-y;
   this->sum    = t;
   ++counter;
   return *this;
}
#pragma clang optimize on

#endif
