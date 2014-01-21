// Adapted from D0 Experiment jetcorr/src/RjetCorr.cpp

#include "ErrorTypes.hpp"

namespace jec {

  ErrorTypes::ErrorTypes(unsigned long int first, unsigned long int second) :
    _first(first), _second(second) {};  

  ErrorTypes::operator bool() const {
    return (_first || _second);
  }
  bool ErrorTypes::operator<(const ErrorTypes& rhs) const {
    if (this->_first==rhs._first) return (this->_second < rhs._second);
    else return (this->_first < rhs._first);
  }

  ErrorTypes operator~(const ErrorTypes& err) {
    return ErrorTypes(~err._first, ~err._second);
  }

  ErrorTypes operator&(const ErrorTypes& err1, const ErrorTypes& err2) {
    return ErrorTypes(err1._first & err2._first, err1._second & err2._second);
  }

  ErrorTypes operator|(const ErrorTypes& err1, const ErrorTypes& err2) {
    return ErrorTypes(err1._first | err2._first, err1._second | err2._second);
  }

  ErrorTypes operator^(const ErrorTypes& err1, const ErrorTypes& err2) {
    return ErrorTypes(err1._first ^ err2._first, err1._second ^ err2._second);
  }

}
