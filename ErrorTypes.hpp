// Adapted from D0 Experiment jetcorr/jetcorr/RjetCorr.hpp

#ifndef INC_ERRORTYPES_HPP
#define INC_ERRORTYPES_HPP

//might want to consider to use a std::bitset which is sufficiently long in order to avoid the first/second and 32 bit/64 bit confusion?
namespace jec {

  class ErrorTypes {
  public:
    ErrorTypes(unsigned long int first = 0L, unsigned long int second = 0L);
    ~ErrorTypes(){};    

    operator bool() const;
    bool operator<(const ErrorTypes& rhs) const;

    friend ErrorTypes operator~(const ErrorTypes& err);
    friend ErrorTypes operator&(const ErrorTypes& err1, const ErrorTypes& err2);
    friend ErrorTypes operator|(const ErrorTypes& err1, const ErrorTypes& err2);
    friend ErrorTypes operator^(const ErrorTypes& err1, const ErrorTypes& err2);

  private:
    // Each unsigned long int is 4 bytes or 32 bits long,
    // therefore the current implementation can hold 64 error bits
    // It is possible to add more internally without changing the interface
    unsigned long int _first;
    unsigned long int _second;
  };

  ErrorTypes operator~(const ErrorTypes& err);
  ErrorTypes operator&(const ErrorTypes& err1, const ErrorTypes& err2);
  ErrorTypes operator|(const ErrorTypes& err1, const ErrorTypes& err2);
  ErrorTypes operator^(const ErrorTypes& err1, const ErrorTypes& err2);

} /* namespace jes */

#endif /* INC_ERRORTYPES_HPP */
