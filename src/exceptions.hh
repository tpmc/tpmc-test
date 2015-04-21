#ifndef TPMC_TEST_EXCEPTIONS_HH
#define TPMC_TEST_EXCEPTIONS_HH

#include <exception>
#include <string>

namespace tpmc_test
{
  class TpmcTestException : public std::exception
  {
  public:
    explicit TpmcTestException(const std::string& message)
        : message_(message)
    {
    }

    virtual const char* what() const throw() { return message_.c_str(); }
    virtual ~TpmcTestException() {}

  private:
    std::string message_;
  };

  struct IllegalArgumentException : public TpmcTestException {
    IllegalArgumentException(const std::string& message)
        : TpmcTestException(message)
    {
    }
  };

  struct NotImplementedException : public TpmcTestException {
    NotImplementedException(const std::string& message)
        : TpmcTestException(message)
    {
    }
  };
}

#endif // TPMC_TEST_EXCEPTIONS_HH
