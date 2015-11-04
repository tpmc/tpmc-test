#ifndef TPMC_TEST_UTILITY_HH
#define TPMC_TEST_UTILITY_HH

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <cstdlib>
#include "exceptions.hh"
#include <boost/filesystem.hpp>

namespace tpmc_test {

  // kill prefix / suffix whitespace
  inline std::string trim(std::string& str)
  {
    str.erase(0, str.find_first_not_of(' '));
    str.erase(str.find_last_not_of(' ')+1);
    return str;
  }

  std::vector<std::string> readFile(std::string fname)
  {
    std::ifstream file;
    file.open(fname);
    if (! file)
    {
      throw TpmcTestException("failed to open file \"" + fname + "\"");
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line))
      lines.push_back(line);
    return lines;
  }

  class ini_value : public std::string
  {
  public:
    ini_value() = default;
    ini_value(std::string value) : std::string(value) {}

    double to_double() const {
      return std::atof(std::string::c_str());
    }
    unsigned int to_uint() const {
      return std::atoi(std::string::c_str());
    }
    int to_int() const {
      return std::atoi(std::string::c_str());
    }
    std::vector<ini_value> to_vector() const {
      std::stringstream ss(*this);
      std::istream_iterator<std::string> begin(ss);
      std::istream_iterator<std::string> end;
      std::vector<ini_value> vec(begin, end);
      return vec;
    }
  };

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

  inline boost::filesystem::path parseCMDLineParameters(int argc, char** argv)
  {
    boost::filesystem::path application(argv[0]);
    application = application.filename();
    std::string testset = argv[1];
    std::string testdir = testset + "_tests";
    if (argc < 2
      or (std::string("basic") != testset and std::string("full") != testset))
    {
      std::cerr << "Error: no testtype provided.\nUsage:\t"
                << application << " [basic|full]\n";
      exit(-1);
    }
    return boost::filesystem::path(TOSTRING(TPMC_SRC_DIR)) / testdir / application.replace_extension(".ini");
  }

#undef STRINGIFY
#undef TOSTRING

}

#endif
