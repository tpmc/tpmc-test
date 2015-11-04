#ifndef TPMC_TEST_IO_HH
#define TPMC_TEST_IO_HH

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <cstdlib>
#include "exceptions.hh"

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

  std::map<std::string,ini_value> parseIniFile(std::string fname)
  {
    auto lines = readFile(fname);
    std::map<std::string,ini_value> param;
    // parse content
    for (const auto & line : lines)
    {
      std::size_t pos;
      pos = line.find('=');
      if (pos != std::string::npos)
      {
        std::string key   = line.substr(0, pos);
        ini_value   value = line.substr(pos+1, line.size()-pos-1);
        param[trim(key)] = trim(value);
      }
    }
    return param;
  }

  std::pair<std::string,std::string> pathInfo(std::string fname)
  {
    std::string dir = "";
    std::size_t pos = fname.rfind('/');
    if (pos != std::string::npos)
    {
        dir   = fname.substr(0, pos+1);
        fname = fname.substr(pos+1, fname.size()-pos-1);
    }
    return std::pair<std::string,std::string>(dir,fname);
  }

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

  inline std::string parseCMDLineParameters(int argc, char** argv)
  {
    std::string application = tpmc_test::pathInfo(argv[0]).second;
    if (argc < 2
      or (std::string("basic") != argv[1] and std::string("full") != argv[1]))
    {
      std::cerr << "Error: no testtype provided.\nUsage:\t"
                << application << " [basic|full]\n";
      exit(-1);
    }
    return std::string(TOSTRING(TPMC_SRC_DIR)) +
      argv[1] + "_tests/" + application + ".ini";
  }

#undef STRINGIFY
#undef TOSTRING

}

#endif
