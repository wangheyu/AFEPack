/**
 * @file   Miscellaneous.cpp
 * @author Robert Lie
 * @date   Sun Nov  4 08:53:21 2007
 * 
 * @brief  
 * 
 * 
 */

#include "Miscellaneous.h"

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/char_traits.hpp>
#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/operations.hpp>

#define LIB_PATH "AFEPACK_TEMPLATE_PATH"

AFEPACK_OPEN_NAMESPACE

namespace {

  class comments_input_filter : public boost::iostreams::input_filter {
  public:
    explicit comments_input_filter(char comment_char = '#')
      : comment_char_(comment_char), skip_(false) {}

    template <typename Source>
    int get(Source& src)
    {
      int c;
      while (true) {
        if ((c = boost::iostreams::get(src)) == EOF || 
            (c == boost::iostreams::WOULD_BLOCK))
          break;
        skip_ = c == comment_char_ ?
          true :
          c == '\n' ?
          false :
          skip_;
        if (!skip_)
          break;
      }
      return c;
    }

    template <typename Source>
    void close(Source&) { skip_ = false; }
  private:
    char comment_char_;
    bool skip_;
  };

}

void OpenFilteredStream(const std::string& filename,
                        filtering_istream& is)
{
  is.push(comments_input_filter());
  is.push(boost::iostreams::file_source(filename.c_str()));
}




void OpenAFEPackLibraryFile(const std::string& file, 
                            filtering_istream& is)
{
  OpenFilteredStream(file, is);
  if (is.good()) {
    std::cerr << "AFEPack library file opened: "
              << file
              << std::endl;
  } else {
    std::cerr << "failed to open AFEPack library file: "
              << file
              << std::endl;
    abort();
  }
}

void LoadLibraryFunction(dlhandle_t& handle,
                         const std::string& sym,
                         dlhandle_t& fun_ptr)
{
  fun_ptr = dlsym(handle, sym.c_str());
  if (fun_ptr == NULL) {
    std::cerr << dlerror() << std::endl;
    abort();
  }
}

void ExpandString(std::string& str)
{
  wordexp_t result;
  switch (wordexp(str.c_str(), &result, 0)) {
  case 0: 
    break;
  case WRDE_NOSPACE: 
    wordfree(&result);
  default:
    std::cerr << "word expansion error." << std::endl;
    abort();
  };
  str = result.we_wordv[0];
  wordfree(&result);
}



void StringToWord(const std::string& str, const char& c, std::vector<std::string>& result)
{
  size_t i, j, n;
  result.clear();
  n = str.length();
  i = 0;
  do {
    j = str.find(c, i);
    if (j == std::string::npos) {
      result.push_back(str.substr(i, n-i));
      break;
    }
    result.push_back(str.substr(i, j-i));
    i = j + 1;
  } while (1);
}



void CombineString(const std::vector<std::string>& prefix, 
                   const std::vector<std::string>& suffix, 
                   std::vector<std::string>& result)
{
  size_t m = prefix.size();
  size_t n = suffix.size();
  result.resize(m*n);
  for (size_t i = 0;i < m;i ++)
    for (size_t j = 0;j < n;j ++)
      result[i*m + j] = prefix[i] + suffix[j];
}



std::string FindAFEPackLibraryFilePath(const std::string& filename)
{
  char * buffer = getenv(LIB_PATH);
  std::string lib_path;
  if (buffer == NULL)
    lib_path = ".";
  else {
    lib_path = buffer;
    lib_path += ":.";
  }
  std::vector<std::string> lib_path_vec;
  StringToWord(lib_path, ':', lib_path_vec);
  size_t n_path = lib_path_vec.size();
  size_t i = 0;
  for (;i < n_path;i ++) {
    std::string file = lib_path_vec[i] + "/" + filename;
    ExpandString(file);
    int filedes;
    if ((filedes = open(file.c_str(), O_RDONLY)) != -1) {
      close(filedes);
      std::cerr << "AFEPack library file found: "
                << file
                << std::endl;
      return lib_path_vec[i];
    }
  }
  if (i == n_path) {
    std::cerr << "AFEPack library file "
              << filename
              << " is not found in the following path:";
    for (i = 0;i < n_path;i ++) {
      std::string file = lib_path_vec[i];
      ExpandString(file);
      std::cerr << "\n\t" << file;
    }
    std::cerr << std::endl;
    abort();
  }
  // to remove compiling warning
  return "NULL";
}



dlhandle_t AFEPackDLOpen(const std::string& filename)
{
  if (filename.size() == 0) return NULL;

  std::string temp(filename);
  ExpandString(temp);
	
  std::cerr << "Opening shared library " << temp << " ..." << std::endl;
  dlhandle_t handle = dlopen(temp.c_str(), 1);
  if (handle == NULL) {
    std::cerr << "\ttried " << temp << ": failed" << std::endl;
  }
  else {
    std::cerr << "\ttried " << temp << ": success" << std::endl;
    return handle;
  }
	
  if (temp[temp.size()-3] == '.' && temp[temp.size()-2] == 'l'
      && temp[temp.size()-1] == 'a') { // I think this is a GNU Libtool shared library file
    std::ifstream is(temp.c_str());
    if (is.good()) {
      std::cerr << "\tI think this is a GNU libtool shared library while I can't open the file"
                << std::endl;
      abort();
    }
    char buffer[80];
    do {
      if (!is.good()) {
        std::cerr << "\tI think this is a GNU libtool shared library accoring the extension of the filename while it look not such a file"
                  << std::endl;
        abort();
      }
      is.getline(buffer, 80);
      if (strstr(buffer, "dlname")) {
        char * start = strchr(buffer, '\'');
        start ++;
        char * end = strrchr(buffer, '\'');
        int last = temp.size();
        int current = temp.rfind("/") + 1;
        temp.replace(current, last - current, start, (end - start)/sizeof(char));
        handle = dlopen(temp.c_str(), 1);
        if (handle == NULL) {
          std::cerr << "\ttried " << temp << ": failed" << std::endl;
          abort();
        }
        else {
          std::cerr << "\ttried " << temp << ": success" << std::endl;
          return handle;
        }
      }
    } while (1);
  }
	
  {
    int last = temp.size();
    int current = temp.rfind("so");
    temp.replace(current, last - current, "la");
    std::ifstream is(temp.c_str());
    if (!is) {
      std::cerr << "\ttried " << temp << ": can't open file" << std::endl;
      abort();
    }
    char buffer[80];
    do {
      if (!is.good()) {
        std::cerr << "\ttried " << temp << ": end of file" << std::endl;
        abort();
      }
      is.getline(buffer, 80);
      if (strstr(buffer, "dlname")) {
        char * start = strchr(buffer, '\'');
        start ++;
        char * end = strrchr(buffer, '\'');
        last = temp.size();
        current = temp.rfind("/") + 1;
        temp.replace(current, last - current, start, (end - start)/sizeof(char));
        handle = dlopen(temp.c_str(), 1);
        if (handle == NULL) {
          std::cerr << "\ttried " << temp << ": failed" << std::endl;
          abort();
        }
        else {
          std::cerr << "\ttried " << temp << ": success" << std::endl;
          return handle;
        }
      }
    } while (1);
  }
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

