#ifndef PSTRING_HH
#define PSTRING_HH

#include <string>
#include <vector>

/* A note about generalized indices:
 *   This class is meant to mimic Python's string functionality.  In particular,
 *   it uses what I call "generalized indices", meaning that any integer
 *   (positive or negative) is a valid index.  The mapping from a generalized
 *   index 'i' to an actual string position, for a string of length 'n', is:
 *     i -> 0           for i < -n,
 *     i -> i + n       for -n <= i < 0,
 *     i -> i           for 0 <= i < n,
 *     i -> n - 1       for i >= n.
 */

/**********************************************************************
 * pstring
 *
 * A wrapper around std::string, adding a lot of convenient functionality.
 * The class interface is modeled off Python's builtin string type.
 *********************************************************************/
class pstring : public std::string {
public:
    pstring();
    pstring(const std::string& s, size_type pos = 0, size_type n = npos);
    pstring(const char* s);
    pstring(const char* s, size_type n);
    pstring(char c, size_type n = 1);

    /* Assignment */
    pstring& operator=(const std::string& s);
    pstring& operator=(const char* s);
    pstring& operator=(char c);

    /* Concatenate strings */
    pstring& operator+=(const std::string& s);
    pstring& operator+=(const char* s);
    pstring& operator+=(char c);

    /* Return the character at the generalized index 'i' */
    char get(int i) const;

    /* Return the substring between generalized indices 'i' and 'j' */
    pstring slice(int i) const;
    pstring slice(int i, int j) const;

    /* Return the first or last 'n' characters of the string */
    pstring head(unsigned int n = 1) const;
    pstring tail(unsigned int n = 1) const;

    /* Strip whitespace from the ends of the string */
    pstring rstrip() const;
    pstring lstrip() const;
    pstring strip() const;
    void strip(bool);

    /* Strip quotes from a string */
    pstring stripquotes() const;
    void stripquotes(bool);

    /* Case conversion methods */
    pstring capitalize() const;
    pstring lower() const;
    pstring upper() const;

    /* Test whether the string contains the specified characters */
    bool contains(char c) const;
    bool contains(const char* s) const;
    bool contains(const std::string& s) const;

    /* Test whether the string starts or ends with the specified character or substring */
    bool startswith(char c) const;
    bool startswith(const char* s) const;
    bool startswith(const std::string& s) const;
    bool endswith(char c) const;
    bool endswith(const char* s) const;
    bool endswith(const std::string& s) const;

    /* Count the number of occurrences of the specified character or substring */
    int count(char c) const;
    int count(const char* c) const;
    int count(const std::string& s) const;

    /* Read up until the first newline character, starting from position 'm' */
    pstring firstline(unsigned int m = 0) const;

    /* Split the string into whitespace-separated components */
    std::vector<pstring> split() const;

    /* Split the string into components separated by 'sep' */
    std::vector<pstring> split(char sep, bool strip_results = true) const;

    /* Automatic conversion to C-style strings */
//    operator const char*() const;

    /* Numeric conversions */
    operator bool() const;
    operator long() const;
    operator double() const;
};

pstring str(bool b);
pstring str(int i);
pstring str(float f);
pstring str(double d);
pstring str(double d, int precision);

#endif // PSTRING_HH
