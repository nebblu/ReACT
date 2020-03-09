
# include <config.h>

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pstring.h"

using namespace std;

/* Static buffer for numeric conversions */
static char convbuf[64];

/* Clamp the index 'i' between 'min' and 'max' */
void clamp(int& i, int min, int max) {
    if(i < min)
        i = min;
    else if(i > max)
        i = max;
}


/*******************************************************************************
 * pstring
 ******************************************************************************/

pstring::pstring() {
}

pstring::pstring(const string& s, size_type pos, size_type n)
    : string(s, pos, n)
{
}

pstring::pstring(const char* s)
    : string(s)
{
}

pstring::pstring(const char* s, size_type n)
    : string(s, n)
{
}

pstring::pstring(char c, size_type n)
    : string(n, c)
{
}

pstring& pstring::operator=(const string& s) {
    this->assign(s);
    return *this;
}

pstring& pstring::operator=(const char* s) {
    this->assign(s);
    return *this;
}

pstring& pstring::operator=(char c) {
    this->assign(1, c);
    return *this;
}

pstring& pstring::operator+=(const string& s) {
    this->append(s);
    return *this;
}

pstring& pstring::operator+=(const char* s) {
    this->append(s);
    return *this;
}

pstring& pstring::operator+=(char c) {
    this->append(1, c);
    return *this;
}

char pstring::get(int i) const {
    if(i < 0)
        i += this->length();
    clamp(i, 0, this->length());
    return this->operator[](i);
}

pstring pstring::slice(int i) const {
    if(i < 0)
        i += this->length();
    clamp(i, 0, this->length());
    return pstring(*this, i, npos);
}

pstring pstring::slice(int i, int j) const {
    int len = this->length();
    if(i < 0)
        i += len;
    clamp(i, 0, this->length());

    if(j < 0)
        j += len;
    clamp(j, i, this->length());

    return pstring(*this, i, j - i);
}

pstring pstring::head(unsigned int n) const {
    return pstring(*this, 0, n);
}

pstring pstring::tail(unsigned int n) const {
    size_type len = this->length();
    if(n >= len)
        return *this;
    else
        return pstring(*this, len - n, npos);
}

pstring pstring::rstrip() const {
    size_type n = this->length();
    while(n > 0 && isspace((*this)[n - 1]))
        n--;

    return pstring(*this, 0, n);
}

pstring pstring::lstrip() const {
    size_type m = 0;
    while(m < this->length() && isspace((*this)[m]))
        m++;

    return pstring(*this, m, npos);
}

pstring pstring::strip() const {
    size_type n = length();
    while(n > 0 && isspace((*this)[n - 1]))
        n--;

    size_type m = 0;
    while(m < n && isspace((*this)[m]))
        m++;

    return pstring(*this, m, n - m);
}

pstring pstring::stripquotes() const {
    if(this->length() >= 2 && this->startswith('"') && this->endswith('"'))
        return this->slice(1, -1);
    else
        return *this;
}

void pstring::strip(bool) {
    size_type n = this->length();
    while(n > 0 && isspace((*this)[n - 1]))
        n--;
    if(n != this->length())
        this->erase(n, npos);

    size_type m = 0;
    while(m < n && isspace((*this)[m]))
        m++;
    if(m != 0)
        this->erase(0, m);
}

void pstring::stripquotes(bool) {
    *this = this->stripquotes();
}

pstring pstring::capitalize() const {
    pstring s = *this;
    if(!s.empty())
        s[(unsigned int)0] = toupper(s[(unsigned int)0]);
    return s;
}

pstring pstring::lower() const {
    pstring s = *this;
    for(unsigned int i = 0; i < s.length(); i++)
        s[i] = tolower(s[i]);
    return s;
}

pstring pstring::upper() const {
    pstring s = *this;
    for(unsigned int i = 0; i < s.length(); i++)
        s[i] = toupper(s[i]);
    return s;
}

bool pstring::contains(char c) const {
    return this->find(c) != npos;
}

bool pstring::contains(const char* s) const {
    return this->find(s) != npos;
}

bool pstring::contains(const string& s) const {
    return this->find(s) != npos;
}

bool pstring::startswith(char c) const {
    return !this->empty() && this->get(0) == c;
}

bool pstring::startswith(const char* s) const {
    return this->compare(0, strlen(s), s) == 0;
}

bool pstring::startswith(const string& s) const {
    return this->compare(0, s.length(), s) == 0;
}

bool pstring::endswith(char c) const {
    return !this->empty() && this->get(-1) == c;
}

bool pstring::endswith(const char* s) const {
    size_type n = strlen(s);
    return (n <= this->length()) && this->compare(this->length() - n, n, s) == 0;
}

bool pstring::endswith(const string& s) const {
    size_type n = s.length();
    return (n <= this->length()) && this->compare(this->length() - n, n, s) == 0;
}

int pstring::count(char c) const {
    int num = 0;
    for(size_type n = 0; (n = this->find(c, n)) != npos; n++)
        num++;
    return num;
}

int pstring::count(const char* s) const {
    int num = 0;
    for(size_type n = 0; (n = this->find(s, n)) != npos; n++)
        num++;
    return num;
}

int pstring::count(const string& s) const {
    int num = 0;
    for(size_type n = 0; (n = this->find(s, n)) != npos; n++)
        num++;
    return num;
}

pstring pstring::firstline(unsigned int m) const {
    if(m > this->length())
        return "";
    else
        return pstring(*this, m, this->find('\n', m) - m);
}

vector<pstring> pstring::split() const {
    vector<pstring> results;

    size_type i = 0, j = 0;
    while(true) {
        /* Skip whitespace */
        while(i < length() && isspace(at(i)))
            i++;
        if(i >= length())
            break;

        /* Read non-whitespace characters */
        j = i + 1;
        while(j < length() && !isspace(at(j)))
            j++;
        results.push_back(substr(i, j - i));

        i = j + 1;
    }
    return results;
}

vector<pstring> pstring::split(char sep, bool strip_results) const {
    vector<pstring> results;

    size_type i = 0, j = 0;
    while(true) {
        /* Find the next separator (or the end of the string) */
        while(j < length() && at(j) != sep)
            j++;

        /* Add the substring between the last separator and this one */
        if(strip_results)
            results.push_back(pstring(substr(i, j - i)).strip());
        else
            results.push_back(substr(i, j - i));

        if(j == length())
            break;
        i = ++j;
    }
    return results;
}

//pstring::operator const char*() const {
//    return this->c_str();
//}

pstring::operator bool() const {
    pstring small = this->strip().lower();
    if(small == "true")
        return true;
    else if(small == "false")
        return false;
    else {
        fprintf(stderr, "pstring::tobool: invalid boolean string: %s\n", small.c_str());
        return false;
    }
}

pstring::operator long() const {
    errno = 0;
    long l = strtol(this->c_str(), NULL, 10);
    if(errno != 0)
        fprintf(stderr, "pstring::tolong: error converting %s:\n  %s\n", this->c_str(), strerror(errno));
    return l;
}

pstring::operator double() const {
    errno = 0;
    double d = strtod(this->c_str(), NULL);
    if(errno != 0)
        fprintf(stderr, "pstring::todouble: error converting %s:\n  %s\n", this->c_str(), strerror(errno));
    return d;
}


/* Numeric conversion routines */

pstring str(bool b) {
    static pstring t = "true";
    static pstring f = "false";
    return b ? t : f;
}

pstring str(int i) {
    snprintf(convbuf, 64, "%d", i);
    return pstring(convbuf);
}

pstring str(float f) {
    snprintf(convbuf, 64, "%f", f);
    return pstring(convbuf);
}

pstring str(double d) {
    snprintf(convbuf, 64, "%f", d);
    return pstring(convbuf);
}

pstring str(double d, int precision) {
    snprintf(convbuf, 64, "%#.*g", precision, d);
    return pstring(convbuf);
}
