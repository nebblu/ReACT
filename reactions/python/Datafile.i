%{
#include "Datafile.h"
%}

class Datafile {
public:
    Datafile();
    Datafile(const char* filename);
    virtual ~Datafile();

    int NumColumns() const;
    int NumRows() const;
    const char* GetFilename();

    bool Read(const char* filename);

    const char* GetHeader();
    const char* GetFooter();

    array GetColumn(int j);
    array GetRow(int i) const;

    bool Write(const char* filename);
    bool Write() const;

    void SetHeader(const char* s);
    void SetFooter(const char* s);

    void AddColumn(const array& Y);
    void InsertColumn(const array& Y, int j);
    void AddRow(const array& X);
    void InsertRow(const array& X, int i);
};
