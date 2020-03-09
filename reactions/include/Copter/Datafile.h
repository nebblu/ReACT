#ifndef DATAFILE_H
#define DATAFILE_H

#include <string>
#include <vector>

#include "Common.h"
#include "array.h"
#include "pstring.h"

/* Datafile
 *
 * A support class for reading data from file, or writing data to file.  A
 * "data set" in this context is a 2-D array of real numbers, arranged into
 * columns.  */

class Datafile {
public:
    /* Create an empty data set. */
    Datafile();

    /* Read data from file. */
    Datafile(const char* filename);

    /* Read data from an open file stream.  Closes stream when finished. */
    Datafile(FILE* fp);

    virtual ~Datafile();

    int NumColumns() const { return ncol; }
    int NumRows() const { return nrow; }
    const char* GetFilename() { return fname; }

    /***** Read mode *****/

    /* Read data from file. */
    bool Read(const char* filename);

    /* Read data from an open file stream.  Closes stream when finished. */
    bool Read(FILE* fp);

    /* Retrieve file metadata */
    const char* GetHeader() const { return header.c_str(); }
    const char* GetFooter() const { return footer.c_str(); }

    /* Retrieve data */
    const array& GetColumn(int j) const;
    array GetRow(int i) const;


    /***** Write mode *****/

    /* Write data to file. */
    bool Write(const char* filename);

    /* Write data to an open file stream.  Closes stream when finished. */
    bool Write(FILE* fp) const;

    /* Write data back to the same file that data was read from. */
    bool Write() const;

    /* Set file metadata */
    void SetHeader(const char* s) { header = s; }
    void SetFooter(const char* s) { footer = s; }

    /* Add Y as the last column of data. */
    void AddColumn(const array& Y);

    /* Insert Y as column j of the data set. */
    void InsertColumn(const array& Y, int j);

    /* Add X as the last row of data. */
    void AddRow(const array& X);

    /* Insert X as row i of the data set. */
    void InsertRow(const array& X, int i);

protected:
    const char* fname;
    pstring header;
    pstring footer;

    std::vector<array> columns;
    int ncol;
    int nrow;
};

#endif // DATAFILE_H
