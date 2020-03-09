#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cassert>
#include <cstdio>

#include "Datafile.h"
#include "pstring.h"


Datafile::Datafile() {
    fname = "";
    ncol = 0;
    nrow = 0;
}

Datafile::Datafile(const char* filename) {
    fname = filename;
    ncol = 0;
    nrow = 0;
    Read(fname);
}

Datafile::Datafile(FILE* fp) {
    fname = "";
    ncol = 0;
    nrow = 0;
    Read(fp);
}

Datafile::~Datafile() {
}

bool Datafile::Read(const char* filename) {
    fname = filename;
    FILE* fp = fopen(fname, "r");
    return Read(fp);
}

bool Datafile::Read(FILE* fp) {
    header = footer = "";
    columns.clear();
    ncol = nrow = 0;

    if(!fp) {
        warning("Datafile: could not read from file %s\n", fname);
        return false;
    }

    char buf[1024];
    std::vector<pstring> fields;
    while(fgets(buf, sizeof(buf), fp) != NULL) {
        pstring line = buf;
        if(line.empty() || line.startswith('#')) {
            /* Add comment lines to header or footer, depending on whether or
             * not we've reached the first line of data */
            if(nrow == 0)
                header += line;
            else
                footer += line;
            continue;
        }

        /* Split line into fields */
        fields = line.split();
        if(ncol == 0) {
            ncol = (int)fields.size();
            columns.resize(ncol);
        }
        if((int)fields.size() != ncol) {
            warning("Datafile: expecting %d fields\n  line = '%s'", ncol, line.c_str());
            continue;
        }

        /* Add fields to column vectors */
        for(int j = 0; j < ncol; j++)
            columns[j].push_back((double)fields[j]);

        nrow++;
    }

    fclose(fp);
    return true;
}

const array& Datafile::GetColumn(int j) const {
    assert(1 <= j && j <= ncol);
    return columns[j-1];
}

array Datafile::GetRow(int i) const {
    assert(ncol > 0 && 1 <= i && i <= nrow);
    array X(ncol);
    for(int j = 0; j < ncol; j++)
        X[j] = columns[j][i-1];
    return X;
}

bool Datafile::Write(const char* filename) {
    fname = filename;
    return Write();
}

bool Datafile::Write() const {
    FILE* fp = fopen(fname, "w");
    return Write(fp);
}

bool Datafile::Write(FILE* fp) const {
    if(!fp) {
        warning("Datafile: could not write to file %s\n", fname);
        return false;
    }

    if(!header.empty()) {
        if(header.endswith('\n'))
            fprintf(fp, "%s", header.c_str());
        else
            fprintf(fp, "%s\n", header.c_str());
    }

    if(nrow > 0 && ncol > 0) {
        for(int i = 0; i < nrow; i++) {
            for(int j = 0; j < ncol-1; j++)
                fprintf(fp, "%e ", columns[j][i]);
            fprintf(fp, "%e\n", columns[ncol-1][i]);
        }
    }

    if(!footer.empty()) {
        if(footer.endswith('\n'))
            fprintf(fp, "%s", footer.c_str());
        else
            fprintf(fp, "%s\n", footer.c_str());
    }

    fclose(fp);
    return true;
}

void Datafile::AddColumn(const array& Y) {
    if(nrow != 0 && (int)Y.size() != nrow) {
        warning("Datafile: adding column of length %u to data file of size %dx%d\n", Y.size(), nrow, ncol);
        array Ytmp(Y);
        Ytmp.resize(nrow, 0);
        columns.push_back(Ytmp);
    }
    else {
        if(nrow == 0)
            nrow = (int)Y.size();
        columns.push_back(Y);
    }

    ncol++;
}

void Datafile::InsertColumn(const array& Y, int j) {
    if(nrow != 0 && (int)Y.size() != nrow) {
        warning("Datafile: inserting column of length %u in data array of size %dx%d\n", Y.size(), nrow, ncol);
        array Ytmp(Y);
        Ytmp.resize(nrow, 0);
        std::vector<array>::iterator it = columns.begin() + j - 1;
        columns.insert(it, Ytmp);
    }
    else {
        if(j < 0 || j > ncol)
            warning("Datafile: inserting column at position %d in data array of size %dx%d\n", j, nrow, ncol);
        if(j < 0) j = 0;
        if(j > ncol) j = ncol;
        if(nrow == 0)
            nrow = (int)Y.size();
        std::vector<array>::iterator it = columns.begin() + j - 1;
        columns.insert(it, Y);
    }

    ncol++;
}

void Datafile::AddRow(const array& X) {
    if(ncol != 0 && (int)X.size() != ncol)
        warning("Datafile: adding row of length %u to data file of size %dx%d\n", X.size(), nrow, ncol);

    if(ncol == 0)
        ncol = (int)X.size();
    for(int j = 0; j < ncol; j++)
        columns[j].push_back((j < (int)X.size()) ? X[j] : 0);

    nrow++;
}

void Datafile::InsertRow(const array& X, int i) {
    if(ncol != 0 && (int)X.size() != ncol)
        warning("Datafile: adding row of length %u to data file of size %dx%d\n", X.size(), nrow, ncol);
    if(i < 0 || i > nrow)
        warning("Datafile: inserting row at position %d in data array of size %dx%d\n", i, nrow, ncol);

    if(ncol == 0)
        ncol = (int)X.size();
    array::iterator it;
    for(int j = 0; j < ncol; j++) {
        it = columns[j].begin() + i-1;
        columns[j].insert(it, (j < (int)X.size()) ? X[j] : 0);
    }

    nrow++;
}
