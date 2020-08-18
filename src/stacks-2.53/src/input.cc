// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

//
// input.cc -- routines to read various formats of data into the XXX data structure.
//

#include "input.h"

Input::Input() {
    memset(this->line, '\0', max_len);
}

Input::Input(const char *path) {
    memset(this->line, '\0', max_len);

    this->path = string(path);
    //
    // Open the file for reading
    //
    this->fh.open(path, ifstream::in);

    if (this->fh.fail())
        cerr << "Error opening input file '" << path << "'\n";
}

Input::~Input() {
    // Close the file
    this->fh.close();
}

int
parse_tsv(const char *line, vector<string> &parts)
{
    const char  *p, *q;
    string part;

    parts.clear();
    p = line;

    do {
        for (q = p; *q != '\t' && *q != '\0'; q++);
        if (q - p == 0)
            part = "";
        else
            part.assign(p, (q - p));
        parts.push_back(part);

        p = q + 1;
    } while (*q != '\0');

    //for (size_t i = 0; i < parts.size(); i++)
    //    cerr << "Parts[" << i << "]: " << parts[i].c_str() << "\n";
    //cerr << "\n";

    return 0;
}

int
parse_ssv(const char *line, vector<string> &parts)
{
    const char  *p, *q;
    string part;

    parts.clear();
    p = line;

    do {
        for (q = p; *q != ' ' && *q != '\0'; q++);
        if (q - p == 0)
            part = "";
        else
            part.assign(p, (q - p));
        parts.push_back(string(part));

        p = q + 1;
    } while (*q != '\0');

    return 0;
}

int read_line(ifstream &fh, char **line, int *size) {
    char  buf[max_len];
    int   blen, llen;

    memset(*line, 0, *size);
    llen = 0;

    //
    // Make sure we read the entire line.
    //
    do {
        fh.clear();
        fh.getline(buf, max_len);

        blen = strlen(buf);
        if (blen + llen <= (*size) - 1) {
            strcat(*line, buf);
            llen += blen;
        } else {
            *size *= 2;
            llen  += blen;
            *line  = (char *) realloc(*line, *size);
            strcat(*line, buf);
        }
    } while (fh.fail() && !fh.bad() && !fh.eof());

    if (fh.eof() || fh.bad())
        return 0;

    return 1;
}

int read_gzip_line(gzFile &fh, char **line, int *size) {
    char  buf[max_len];
    int   blen, llen;
    bool  eol;

    memset(*line, 0, *size);
    llen = 0;
    eol  = false;

    //
    // Make sure we read the entire line.
    //
    do {
        if (gzgets(fh, buf, max_len) == NULL) break;

        blen = strlen(buf);

        if (blen > 0 && buf[blen - 1] == '\n') {
            eol = true;
            buf[blen - 1] = '\0';
        }

        if (blen + llen <= (*size) - 1) {
            strcat(*line, buf);
            llen += blen;
        } else {
            *size *= 2;
            llen  += blen;
            *line  = (char *) realloc(*line, *size);
            strcat(*line, buf);
        }
    } while (!gzeof(fh) && !eol);

    if (gzeof(fh))
        return 0;

    return 1;
}

bool
is_comment(const char *line)
{
    const char *p = line;

    while (*p != '\0')
        switch(*p) {
        case '#':
            return true;
            break;
        case ' ':
        case '\t':
            p++;
            break;
        default:
            return false;
            break;
        }

    return false;
}
