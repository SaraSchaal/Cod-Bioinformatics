// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2017, Julian Catchen <jcatchen@illinois.edu>
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

#include "gzFasta.h"

#ifdef HAVE_LIBZ

GzFasta::GzFasta(const char *path) : Input(), gz_fh(NULL)
{
    this->open(path);
}

void
GzFasta::open(const char *path)
{
    this->gz_fh = gzopen(path, "rb");
    if (!this->gz_fh) {
        cerr << "Failed to open gzipped file '" << path << "': " << strerror(errno) << ".\n";
        exit(EXIT_FAILURE);
    }
    #if ZLIB_VERNUM >= 0x1240
    gzbuffer(this->gz_fh, libz_buffer_size);
    #endif
    int first = gzgetc(gz_fh);
    if (first == -1) {
        int errnum;
        const char* errstr = gzerror(gz_fh, &errnum);
        cerr << "Error: Failed to read any content from '" << path
            << "' (gzerror: " << errstr << ").\n";
        exit(EXIT_FAILURE);
    } else if (first != '>') {
        cerr << "Error: '" << path << "': not in fasta format (expected '>', got 0x"
             << std::hex << first << " '" << flush << (char) first << "')."
             << " (note: using zlib version " << zlibVersion()
             << "; compiled with zlib version " << ZLIB_VERSION << ")\n" << flush;
        exit(EXIT_FAILURE);
    }
    gzrewind(gz_fh);
}

Seq *
GzFasta::next_seq()
{
    #ifdef DEBUG
    DOES_NOT_HAPPEN; // As this function isn't efficient.
    #endif
    Seq* s = new Seq();
    s->id = new char[id_len];
    if(!next_seq(*s)) {
        delete s;
        return NULL;
    }
    return s;
}

int
GzFasta::next_seq(Seq &s)
{
    this->buf.clear();

    //
    // Check the contents of the line buffer. When we finish reading a FASTA record
    // the buffer will either contain whitespace or the header of the next FAST
    // record.
    //
    while (this->line[0] != '>' && !gzeof(this->gz_fh)) {
        gzgets(this->gz_fh, this->line, max_len);
    }

    if (gzeof(this->gz_fh)) {
        return false;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (len >= 1 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len >= 2 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    //
    // Check if the ID line of the FASTA file has a comment after the ID.
    //
    char* q = this->line + 1;
    ++q;
    while (*q != '\0' && *q != ' ' && *q != '\t')
        ++q;
    if (*q != '\0') {
        // Comment present.
        *q = '\0';
        ++q;
        s.comment.assign(q);
    }
    assert(s.id != NULL);
    strcpy(s.id, this->line + 1);

    //
    // Read the sequence from the file -- keep reading lines until we reach the next
    // record or the end of file.
    //
    gzgets(this->gz_fh, this->line, max_len);

    while (this->line[0] != '>' && !gzeof(this->gz_fh)) {
        len = strlen(this->line);
        if (len >= 1 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
        if (len >= 2 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

        this->buf    += this->line;
        this->line[0] = '\0';
        gzgets(this->gz_fh, this->line, max_len);
    }

    if (gzeof(this->gz_fh)) {
        len = strlen(this->line);
        if (len >= 1 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
        if (len >= 2 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

        this->buf += this->line;
        this->line[0] = '\0';
    }

    s.reserve(buf.length(), false);
    strcpy(s.seq, this->buf.c_str());

    return true;
}

#endif // HAVE_LIBZ
