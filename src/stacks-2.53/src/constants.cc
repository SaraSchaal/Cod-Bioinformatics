#include <regex>
#include <cctype>

#include "constants.h"
#include "utils.h"

using namespace std;

int stacks_handle_exceptions(const exception& e) {
    std::cerr << "Aborted.";
    if (typeid(e) != typeid(std::exception)) {
        std::cerr << " (" << e.what();
        if (typeid(e) == typeid(std::bad_alloc))
            cerr << " -- did you run out of memory?";
        std::cerr << ")";
    }
    std::cerr << "\n";
    return 1;
}

const
map<string, FileT> known_extensions = {
    {".fa", FileT::fasta},
    {".fasta", FileT::fasta},
    {".fa.gz", FileT::gzfasta},
    {".fasta.gz", FileT::gzfasta},
    {".fq", FileT::fastq},
    {".fastq", FileT::fastq},
    {".fq.gz", FileT::gzfastq},
    {".fastq.gz", FileT::gzfastq},
    {".sam", FileT::sam},
    {".bam", FileT::bam},
    {".vcf", FileT::vcf},
    {".vcf.gz", FileT::gzvcf},
    {".calls", FileT::gzvcf}
};

string remove_suffix(FileT type, const string& orig) {
    string file (orig);

    int pos = file.find_last_of(".");

    if ((type == FileT::gzfastq || type == FileT::gzfasta) && file.substr(pos) == ".gz")
        file = file.substr(0, pos);

    pos = file.find_last_of(".");

    if (type == FileT::gzfastq || type == FileT::fastq) {

        if (file.substr(pos) == ".fastq" || file.substr(pos) == ".fq")
            file = file.substr(0, pos);

    } else if (type == FileT::gzfasta || type == FileT::fasta) {

        if (file.substr(pos) == ".fasta" || file.substr(pos) == ".fa")
            file = file.substr(0, pos);
    }

    return file;
}

regex init_file_ext_regex () {
    vector<string> exts;
    for (auto& filet : known_extensions)
        exts.push_back(filet.first);

    stringstream ss;
    ss << "(";
    join(exts, '|', ss);
    ss << ")$";

    string s = ss.str();
    escape_char('.', s);

    return regex(s);
}

FileT guess_file_type (const string& path) {

    static const regex reg = init_file_ext_regex();

    // Apply std::tolower.
    string copy = path;
    for (char& c : copy)
        c = tolower(c);

    smatch m;
    regex_search(copy, m, reg);

    if (m.empty())
        return FileT::unknown;
    else
        return known_extensions.at(m.str());
}

void escape_char(char c, string& s) {
    vector<size_t> dots;
    size_t i = -1;
    while ((i = s.find(c, i+1)) != string::npos)
        dots.push_back(i);

    for(auto j=dots.rbegin(); j!=dots.rend(); ++j)
        s.insert(*j, 1, '\\');
}
