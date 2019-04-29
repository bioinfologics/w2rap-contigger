#include <iostream>
#include <fstream>
#include <unistd.h>
#include <SpectraCn.hpp>
#include "deps/cxxopts/cxxopts.hpp"

struct cmdline_args {
    std::string fasta_file;
    std::string kmers_file;
    std::string out_dir;
    std::string prefix;
};


int dir_exists_and_writable(const std::string &dirname) {
    return (access(dirname.data(), W_OK));
}

int file_exist (const std::string &filename)
{
    std::ifstream file(filename);
    return file.good();
}


struct cmdline_args parse_cmdline_args( int argc,  char* argv[]) {
    struct cmdline_args parsed_args;
    cxxopts::Options options(argv[0], "");
    try {
        options.add_options()
                ("f", "fasta", cxxopts::value(parsed_args.fasta_file))
                ("o,output_dir","output dir", cxxopts::value(parsed_args.out_dir))
                ("p,prefix","output prefix", cxxopts::value(parsed_args.prefix))
                ("h,help","show help message")
                ;

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("f")) {
            if (!file_exist(parsed_args.fasta_file)) {
                std::string error_str("Input file error, Fasta file: "+parsed_args.fasta_file+'\n');
                std::perror(error_str.data());
                throw std::runtime_error("Fasta error");
            }
        }

        if (result.count("output_dir")!=1 or result.count("prefix")!=1) {
            throw std::runtime_error("Please specify output dir and prefix");
        }

        // Validate directories exist
        if (dir_exists_and_writable(parsed_args.out_dir)) {
            std::string error_str("Output directory error: "+parsed_args.out_dir+", ");
            std::perror(error_str.data());
            throw std::runtime_error("Output directory error");
        }

        if (!file_exist(parsed_args.out_dir+"/raw_kmers.data")) {
            std::string error_str("Input file error, Kmers file: "+parsed_args.out_dir+"/raw_kmers.data"+'\n');
            std::perror(error_str.data());
            throw std::runtime_error("Kmers file error");
        }

        return parsed_args;

    } catch (const cxxopts::OptionException& e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }
}

/**
 * @brief
 * This program calculates a spectra-cn from a fasta file and the sorted output of the kmer counts of an assembly
 * Preconditions: raw_kmers.data has to be in sorted order
 */
int main(int argc, char **argv) {
    auto args=parse_cmdline_args(argc,argv);

    SpectraCN::DumpSpectraCN(args.fasta_file, args.out_dir, args.prefix);
    return 0;
}