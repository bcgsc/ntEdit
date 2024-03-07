#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/util.hpp"
#include <argparse/argparse.hpp>
#include <cmath>
#include <iostream>
#include <vector>

#if _OPENMP
#include <omp.h>
#endif

/*
Creating a Bloom filter from the provided genome sequence
*/

// Converting bits to bytes
const unsigned NUM_BITS_PER_BYTE = 8;

/*
Return the genome size of the given fasta file
*/
long long find_genome_size(std::vector<std::string> genome_files, int threads) {
  long long genome_size = 0;
  for (const auto& genome_file : genome_files) {
    btllib::SeqReader reader(
      genome_file, btllib::SeqReader::Flag::LONG_MODE, threads);
    for (const auto record : reader) {
      genome_size += record.seq.length();
    }
  }
  std::cout << "Genome size (bp): " << genome_size << std::endl;
  return genome_size;
}

/*
Approximate the required BF size based on the expected number of elements
Using the formula found in Broder & Mitzenmacher, 2004, 
code adapted from ntHits (https://github.com/bcgsc/ntHits)
*/
long long
get_bf_size(long long num_elements, double num_hashes, double fpr)
{
  double r = -num_hashes / log(1.0 - exp(log(fpr) / num_hashes));
  long long m = ceil(num_elements * r) / NUM_BITS_PER_BYTE;
  return m;
}

int
main(int argc, const char** argv)
{

  argparse::ArgumentParser parser("make_genome_bf");
  parser.add_argument("--genome")
    .nargs(argparse::nargs_pattern::at_least_one)
    .help("Input genome fasta file")
    .required();

  parser.add_argument("-k")
    .help("k-mer size (bp)")
    .required()
    .scan<'u', unsigned>();

  parser.add_argument("--fpr")
    .help("False positive rate for Bloom filter")
    .default_value((double)0.01)
    .scan<'g', double>();

  parser.add_argument("--hashes")
    .help("Number of hash functions")
    .default_value((unsigned)3)
    .scan<'u', unsigned>();

  parser.add_argument("-o")
    .help("Name for output Bloom filter")
    .default_value("genome_bf.bf");

  parser.add_argument("--bf")
    .help("Bloom filter size in bytes (optional)")
    .scan<'d', long long>();

  parser.add_argument("--num_elements")
    .help("Approximate number of elements for Bloom filter (used for calculating Bloom filter size, optional)")
    .scan<'d', long long>();

  parser.add_argument("-t")
    .help("Number of threads")
    .default_value(12U)
    .scan<'u', unsigned>();

  /* Parse the command-line arguments */
  try {
    parser.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << parser;
    std::exit(1);
  }

  std::vector<std::string> genome_files =
    parser.get<std::vector<std::string>>("genome");
  unsigned num_threads = parser.get<unsigned>("t");
  double fpr = parser.get<double>("fpr");
  unsigned k = parser.get<unsigned>("k");
  unsigned hashes = parser.get<unsigned>("hashes");
  std::string out_file = parser.get<std::string>("o");

  std::cout << "Parameters:" << std::endl;
  std::cout << "\t\t--genome ";
  for (const auto& genome : genome_files) {
    std::cout << genome << " ";
  }  
  std::cout << std::endl;
  std::cout << "\t\t-t " << num_threads << std::endl;
  std::cout << "\t\t-k " << k << std::endl;
  std::cout << "\t\t--fpr " << fpr << std::endl;
  std::cout << "\t\t--hashes " << hashes << std::endl;
  std::cout << "\t\t-o " << out_file << std::endl;

  #if _OPENMP
  omp_set_num_threads(num_threads);
  #endif

  /* Calculate BF size based on genome size */
long long bf_size;
  if (parser.is_used("--bf")) {
    bf_size = parser.get<long long>("bf");
    std::cout << "\t\t--bf " << bf_size << std::endl;
  } else if (parser.is_used("--num_elements")) {
    long long num_elements = parser.get<long long>("num_elements");
    std::cout <<"\t\t--num_elements " << num_elements << std::endl;
    bf_size = get_bf_size(num_elements, hashes, fpr);
  }
  else {
    std::cout << "Calculating BF size based on input genome size" << std::endl;
    long long genome_sum = find_genome_size(genome_files, num_threads);
    bf_size = get_bf_size(genome_sum, hashes, fpr);
  }
  std::cout << "BF size (bytes): " << bf_size << std::endl;


  /* Load the BF with k-mers from the genomes  */
  btllib::KmerBloomFilter* bf =
    new btllib::KmerBloomFilter(bf_size, hashes, k);

  for (const auto& genome : genome_files) {
      btllib::log_info("Reading " + genome);
      btllib::SeqReader reader(
          genome, btllib::SeqReader::Flag::LONG_MODE, num_threads);

      #pragma omp parallel
      for (const auto record : reader) {
        if (record.seq.length() >= k) {
          bf->insert(record.seq);
        }
    }
  }

  std::cout << "Bloom filter FPR: " << bf->get_fpr() << std::endl;

  btllib::log_info("Saving Bloom filter");
  bf->save(out_file);
  btllib::log_info("Done!");
  delete bf;
}