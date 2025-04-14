#ifndef UCR_PARSING_H
#define UCR_PARSING_H

/**
 * @file ucr_parsing.h
 * @brief Header file declaring methods to read in files from TSV format and specifically the format of UCR Time Series Archive.
 *
 */

#include <string>
#include <vector>

/**
 * @brief ucr_parsing namespace contains functions to parse datasets and file reading capabilities.
 */
namespace ucr_parsing {

  /**
   * @brief parse_tsv function takes a filename and a max_lines count (negative values will read the entire file).
   * @param filename is the name of the file, including filepath if not in present working directory and including file extension.
   * @param max_lines is the maximum number of lines to be read from the TSV file, negative will read all lines.
   * @return The series read from the file. 
   * The function skips any empty or NaN values in the file.
   */
  std::vector<double> parse_tsv(std::string filename, int max_lines);
  /**
   * @brief parse_ucr_tsv function takes a filename and a max_lines count (negative values will read the entire file).
   * @param filename is the name of the file, including filepath if not in present working directory and including file extension.
   * @param max_lines is the maximum number of lines to be read from the TSV file, negative will read all lines.
   * @return The series read from the file. 
   * The function skips any empty or NaN values in the file and skips the first entry of every line as this denotes the 'class' of the line in the UCR Time Series Archive.
   */
  std::vector<double> parse_ucr_tsv(std::string filename, int max_lines);

  /**
   * @brief DatasetType is an enum specific to the UCR Time Series Archive, datasets contain a TEST and TRAIN tsv, this enum specifies whether you want the series from TEST, TRAIN or TEST and TRAIN.
   */
  enum DatasetType { TEST, TRAIN, TRAIN_APPEND_TEST };

  /**
   * @brief parse_ucr_dataset function takes a datasets name, the location of all UCR datasets and which files to be read and returns the Time Series in the file(s).
   * @param dataset_name is the name of the dataset, the name of the folder the dataset is located in and the prefix of the name of the files.
   * @param dataset_loc is the folder location of the datasets from the present working directory, as it points to a folder, the string should end in /.
   * @param type is the files from the dataset that should be read, options TEST, TRAIN or TEST_APPEND_TRAIN specify the Time Series returned based on TEST and TRAIN present in directory.
   * @return The series read from the file.
   * The function skips any empty or NaN values in the file and skips the first entry of every line as this denotes the 'class' of the line in the UCR Time Series Archive.
   * Files inside the folder of the dataset should be of the format <dataset_name>_TEST.tsv or <dataset_name>_TRAIN.tsv.
   */
  std::vector<double> parse_ucr_dataset(std::string dataset_name, std::string dataset_loc, DatasetType type);

  /**
   * @brief parse_folder_names function takes a directory_path and returns the titles of any folders present at that location.
   * @param directory_path is the name of the file, including filepath if not in present working directory and including file extension.
   * @return The names of any folders in that directory.
   * The function skips any empty or NaN values in the file and skips the first entry of every line as this denotes the 'class' of the line in the UCR Time Series Archive.
   */
  std::vector<std::string> parse_folder_names(std::string directory_path);
};

#endif
