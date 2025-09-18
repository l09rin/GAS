/**
 * @file record.h
 * @brief Declaration of the record class for reading and splitting text lines.
 */

#pragma once

#include <cstddef>
#include <fstream>

/**
 * @class record
 * @brief A text record that stores a line and its split words.
 *
 * The class manages a dynamically allocated character buffer (`line`)
 * and an array of words (`word`) extracted from that line.
 */
class record{
 public:
  /// Pointer to a dynamically allocated C-string representing the line.
  char *line = NULL ;
  /// Array of dynamically allocated C-strings representing the words.
  char **word = NULL ;
  /// Number of characters in the line (excluding null terminator).
  int char_number ;
  /// Number of words currently stored in the record.
  int words_number ;

  record( int numberOfChar = 0 ) ;
  record( const char *stringa ) ;
  record( const record &other_rec ) ;
  ~record() ;

  bool getrecord( std::ifstream &infile ) ;
  void split( void ) ;
  void split( char separation ) ;
};
