/**
 * @file record.cpp
 * @brief Implementation of the record class for handling text lines and word splitting.
 *
 * The record class manages a dynamically allocated character buffer (`line`)
 * and an array of words (`word`) extracted from that line. It supports
 * reading lines from files and splitting them into tokens.
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "record.h"

/**
 * @brief Construct a record with a line of random characters.
 *
 * Allocates a buffer of size @p numberOfChar. The buffer is initialized but
 * left unfilled; the terminating null character is set.
 *
 * @param numberOfChar Number of characters to allocate for the line.
 */
record::record( int numberOfChar ) {
  char_number = numberOfChar ;
  if( char_number > 0 ) {
    line = new char [char_number+1] ;
    line[char_number] = '\0' ;
  }else{
    line = NULL ;
  }
  words_number = 0 ;
  word = NULL ;
}



/**
 * @brief Construct a record from a C-style string.
 *
 * Allocates memory and copies the given string into the record's buffer.
 *
 * @param stringa Null-terminated C string.
 */
record::record( const char *stringa ) {
  char_number = 0 ;
  if( stringa != NULL ) {
    while( stringa[char_number] != '\0' ) char_number ++ ;
    line = new char [char_number+1] ;
    line[char_number] = '\0' ;
    for( int i=0 ; i< char_number ; i++ ) line[i] = stringa[i] ;
  }else{
    line = NULL ;
  }
  words_number = 0 ;
  word = NULL ;
}


/**
 * @brief Copy constructor.
 *
 * Creates a deep copy of another record, including line buffer and word array.
 *
 * @param other_rec Record to copy.
 */
record::record( const record &other_rec ) {
  char_number = other_rec.char_number ;
  if( char_number > 0 ) {
    line = new char [char_number+1] ;
    line = strcpy( line, other_rec.line ) ;
  }else{
    line = NULL ;
  }
  words_number = 0 ;
  word = NULL ;
}


/**
 * @brief Destructor.
 *
 * Frees allocated memory for the line and words.
 */
record::~record() {
  char_number = 0 ;
  delete[] line ;
  line = NULL ;
  for( int i=0; i<words_number; i++ ) {
    delete[] word[i] ;
    word[i] = NULL ;
  }
  words_number = 0 ;
  free( word ) ;
  word = NULL ;
}


/**
 * @brief Read a record (line) from an input file stream.
 *
 * Reads a line of text into the record's buffer, resizing if necessary.
 *
 * @param infile Input file stream.
 * @return true if a line was successfully read, false on EOF.
 */
bool record::getrecord( std::ifstream &infile ) {
  char c ;
  char *read_line = NULL ;
  delete[] line ;
  line = NULL ;
  char_number = 0 ;

  if( infile.eof() ) {
    return 0 ;
  }
  infile.get( c ) ;
  if( infile.eof() ) {
    return 0 ;
  }else{
    int i = 0 ;
    while( c != '\n' && (!infile.eof()) ) {
      if( i < char_number ) {
	read_line[i] = c ;
	i++ ;
	infile.get( c ) ;
      }else{
	char_number += 100 ;
	read_line = (char*)realloc( read_line, char_number*sizeof(char) ) ;
      }
    }
    char_number = i ;
    read_line = (char*)realloc( read_line, char_number*sizeof(char) ) ;
    line = new char [char_number+1] ;
    for( int j=0; j<char_number; j++ ) line[j] = read_line[j] ;
    line[char_number] = '\0' ;
    free( read_line ) ;
    read_line = NULL ;
    return 1 ;
  }
}


/**
 * @brief Split the line into words using whitespace as separator.
 *
 * Allocates the `word` array to store pointers to each token.
 * Overwrites existing `word` data if already allocated.
 */
void record::split( void ) {
  int word_length = 0, reading_word = 0 ;

  if( char_number > 0 ) {
    if( words_number > 0 ) {
      for( int i=0; i<words_number; i++ ) {
	delete[] word[i] ;
	word[i] = NULL ;
      }
      words_number = 0 ;
      free( word ) ;
      word = NULL ;
    }
    if( line[0] != ' ' && line[0] != '\t' ) words_number = 1 ;
    for( int i=1; i<char_number; i++ ) {
      if( ( line[i] != ' ' && line[i] != '\t' ) && ( line[i-1] == ' ' || line[i-1] == '\t' ) ) words_number++ ;
    }
    word = (char**)calloc( words_number, sizeof(char*) ) ;

    for( int i=0; i<char_number; i++ ) {
      if( line[i] != ' ' && line[i] != '\t' ) {
	for( word_length = 0; ( (word_length+i<char_number) && ( line[i+word_length]!=' ' && line[i+word_length]!='\t' ) ); word_length++ ) ;
	word[reading_word] = new char [word_length+1] ;
	for( int j=i; j<i+word_length; j++ ) word[reading_word][j-i] = line[j] ;
	word[reading_word][word_length] = '\0' ;
	reading_word ++ ;
	i = i + word_length - 1 ;
      }
    }
  }
}

/**
 * @brief Split the line into words using a custom separator character.
 *
 * Allocates the `word` array with dynamically allocated strings.
 * Overwrites existing `word` data if already allocated.
 *
 * @param separation Separator character (e.g., space, comma, tab).
 */
void record::split( char separator ) {
  int word_length = 0, reading_word = 0 ;

  if( char_number > 0 ) {
    if( words_number > 0 ) {
      for( int i=0; i<words_number; i++ ) {
	delete[] word[i] ;
	word[i] = NULL ;
      }
      words_number = 0 ;
      free( word ) ;
      word = NULL ;
    }
    if( line[0] != separator ) words_number = 1 ;
    for( int i=1; i<char_number; i++ ) {
      if( ( line[i] != separator ) && ( line[i-1] == separator ) ) words_number++ ;
    }
    word = (char**)calloc( words_number, sizeof(char*) ) ;

    for( int i=0; i<char_number; i++ ) {
      if( line[i] != separator ) {
	for( word_length = 0; ( (word_length+i<char_number) && ( line[i+word_length]!=separator ) ) ; word_length++ ) ;
	word[reading_word] = new char [word_length+1] ;
	for( int j=i ; j<i+word_length ; j++ ) word[reading_word][j-i] = line[j] ;
	word[reading_word][word_length] = '\0' ;
	reading_word ++ ;
	i = i + word_length - 1 ;
      }
    }
  }
}
