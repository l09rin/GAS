/**
 * @brief Rename a file by prompting the user for a new name.
 *
 * Uses `mv` command to rename the file at @p old_name.
 *
 * @param old_name Current file name/path.
 * @return Newly allocated string containing the new file name.
 */
char *rename(const char *old_name) {
  printf(" Insert a name for the file:_");
  fflush(0);
  char c, *new_name=NULL, *mv_command=NULL;
  int i=0, n_charNew=0, n_charOld=0, flag=0, equal=1;

  while(flag==0) {
    if(*(old_name+n_charOld)=='\0') {
      flag = 1;
    }else{
      n_charOld++;
    }
  }
  while((c=getchar())!='\n') {
    n_charNew++;
    new_name = (char*)realloc(new_name, n_charNew*sizeof(char));
    check_memalloc(new_name, "new_name in function rename()");
    *(new_name+n_charNew-1) = c;
  }
  *(new_name+n_charNew) = '\0';

  if(n_charOld==n_charNew) {
    for(i=0; i<n_charNew; i++) {
      if(*(old_name+i)!=(*(new_name+i))) equal = 0;
    }
  }else{
    equal = 0;
  }

  if(equal==0) {
    mv_command = (char*)calloc((9+n_charOld+n_charNew), sizeof(char));
    check_memalloc(mv_command, "mv_command in function rename()");
    sprintf(mv_command, "mv -f %s \"%s\"", old_name, new_name);
    flag = system(mv_command);
    free(mv_command);
  }else{
    flag = 0;
  }
  
  if(flag==0) {
    return new_name;
  }else{
    printf("*** ERROR : mv in function rename\n");
    return NULL;
  }
}
/*********************************************************************/

/**
 * @brief Ask the user a multiple-choice question.
 *
 * Displays the question and available options, waits for input,
 * and returns the selected character.
 *
 * @param question The prompt string.
 * @param options Valid options (string).
 * @param separator Character to separate options in prompt.
 * @return The chosen character (must be in @p options).
 */
char request(const char *question, const char *options, char separator) {
  char *opzioni, scelta='\0';
  int i=0, Nop=0, flag=0;
  opzioni = (char*)calloc(1, sizeof(char));
  check_memalloc(opzioni, "options in function request()");
  
  while((*(options+i))!='\0') {
    if((*(options+i))!=separator) {
      Nop++;
      opzioni = (char*)realloc(opzioni, Nop*sizeof(char));
      check_memalloc(opzioni, "options in function request()");
      *(opzioni+Nop-1) = *(options+i);
    }
    i++;
  }

  do{
    printf(" %s\n [", question);
    for(i=0; i<Nop; i++) {
      printf("%c/", (*(opzioni+i)));
    }
    printf("\b] ");
    i = scanf("%c", &scelta);                  /* Gli invio dati da tastiera dopo gli inserimenti restano nel buffer di lettura.              */
                                              /* scanf li ignora ma getchar no. Con %*c dopo il carattere letto lo elimino dal buffer (?).   */
    while(getchar() != '\n');   // Elimina il precedente invio dal buffer di lettura
    for(i=0; i<Nop; i++) {                    /* In scanf gli spazi significano ignora carattere?                                            */
      if(scelta==(*(opzioni+i))) flag=1;
    }
    if(flag==0) {
      scelta='\0';
    }
                                       
  }while(scelta=='\0');               
  free(opzioni);
  return scelta;
}
/*********************************************************************/


/**
 * @brief Check if a memory allocation was successful.
 *
 * Exits the program if the pointer is null.
 *
 * @tparam T Pointer type.
 * @param root Pointer returned by new/malloc.
 * @param point Debug string for identifying the allocation location.
 * @return 0 if allocation succeeded, never returns otherwise.
 */
template < class T >
inline int check_memalloc( T *root , const char *error ) {
    if ( root == NULL ) {
      printf( "*** ERROR : dynamic memory allocation failed !! -> %s\n" , error ) ; fflush(0) ;
      return 1 ;
    }else return 0 ;
}
/*********************************************************************/


/**
 * @brief Convert a string to a parameter of type T.
 *
 * Generic template version: errors out for unsupported types.
 *
 * @tparam T Type to convert to.
 * @param param_val Input string.
 * @param val Pointer to variable to store the converted value.
 */
template <typename T>
void get_param_from_string( const char *param_val , T *val ) {
  std::cout << "*** ERROR: Unknown type in function get_param_from_string()" << std::endl ;
  exit(EXIT_FAILURE) ;
}
template <>
void get_param_from_string( const char *param_val , double *val ) {
  sscanf( param_val , "%lf" , val ) ;
}
template <>
void get_param_from_string( const char *param_val , int *val ) {
  sscanf( param_val , "%d" , val ) ;
}
template <>
void get_param_from_string( const char *param_val , float *val ) {
  sscanf( param_val , "%f" , val ) ;
}
template <>
void get_param_from_string( const char *param_val , long *val ) {
  sscanf( param_val , "%ld" , val ) ;
}
template <>
void get_param_from_string( const char *param_val , unsigned int *val ) {
  sscanf( param_val , "%u" , val ) ;
}
template <>
void get_param_from_string( const char *param_val , char *val ) {
  sscanf( param_val , "%s" , val ) ;
}
/*********************************************************************/


/**
 * @brief Print the name and value of a parameter.
 *
 * @tparam T Parameter type.
 * @param param_name Name of the parameter (for printing).
 * @param param_val Pointer to the parameter value.
 */
template <typename T>
void print_parameter( const char *param_name , T *param_val ) {
  std::cout << " Parameter " << param_name << " : " << *param_val << std::endl ;
}
template <>
void print_parameter( const char *param_name , char *param_val ) {
  std::cout << " Parameter " << param_name << " : " << param_val << std::endl ;
}
/*********************************************************************/

/**
 * @brief Read a parameter from a configuration file.
 *
 * Searches for @p param_name in the input file @p input_script_name
 * and converts its value into type T.
 *
 * @tparam T Parameter type.
 * @param param_name Name of the parameter to search for.
 * @param input_script_name Path to the input configuration file.
 * @param param_val Pointer to the variable to store the result.
 */
template <typename T>
void read_parameter( const char *param_name , const char *input_script_name , T *param_val ) {
  std::ifstream config_file ;

  config_file.open( input_script_name ) ;
  if( config_file.fail() ) {
    std::cout << "\n  ERROR in function read_parameter() : parameters file not found!" << std::endl ;
    exit( EXIT_FAILURE ) ;
  }else{

    record *file_rec = NULL ;
    file_rec = new record ;
    while( file_rec->getrecord( config_file ) ) {
      file_rec->split() ;
      if( file_rec->words_number > 0 ) {
	if( file_rec->word[0][0] != '#' ) {
	  if( strcmp( file_rec->word[0], param_name ) == 0 ) get_param_from_string<T>( file_rec->word[1] , param_val ) ;
	}
      }
      delete file_rec ;
      file_rec = new record ;
    }
    delete file_rec ;
    file_rec = NULL ;
    config_file.close() ;
    print_parameter<T>( param_name , param_val ) ;
  }
}
/*********************************************************************/

/**
 * @brief Return the class name of an object given its pointer.
 *
 * Generic version returns "unknown". Specializations return
 * type names for supported classes.
 *
 * @tparam T Class type.
 * @param ptr Pointer to object.
 * @return Class name as string.
 */
template<class T>
const std::string ClassName(const T *ptr) {
  std::string s( "unknown" ) ;
  return s ;
}

template<>
const std::string ClassName(const particle_3D *ptr) {
  std::string s( "particle_3D" ) ;
  return s ;
}

template<>
const std::string ClassName(const particle_2D *ptr) {
  std::string s( "particle_2D" ) ;
  return s ;
}
/*********************************************************************/
