/* Written 2012 by Matthias S. Benkmann
 *
 * The author hereby waives all copyright and related rights to the contents
 * of this example file (example.cpp) to the extent possible under the law.
 */

/**
 * @file
 * @brief Small demo of The Lean Mean C++ Option Parser.
 *
 * @include example.cpp
 */

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <vector>

#include "../../optionparser/src/optionparser.h"

int
example( int argc, char const * argv[] )
{
  enum optionIndex
  {
    UNKNOWN,
    HELP,
    PLUS
  };
  const option::Descriptor usage[] = {
    { UNKNOWN,
      0,
      "",
      "",
      option::Arg::None,
      "USAGE: example [options]\n\n"
      "Options:" },
    { HELP,
      0,
      "",
      "help",
      option::Arg::None,
      "  --help  \tPrint usage and exit." },
    { PLUS,
      0,
      "p",
      "plus",
      option::Arg::None,
      "  --plus, -p  \tIncrement count." },
    { UNKNOWN,
      0,
      "",
      "",
      option::Arg::None,
      "\nExamples:\n"
      "  example --unknown -- --this_is_no_option\n"
      "  example -unk --plus -ppp file1 file2\n" },
    { 0, 0, 0, 0, 0, 0 } };

  argc -= ( argc > 0 );
  argv += ( argc > 0 );  // skip program name argv[0] if present
  option::Stats stats( usage, argc, argv );
  std::vector< option::Option > options( stats.options_max );
  std::vector< option::Option > buffer( stats.buffer_max );
  option::Parser parse( usage, argc, argv, &options[0], &buffer[0] );

  if( parse.error() ) return 1;

  if( options[HELP] || argc == 0 )
  {
    option::printUsage( std::cout, usage );
    return 0;
  }

  std::cout << "--plus count: " << options[PLUS].count() << "\n";

  for( option::Option * opt = options[UNKNOWN]; opt; opt = opt->next() )
    std::cout << "Unknown option: " << std::string( opt->name, opt->namelen )
              << "\n";

  for( int i = 0; i < parse.nonOptionsCount(); ++i )
    std::cout << "Non-option #" << i << ": " << parse.nonOption( i ) << "\n";
  return 0;
}

TEST( testOptionParser, example )
{
  int const argc = 3;
  char const * argv[] = { "testOptionParser", "--plus", "--plus" };
  example( argc, argv );
}

struct Arg : public option::Arg
{
  static void
  printError( const char * msg1,
              const option::Option & opt,
              const char * msg2 )
  {
    fprintf( stderr, "%s", msg1 );
    fwrite( opt.name, opt.namelen, 1, stderr );
    fprintf( stderr, "%s", msg2 );
  }

  static option::ArgStatus
  Unknown( const option::Option & option, bool msg )
  {
    if( msg ) printError( "Unknown option '", option, "'\n" );
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus
  Required( const option::Option & option, bool msg )
  {
    if( option.arg != 0 ) return option::ARG_OK;

    if( msg ) printError( "Option '", option, "' requires an argument\n" );
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus
  NonEmpty( const option::Option & option, bool msg )
  {
    if( option.arg != 0 && option.arg[0] != 0 ) return option::ARG_OK;

    if( msg ) printError( "Option '", option, "' requires a non-empty argument\n" );
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus
  Numeric( const option::Option & option, bool msg )
  {
    char * endptr = 0;
    if( option.arg != 0 && strtol( option.arg, &endptr, 10 ) )
    {
    };
    if( endptr != option.arg && *endptr == 0 ) return option::ARG_OK;

    if( msg ) printError( "Option '", option, "' requires a numeric argument\n" );
    return option::ARG_ILLEGAL;
  }
};

enum optionIndex
{
  UNKNOWN,
  HELP,
  OPTIONAL,
  REQUIRED,
  NUMERIC,
  NONEMPTY
};
const option::Descriptor usage[] = {
  { UNKNOWN,
    0,
    "",
    "",
    Arg::Unknown,
    "USAGE: example_arg [options]\n\n"
    "Options:" },
  { HELP, 0, "", "help", Arg::None, "  \t--help  \tPrint usage and exit." },
  { OPTIONAL,
    0,
    "o",
    "optional",
    Arg::Optional,
    "  -o[<arg>], \t--optional[=<arg>]"
    "  \tTakes an argument but is happy without one." },
  { REQUIRED,
    0,
    "r",
    "required",
    Arg::Required,
    "  -r <arg>, \t--required=<arg>  \tMust have an argument." },
  { NUMERIC,
    0,
    "n",
    "numeric",
    Arg::Numeric,
    "  -n <num>, \t--numeric=<num>  \tRequires a number as argument." },
  { NONEMPTY,
    0,
    "1",
    "nonempty",
    Arg::NonEmpty,
    "  -1 <arg>, \t--nonempty=<arg>"
    "  \tCan NOT take the empty string as argument." },
  { UNKNOWN,
    0,
    "",
    "",
    Arg::None,
    "\nExamples:\n"
    "  example_arg --unknown -o -n10 \n"
    "  example_arg -o -n10 file1 file2 \n"
    "  example_arg -nfoo file1 file2 \n"
    "  example_arg --optional -- file1 file2 \n"
    "  example_arg --optional file1 file2 \n"
    "  example_arg --optional=file1 file2 \n"
    "  example_arg --optional=  file1 file2 \n"
    "  example_arg -o file1 file2 \n"
    "  example_arg -ofile1 file2 \n"
    "  example_arg -unk file1 file2 \n"
    "  example_arg -r -- file1 \n"
    "  example_arg -r file1 \n"
    "  example_arg --required \n"
    "  example_arg --required=file1 \n"
    "  example_arg --nonempty= file1 \n"
    "  example_arg --nonempty=foo --numeric=999 --optional=bla file1 \n"
    "  example_arg -1foo \n"
    "  example_arg -1 -- \n"
    "  example_arg -1 \"\" \n" },
  { 0, 0, 0, 0, 0, 0 } };

int
example_arg( int argc, char const * argv[] )
{
  argc -= ( argc > 0 );
  argv += ( argc > 0 );  // skip program name argv[0] if present
  option::Stats stats( usage, argc, argv );

#if 0
    // GCC supports C99 VLAs for C++ with proper constructor calls.
  option::Option options[stats.options_max];
  option::Option buffer[stats.buffer_max];
#else
  // use calloc() to allocate 0-initialized memory. It's not the same
  // as properly constructed elements, but good enough. Obviously in an
  // ordinary C++ program you'd use new[], but this file demonstrates that
  // TLMC++OP can be used without any dependency on the C++ standard library.
  option::Option * options =
    (option::Option *)calloc( stats.options_max, sizeof( option::Option ) );
  option::Option * buffer =
    (option::Option *)calloc( stats.buffer_max, sizeof( option::Option ) );
#endif

  option::Parser parse( usage, argc, argv, options, buffer );

  if( parse.error() ) return 1;

  if( options[HELP] || argc == 0 )
  {
    int columns = getenv( "COLUMNS" ) ? atoi( getenv( "COLUMNS" ) ) : 80;
    option::printUsage( fwrite, stdout, usage, columns );
    return 0;
  }

  for( int i = 0; i < parse.optionsCount(); ++i )
  {
    option::Option & opt = buffer[i];
    fprintf( stdout, "Argument #%d is ", i );
    switch( opt.index() )
    {
      case HELP:
        // not possible, because handled further above and exits the program
      case OPTIONAL:
        if( opt.arg )
          fprintf( stdout, "--optional with optional argument '%s'\n", opt.arg );
        else
          fprintf( stdout, "--optional without the optional argument\n" );
        break;
      case REQUIRED:
        fprintf( stdout, "--required with argument '%s'\n", opt.arg );
        break;
      case NUMERIC:
        fprintf( stdout, "--numeric with argument '%s'\n", opt.arg );
        break;
      case NONEMPTY:
        fprintf( stdout, "--nonempty with argument '%s'\n", opt.arg );
        break;
      case UNKNOWN:
        // not possible because Arg::Unknown returns ARG_ILLEGAL
        // which aborts the parse with an error
        break;
    }
  }

  for( int i = 0; i < parse.nonOptionsCount(); ++i )
    fprintf( stdout, "Non-option argument #%d is %s\n", i, parse.nonOption( i ) );

  return 0;
}

TEST( testOptionParser, example_arg )
{
  int const argc = 3;
  char const * argv[] = { "testOptionParser", "--plus", "--plus" };
  example_arg( argc, argv );
}
