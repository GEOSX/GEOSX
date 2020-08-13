/* Written 2012 by Matthias S. Benkmann
 *
 * The author hereby waives all copyright and related rights to the contents
 * of this example file (testparse.cpp) to the extent possible under the law.
 */

/**
 * @file
 * @brief Test program for option::Stats and option::Parser.
 *
 * @note
 * This program is for developing TLMC++OP. It is neither an example nor a functionality test.
 * Do not worry if it doesn't compile or run on your platform.
 *
 */

/**
 * @mainpage
 * @copydetails optionparser.h
 */

#include "gtest/gtest.h"

#include <assert.h>
#include <stdio.h>

#include "../../optionparser/src/optionparser.h"

using option::ArgStatus;
using option::Descriptor;
using option::Option;
using option::Parser;
using option::Stats;

struct Arg : public option::Arg
{
  static ArgStatus
  Required( const Option & option, bool )
  {
    return option.arg == 0 ? option::ARG_ILLEGAL : option::ARG_OK;
  }
  static ArgStatus
  Empty( const Option & option, bool )
  {
    return ( option.arg == 0 || option.arg[0] == 0 ) ? option::ARG_OK
                                                     : option::ARG_IGNORE;
  }
};

char *
gettext( const char * msgid )
{
  return (char *)msgid;
}

const Descriptor empty_usage[] = { { 0, 0, 0, 0, 0, 0 } };

const Descriptor minimal_usage[] =    //
  { { 0, 0, "x", "", Arg::None, 0 },  //
    { 0, 0, 0, 0, 0, 0 } };

const Descriptor optional_usage[] =       //
  { { 0, 0, "f", "", Arg::Required, 0 },  //
    { 0, 0, 0, 0, 0, 0 } };

const Descriptor gettext_usage[] =                                  //
  { { 0, 0, "f", "", Arg::Required, gettext( "This is a test" ) },  //
    { 0, 0, 0, 0, 0, 0 } };

enum OptionIndex
{
  UNKNOWN,
  FOO,
  VERBOSE,
  X,
  ABBREVIATE,
  EMPTY
};
enum OptionType
{
  UNUSED = 0,
  DISABLED = 1,
  ENABLED = 2
};

const Descriptor multi_usage[] =                           //
  { { UNKNOWN, 0, "", "", Arg::None, 0 },                  // UNKNOWN option catch-all
    { FOO, ENABLED, "", "enable-foo", Arg::None, 0 },      // FOO enable
    { FOO, DISABLED, "", "disable-foo", Arg::None, 0 },    // FOO disable
    { VERBOSE, 0, "v", "verbose", Arg::None, 0 },          // VERBOSE (counted option)
    { X, 0, "X", "X", Arg::Required, 0 },                  // -X<arg>, -X <arg>, -X=<arg>, --X=<arg>
    { ABBREVIATE, 0, "", "abbreviate-me", Arg::None, 0 },  // ABBREVIATE
    { EMPTY, 0, "", "emptyarg", Arg::Empty, 0 },           // EMPTY (ignores arguments that are not "")
    { 0, 0, 0, 0, 0, 0 } };

const char * empty_args[] = { 0 };
const char * non_options[] = { "1", "2", "3", (const char *)-1 };
const char * unknown_option[] = { "--unknown", "nonoption", 0 };
const char * lone_minus[] = { "-f", "-", "-", 0 };
const char * lone_doubleminus[] = { "--", 0 };

// NOTE: (const char*) -1 is used to cause a segfault should this element be dereferenced.
// 0 is not used here, because the parser explicitly checks for 0 which could mask bugs.
// If a number of arguments >= 0 is passed, the parser is supposed to honor that and never
// dereference an element beyond the last.
const char * multi1[] = { "--enable-foo",
                          "--unknown1",
                          "-u",
                          "-vX",
                          "xyzzy",
                          "--",
                          "--strangefilename",
                          (const char *)-1 };
const char * multi2[] = { "-vvXfoo",
                          "-X",
                          "bar",
                          "-X=foobar",
                          "-X",
                          "",
                          "--disable-foo",
                          "-v",
                          (const char *)-1 };
const char * multi3[] = { "-abbr",
                          "-abb",
                          "--emptyarg",
                          "-verbose",
                          "--emptyarg",
                          "",
                          "--emptyarg=",
                          "nonoption1",
                          "nonoption2",
                          (const char *)-1 };

const char * illegal[] = { "-X", 0 };
const char * reorder[] =
  { "-X", "--", "-", "-X", "--", "foo", "-v", "--", "bar", "--", 0 };
const char * reorder2[] = { "-X", "--", "-", "-X", "--", "-", 0 };

int
count( const char ** args )
{
  for( int c = 0;; ++c )
    if( args[c] == (const char *)-1 ) return c;
}

bool
eq( const char * s1, const char * s2 )
{
  if( s1 == s2 ) return true;

  if( s1 == 0 || s2 == 0 ) return false;

  while( *s1 != 0 && *s2 != 0 )
  {
    ++s1;
    ++s2;
  }

  return *s1 == *s2;
}

TEST( testparse, testparse )
{
  {
    Stats stats( empty_usage, -1, empty_args );
    stats.add( empty_usage, 0, empty_args );
    assert( stats.buffer_max == 1 );
    assert( stats.options_max == 1 );
    Option buffer[1];
    Option options[1];
    Parser parse( empty_usage, 99, empty_args, options, buffer );
    parse.parse( empty_usage, -1, empty_args, options, buffer );
    assert( parse.optionsCount() == 0 );
    assert( parse.nonOptionsCount() == 0 );
    assert( !buffer[0] );
    assert( !options[0] );
    assert( buffer[0].count() == 0 );
    assert( parse.nonOptions() == 0 );

    stats.add( empty_usage, 3, non_options );
    assert( stats.buffer_max == 1 );
    assert( stats.options_max == 1 );
    parse.parse( empty_usage, 3, non_options, options, buffer );
    assert( parse.optionsCount() == 0 );
    assert( parse.nonOptionsCount() == 3 );
    assert( !buffer[0] );
    assert( !options[0] );
    assert( parse.nonOptions() == &non_options[0] );

    stats.add( minimal_usage, -1, unknown_option );
    assert( stats.buffer_max == 1 );
    assert( stats.options_max == 2 );
    parse.parse( minimal_usage, -1, unknown_option, options, buffer );
    assert( parse.optionsCount() == 0 );
    assert( parse.nonOptionsCount() == 1 );
    assert( !buffer[0] );
    assert( !options[0] );
    assert( parse.nonOptions() == &unknown_option[1] );
  }
  {
    Stats stats( gettext_usage, -1, lone_minus );
    Stats stats2;
    stats2.add( gettext_usage, -1, lone_minus );
    assert( stats.buffer_max == 2 );
    assert( stats.options_max == 2 );
    assert( stats2.buffer_max == 2 );
    assert( stats2.options_max == 2 );
    Option buffer[2];
    Option options[2];
    Parser parse;
    parse.parse( gettext_usage, -1, lone_minus, options, buffer );
    assert( parse.optionsCount() == 1 );
    assert( parse.nonOptionsCount() == 1 );
    assert( parse.nonOptions() == &lone_minus[2] );
    assert( options[0] );
    assert( buffer[0] );
    assert( options[0].count() == 1 );
    assert( options[0].isFirst() );
    assert( options[0].isLast() );
    assert( options[0].first() == options[0] );
    assert( options[0].last() == options[0] );
    assert( options[0].prevwrap() == &options[0] );
    assert( options[0].nextwrap() == &options[0] );
    assert( options[0].prev() == 0 );
    assert( options[0].next() == 0 );
    assert( options[0].desc == &gettext_usage[0] );
    assert( eq( options[0].name, "f" ) );
    assert( eq( options[0].arg, "-" ) );
  }
  {
    Stats stats( optional_usage, -1, lone_minus );
    Stats stats2;
    stats2.add( optional_usage, -1, lone_minus );
    assert( stats.buffer_max == 2 );
    assert( stats.options_max == 2 );
    assert( stats2.buffer_max == 2 );
    assert( stats2.options_max == 2 );
    Option buffer[2];
    Option options[2];
    Parser parse;
    parse.parse( optional_usage, -1, lone_minus, options, buffer );
    assert( parse.optionsCount() == 1 );
    assert( parse.nonOptionsCount() == 1 );
    assert( parse.nonOptions() == &lone_minus[2] );
    assert( options[0] );
    assert( buffer[0] );
    assert( options[0].count() == 1 );
    assert( options[0].isFirst() );
    assert( options[0].isLast() );
    assert( options[0].first() == options[0] );
    assert( options[0].last() == options[0] );
    assert( options[0].prevwrap() == &options[0] );
    assert( options[0].nextwrap() == &options[0] );
    assert( options[0].prev() == 0 );
    assert( options[0].next() == 0 );
    assert( options[0].desc == &optional_usage[0] );
    assert( eq( options[0].name, "f" ) );
    assert( eq( options[0].arg, "-" ) );
  }
  {
    Stats stats;
    stats.add( minimal_usage, -1, lone_doubleminus );
    assert( stats.buffer_max == 1 );
    assert( stats.options_max == 2 );
    Option buffer[1];
    Option options[2];
    Parser parse( minimal_usage, -1, lone_doubleminus, options, buffer );
    assert( parse.optionsCount() == 0 );
    assert( parse.nonOptionsCount() == 0 );
    assert( !buffer[0] );
    assert( !options[0] );
    assert( parse.nonOptions() == 0 );
  }
  {
    Stats stats;
    stats.add( multi_usage, count( multi1 ), multi1, 4, true );
    assert( stats.buffer_max == 6 );
    assert( stats.options_max == 7 );
    stats.add( multi_usage, count( multi2 ), multi2, 4, true );
    assert( stats.buffer_max == 14 );
    assert( stats.options_max == 7 );
    stats.add( multi_usage, count( multi3 ), multi3, 4, true );
    assert( stats.buffer_max == 22 );
    assert( stats.options_max == 7 );
    Option buffer[22];
    Option options[7];
    assert( options[FOO].last()->type() == UNUSED );
    assert( options[ABBREVIATE].count() == 0 );
    Parser parse;
    assert( !parse.error() );

    parse.parse( multi_usage, count( multi1 ), multi1, options, buffer, 4, true );
    assert( !parse.error() );
    assert( parse.optionsCount() == 5 );
    assert( parse.nonOptionsCount() == 1 );
    assert( eq( parse.nonOptions()[0], "--strangefilename" ) );
    assert( options[FOO].last()->type() == ENABLED );
    assert( eq( options[FOO].last()->name, "--enable-foo" ) );
    assert( options[FOO].last()->arg == 0 );
    assert( options[UNKNOWN].count() == 2 );
    assert( eq( options[UNKNOWN].first()->name, "--unknown1" ) );
    assert( eq( options[UNKNOWN].last()->name, "u" ) );
    assert( options[UNKNOWN].first()->arg == 0 );
    assert( options[UNKNOWN].last()->arg == 0 );
    assert( options[VERBOSE].count() == 1 );
    assert( options[VERBOSE].arg == 0 );
    assert( options[VERBOSE].name[0] == 'v' && options[VERBOSE].namelen == 1 );
    assert( eq( options[X].arg, "xyzzy" ) );
    assert( eq( options[X].name, "X" ) );

    parse.parse( multi_usage, count( multi2 ), multi2, options, buffer, 4, true );
    assert( !parse.error() );
    assert( parse.optionsCount() == 13 );
    assert( parse.nonOptionsCount() == 1 );
    assert( eq( parse.nonOptions()[0], "--strangefilename" ) );
    assert( options[FOO].last()->type() == DISABLED );
    assert( options[FOO].last()->arg == 0 );
    assert( options[UNKNOWN].count() == 2 );
    assert( eq( options[UNKNOWN].first()->name, "--unknown1" ) );
    assert( eq( options[UNKNOWN].last()->name, "u" ) );
    assert( options[VERBOSE].count() == 4 );
    assert( options[X].count() == 5 );
    const char * Xargs[] = { "xyzzy", "foo", "bar", "foobar", "", "sentinel" };
    const char ** Xarg = &Xargs[0];
    for( Option * Xiter = options[X]; Xiter != 0; Xiter = Xiter->next() )
      assert( eq( Xiter->arg, *Xarg++ ) );

    assert( !options[ABBREVIATE] );
    parse.parse( multi_usage, count( multi3 ), multi3, options, buffer, 4, true );
    assert( !parse.error() );
    assert( parse.optionsCount() == 21 );
    assert( parse.nonOptionsCount() == 2 );
    assert( eq( parse.nonOptions()[0], "nonoption1" ) );
    assert( eq( parse.nonOptions()[1], "nonoption2" ) );
    assert( options[ABBREVIATE] );
    assert( options[EMPTY].count() == 3 );
    assert( options[EMPTY].first()->arg == 0 );
    assert( eq( options[EMPTY].last()->arg, "" ) );
    assert( eq( options[EMPTY].last()->prev()->arg, "" ) );
    assert( options[FOO].last()->type() == DISABLED );
    assert( options[UNKNOWN].count() == 5 );
    assert( eq( options[UNKNOWN].first()->name, "--unknown1" ) );
    assert( options[UNKNOWN].first()->arg == 0 );
    assert( eq( options[UNKNOWN].last()->name, "b" ) );
    assert( options[VERBOSE].count() == 5 );
    assert( options[X].count() == 5 );
    Xarg = &Xargs[0];
    for( Option * Xiter = options[X]; Xiter != 0; Xiter = Xiter->next() )
      assert( eq( Xiter->arg, *Xarg++ ) );

    for( Option * opt = buffer[0]; *opt; ++opt )
      if( opt->desc->check_arg != Arg::Required &&
          opt->desc->check_arg != Arg::Empty )
        assert( opt->arg == 0 );
  }
  {
    Option buffer[2];
    Option options[20];
    Parser parse;
    assert( !parse.error() );
    parse.parse( multi_usage, -1, illegal, options, buffer, 0, false, 2 );
    assert( parse.error() );
  }
  {
    Stats stats( multi_usage, count( multi3 ), multi3, 0, true );
    const int bufmax = 3;
    Option buffer[bufmax];
    Option options[100];
    assert( !options[ABBREVIATE] );
    Parser parse( multi_usage, count( multi3 ), multi3, options, buffer, 4, true, bufmax );
    assert( !parse.error() );
    assert( parse.optionsCount() == bufmax );
    assert( parse.nonOptionsCount() == 2 );
    assert( eq( parse.nonOptions()[0], "nonoption1" ) );
    assert( eq( parse.nonOptions()[1], "nonoption2" ) );
    assert( options[ABBREVIATE] );
    assert( options[UNKNOWN].count() ==
            2 );  // because of buxmax the 2nd 'b' cannot be stored
    assert( options[UNKNOWN].first()->name[0] == 'a' &&
            options[UNKNOWN].first()->namelen == 1 );
    assert( options[UNKNOWN].first()->arg == 0 );
    assert( eq( options[UNKNOWN].last()->name, "bb" ) );
  }
  {
    Stats stats( true, multi_usage, -1, reorder );
    Option buffer[100];
    Option options[100];
    Parser parse( true, multi_usage, -1, reorder, options, buffer );
    assert( !parse.error() );
    assert( parse.optionsCount() == 3 );
    assert( parse.nonOptionsCount() == 4 );
    assert( parse.nonOptions() == &reorder[6] );
  }
  {
    Option buffer[10];
    Option options[10];
    Parser parse( true, multi_usage, 666, reorder2, options, buffer, 0, false, 10 );
    assert( !parse.error() );
    assert( parse.optionsCount() == 2 );
    assert( parse.nonOptionsCount() == 2 );
  }

  fprintf( stdout, "All tests passed.\n" );
}
