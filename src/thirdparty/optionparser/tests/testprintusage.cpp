/* Written 2012 by Matthias S. Benkmann
 *
 * The author hereby waives all copyright and related rights to the contents
 * of this example file (testprintusage.cpp) to the extent possible under the law.
 */

/**
 * @file
 * @brief Test program for the option::printUsage() function.
 *
 * @note
 * This program is for developing TLMC++OP. It is neither an example nor a functionality test.
 * Do not worry if it doesn't compile or run on your platform.
 */

#include <cstdio>
#include <iostream>
#include <sstream>
#include <unistd.h>

#include "../../optionparser/src/optionparser.h"

using namespace std;
using option::Arg;
using option::Descriptor;

const Descriptor test_vtabs[] = {
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Cölüümn 1 line ı\vColumn 1 line 2\vColumn 1 line 3  \t\vColumn 2 line 2  "
    "\tColumn 3 line 1\v\vColumn 3 line 3  " },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_columns[] = {
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Column 1 line 1  \t\tColumn 3 line 1\n"
    "Column 1 line 2  \tColumn 2 line 2   \tColumn 3 line 2\n"
    "Column 1 line 3  \t\tColumn 3 line 3" },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_column1[] = {
  { 0, 0, "", "", Arg::None, "11 \t21\v22\v23\t 31\nxx" },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_tables[] = {
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Each table has its own column widths and is not aligned with other "
    "tables." },
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Table 1 Column 1 Line 1 \tTable 1 Column 2 Line 1 \tTable 1 Column 3 Line "
    "1\n"
    "Table 1 Col 1 Line 2 \tTable 1 Col 2 Line 2 \tTable 1 Col 3 Line 2" },
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Table 1 Col 1 Line 3 \tTable 1 Col 2 Line 3 \tTable 1 Column 3 Line 3\n"
    "Table 1 Col 1 Line 4 \tTable 1 Column 2 Line 4 \tTable 1 Column 3 Line 4" },
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0, 0, "", "", Arg::None, "This is the only line of table 2." },
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0,
    0,
    "",
    "",
    Arg::None,
    "This is the very long 1st line of table 3. It is more than 80 characters "
    "in length and therefore needs to be wrapped. In fact it is so long that it "
    "needs to be wrapped multiple times to fit into a normal 80 characters "
    "terminal.\v"
    "This is the very long 2nd line of table 3. It is more than 80 characters "
    "in length and therefore needs to be wrapped. In fact it is so long that it "
    "needs to be wrapped multiple times to fit into a normal 80 characters "
    "terminal.\v"
    "This is a reasonably sized line 3 of table 3." },
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Table 4:\n"
    "  \tTable 4 C 1 L 1 \tTable 4 C 2 L 1 \tTable 4 C 3 L 1\n"
    "\tTable 4 C 1 L 2 \tTable 4 C 2 L 2 \tTable 4 C 3 L 2" },
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0, 0, "", "", Arg::None, "This is the only line of table 5" },
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Table 6 C 1 L 1 \tTable 6 C 2 L 1 \tTable 6 C 3 L 1\n"
    "Table 6 C 1 L 2 \tTable 6 C 2 L 2 \tTable 6 C 3 L 2" },
  { 0, 0, "", "", Arg::None, 0 },  // table break
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Table 7 Column 1 Line 1 \tTable 7 Column 2 Line 1 \tTable 7 Column 3 Line "
    "1\n"
    "Table 7 Column 1 Line 2 \tTable 7 Column 2 Line 2 \tTable 7 Column 3 Line "
    "2\n" },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_nohelp[] = { { 0, 0, "", "", Arg::None, 0 },
                                   { 0, 0, "", "", Arg::None, 0 },
                                   { 0, 0, "", "", Arg::None, 0 },
                                   { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_wide[] = { { 0, 0, "", "", Arg::None, "Roma\t|x漢" },
                                 { 0, 0, "", "", Arg::None, "ｶﾀｶﾅ\t|漢字" },
                                 { 0, 0, "", "", Arg::None, "漢字\t|漢ｶ " },
                                 { 0, 0, "", "", Arg::None, "漢字\t|ｶﾅ 漢字" },
                                 { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_overlong[] = {
  { 0, 0, "", "", Arg::None, "Good \t| Good \t| This is good." },
  { 0,
    0,
    "",
    "",
    Arg::None,
    "Good \t| This is an overlong cell. \t| This is good." },
  { 0, 0, "", "", Arg::None, "Good \t| Good \t| This is good." },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_toomanycolumns[] = {
  { 0,
    0,
    "",
    "",
    Arg::None,
    "This \ttable \thas \ttoo \tmany \tcolumns. \tThe \tlast \tcolumns \tare "
    "\tdiscarded." },
  { 0, 0, "", "", Arg::None, "1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11" },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_ownline[] = {
  { 0,
    0,
    "",
    "",
    Arg::None,
    "1234567890AB\vBA0987654321\tStarts on its own line and is indented "
    "somewhat.\vThis one, too." },
  { 0, 0, 0, 0, 0, 0 } };

const Descriptor test_overflow[] = {
  { 0,
    0,
    "",
    "",
    Arg::None,
    "漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢字漢"
    "字漢字漢字" },
  { 0, 0, 0, 0, 0, 0 } };

void
stderr_write( const char * str, int size )
{
  fwrite( str, size, 1, stderr );
}

struct stderr_writer
{
  void
  write( const char * buf, size_t size ) const
  {
    ::write( 2, buf, size );
  }
};

struct stderr_write_functor
{
  void
  operator()( const char * buf, size_t size )
  {
    ::write( 2, buf, size );
  }
};

int
main()
{
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( stderr_write, test_overflow, 1 );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( stderr_write, test_vtabs );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( stderr_writer(), test_columns );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( write, 2, test_column1 );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( cerr, test_tables );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( fwrite, stderr, test_nohelp );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  ostringstream sst;
  option::printUsage( sst, test_wide, 8 );
  cerr << sst.str();
  fputs( "---------------------------------------------------------------\n",
         stderr );
  stderr_write_functor stderr_write_f;
  option::printUsage( &stderr_write_f, test_overlong, 30 );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( stderr_write, test_toomanycolumns );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  option::printUsage( stderr_write, test_ownline, 20 );
  fputs( "---------------------------------------------------------------\n",
         stderr );
  return 0;
}
