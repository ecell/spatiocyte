/*
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is a part of E-CELL Simulation Environment package
//
//                Copyright (C) 1996-2001 Keio university
//   Copyright (C) 1998-2001 Japan Science and Technology Corporation (JST)
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-CELL is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// E-CELL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with E-CELL -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//END_HEADER
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//	This file is a part of E-CELL2.
//	Original codes of E-CELL1 core were written by Kouichi TAKAHASHI
//	<shafi@e-cell.org>.
//	Some codes of E-CELL2 core are minor changed from E-CELL1
//	by Naota ISHIKAWA <naota@mag.keio.ac.jp>.
//	Other codes of E-CELL2 core and all of E-CELL2 UIMAN are newly
//	written by Naota ISHIKAWA.
//	All codes of E-CELL2 GUI are written by
//	Mitsui Knowledge Industry Co., Ltd. <http://bio.mki.co.jp/>
//
//	Latest version is availabe on <http://bioinformatics.org/>
//	and/or <http://www.e-cell.org/>.
//END_V2_HEADER
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */
/*
 *::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 *	$Id$
 :	$Log$
 :	Revision 1.10  2004/01/21 06:26:30  shafi
 :	a fix from gabor, remove tmp files when exitting abnormally
 :
 :	Revision 1.9  2003/09/27 13:41:07  bgabor
 :	Bugfix.
 :
 :	Revision 1.6  2003/09/27 12:39:15  satyanandavel
 :	more compatibility issues in Windows
 :
 :	Revision 1.5  2003/09/22 04:28:43  bgabor
 :	Fixed a serious undefined reference to my_open_to_read bug in VVector.
 :
 :	Revision 1.4  2003/07/20 06:05:35  bgabor
 :
 :	added support for large files
 :
 :	Revision 1.3  2003/03/18 09:06:36  shafi
 :	logger performance improvement by gabor
 :
 :	Revision 1.2 2003/03/01 02:43:51  shafi
 :	changed envvar name: ECSTMPDIR to VVECTORTMPDIR, changed defult vvector dir: /var/tmp to /tmp
 :
 :	Revision 1.1  2002/04/30 11:21:53  shafi
 :	gabor's vvector logger patch + modifications by shafi
 :
 :	Revision 1.8  2001/10/21 15:27:49  ishikawa
 :	Warn and exit if temprary directory not exist.
 :
 :	Revision 1.7  2001/10/15 17:18:26  ishikawa
 :	improved program interface
 :
 :	Revision 1.6  2001/07/23 23:04:22  naota
 :	Detailed error message.
 :
 :	Revision 1.5  2001/03/23 18:51:17  naota
 :	comment for credit
 :
 :	Revision 1.4  2001/03/08 17:19:55  naota
 :	cp_stop()
 :
 :	Revision 1.3  2001/03/07 18:28:14  naota
 :	default directory is c:/temp on MS-Windows
 :
 :	Revision 1.2  2001/01/19 21:18:42  naota
 :	Can be comiled g++ on Linux.
 :
 :	Revision 1.1  2001/01/10 11:44:46  naota
 :	batch mode running.
 :
 :	Revision 1.1  2000/12/30 14:57:13  naota
 :	Initial revision
 :
//END_RCS_HEADER
 *::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 */
#include "VVector.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <memory.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "osif.h"


#ifndef O_BINARY
// Must be for UNIX.
#define O_BINARY 0
#endif /* O_BINARY */

#ifndef O_LARGEFILE
// Must be for UNIX.
#define O_LARGEFILE 0
#endif /* O_LARGEFILE */


int vvectorbase::_serialNumber = 0;
char const *vvectorbase::_defaultDirectory = NULL;
int vvectorbase::_directoryPriority = 999;
std::vector<char const *> vvectorbase::_tmp_name;
std::vector<int> vvectorbase::_file_desc_read;
std::vector<int> vvectorbase::_file_desc_write;
bool vvectorbase::_atexitSet = false;
vvectorbase::cbfp_t vvectorbase::_cb_error = NULL;
vvectorbase::cbfp_t vvectorbase::_cb_full = NULL;
long vvectorbase::_margin = 10 * 1024; // by K bytes


static void checkDiskFull(char const * const path, int mustCheck)
{
  const long diskFreeMargin = 10L * 1024L; // K bytes
  // const int checkInterval = 1024;
  const int checkInterval = 10;
  static int skipCounter = 0;

  if (mustCheck != 0 || (checkInterval <= skipCounter)) {
    skipCounter = 0;
    errno = 0;
    long kbytes_free = osif_disk_free(path);
	if (kbytes_free < diskFreeMargin) {
                vvectorbase::cbFull();
    }
  } else {
    skipCounter++;
  }
}


vvectorbase::vvectorbase()
{
  if (_defaultDirectory == NULL) {
    char const *envVal = getenv("VVECTORTMPDIR");
    if (envVal != NULL) {
      _defaultDirectory = strdup(envVal);
      _directoryPriority = 3;
    } else
    {
#ifdef	_Windows
      _defaultDirectory = strdup("c:\\temp");
#else
      _defaultDirectory = strdup("/tmp");
#endif	/* _Windows */
      _directoryPriority = 4;
    }
  }
  _myNumber = _serialNumber;
  _serialNumber++;
  if (!_atexitSet) {
    _atexitSet = true;
    if (atexit(removeTmpFile) != 0) {
      fprintf(stderr, "atexit() fails.  errno=%d\n", errno);
      cbError();
    }
  }

}


vvectorbase::~vvectorbase()
{
unlinkfile();

}

void vvectorbase::unlinkfile()
{
#ifndef OPEN_WHEN_ACCESS
	if (0 <= _fdr) {
		close(_fdr);
//		printf("fdr closed\n");
	}
	if (0 <= _fdw) {
		close(_fdw);
//		printf("fdrw closed\n");
	}
#endif /* OPEN_WHEN_ACCESS */
	if (_file_name != NULL) {
	if (unlink(_file_name) != 0)
		{
    		fprintf(stderr, "unlink(%s) failed in VVector.\n", _file_name);
  		}

		free(_file_name);
	}

}


void vvectorbase::setTmpDir(char const * const dirname, int priority)
{
  assert(dirname != NULL);
  assert(dirname[0] != '\0');
#ifdef DEBUG_DONE
  fprintf(stderr, "vvectorbase::setTmpdir(\"%s\", %d)\n",
	  dirname, priority);
#endif /* DEBUG */
  if (priority < _directoryPriority) {
    _directoryPriority = priority;
    if (_defaultDirectory != NULL) {
      free(const_cast<char*>(_defaultDirectory));
    }
    _defaultDirectory = strdup(dirname);
  }
}


void vvectorbase::removeTmpFile()
{
  std::vector<char const *>::iterator iii;

#ifndef OPEN_WHEN_ACCESS
  std::vector<int>::iterator ii;
  for (ii = _file_desc_read.begin(); ii != _file_desc_read.end(); ii++) {

	if (0 <= *ii) {
		close(*ii);
//		printf("fdr closed\n");
	}
}

  for (ii = _file_desc_write.begin(); ii != _file_desc_write.end(); ii++) {

	if (0 <= *ii) {
		close(*ii);
//		printf("fdrw closed\n");
	}
}
#endif /* OPEN_WHEN_ACCESS */


  for (iii = _tmp_name.begin(); iii != _tmp_name.end(); iii++) {

	unlink (*iii);
 }
}


void vvectorbase::initBase(char const * const dirname)
{
  char pathname[256];
  char filename[256];
  if (dirname != NULL) {
    strcpy(pathname, dirname);
  } else {
    strcpy(pathname, _defaultDirectory);
  }
#if defined(__BORLANDC__) || defined(__WINDOWS__) || defined(__MINGW32__)
  if (pathname[strlen(pathname) - 1] != '\\') {
    strcat(pathname, "\\");
  }
#else
  if (pathname[strlen(pathname) - 1] != '/') {
    strcat(pathname, "/");
  }
#endif
  if (osif_is_dir(pathname) == 0) {
    printf("Directory \"%s\" does not exist.\n", pathname);
    exit(1);
  }
  checkDiskFull(pathname, 1);
  sprintf(filename, "vvector-%ld-%04d",
	  osif_get_pid(), _myNumber);
  strcat(pathname, filename);
  _file_name = strdup(pathname);
  _tmp_name.push_back(_file_name);
  _fdw = open(_file_name, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY |O_LARGEFILE, 0600);
 _fdr = open(_file_name, O_RDONLY | O_BINARY | O_LARGEFILE );
_file_desc_write.push_back(_fdr);
_file_desc_write.push_back(_fdw);
  if (_fdw < 0) {
    fprintf(stderr, "open(\"%s\") failed in VVector.\n", _file_name);
    cbError();
    exit(1);
  }

  if (_fdr < 0) {
    fprintf(stderr, "open(\"%s\") failed in VVector err=%s.\n",
	    _file_name, strerror(errno));
    cbError();
    exit(1);
  }

}


void vvectorbase::my_open_to_append()
{
  checkDiskFull(_file_name, 0);
  _fdw = open(_file_name, O_WRONLY | O_APPEND | O_BINARY |O_LARGEFILE);
  if (_fdw < 0) {
    fprintf(stderr, "open(\"%s\") failed in VVector err=%s.\n",
	    _file_name, strerror(errno));
    cbError();
    exit(1);
  }
}


void vvectorbase::my_open_to_read(off_t offset)
{
if (_fdr<0) {
 _fdr = open(_file_name, O_RDONLY | O_BINARY|O_LARGEFILE );
 }
  if (_fdr < 0) {
    fprintf(stderr, "open(\"%s\") failed in VVector err=%s.\n",
	    _file_name, strerror(errno));
    cbError();
    exit(1);
  }
  if (lseek(_fdr, offset, SEEK_SET) == static_cast<off_t>(-1)) {
    fprintf(stderr, "lseek(\"%s\") failed in VVector err=%s.\n",
	    _file_name, strerror(errno));
    assert(0);
  }
}

void vvectorbase::my_close()
{
  assert((0 <= _fdr)&&(0<=_fdw));
  if ((close(_fdr) < 0)||(close(_fdw) < 0)) {
    fprintf(stderr, "close(\"%s\") failed in VVector err=%s.\n",
	    _file_name, strerror(errno));
    cbError();
    exit(1);
  }
  _fdr = -1;
  _fdw = -1;
}


////////////////////////////////////////////////////////////////////////
//	error handlar
////////////////////////////////////////////////////////////////////////
/* static */ void vvectorbase::cbFull()
{
	if (_cb_full != NULL) {
		(*_cb_full)();
	} else {
		fprintf(stderr,
		  "vvector disk full --- return key to continue.\n");
		getchar();
	}
}


/* static */ void vvectorbase::cbError()
{
	if (_cb_full != NULL) {
		(*_cb_full)();
	} else {
		fprintf(stderr, "error in vvector.\n");
		exit(1);
	}
}


////////////////////////////////////////////////////////////////////////
//	for stand alone test
////////////////////////////////////////////////////////////////////////
#ifdef TEST
#include <math.h>


void	my_full_handler()
{
	printf("my_full_handler() : return key to continue\n");
	getchar();
}


void	my_error_handler()
{
	printf("my_error_handler()\n");
	abort();
}


typedef	struct	{
	double	x;
	double	y;
	double	z;
}	test_data_t;

typedef	vvector<test_data_t>	test_vector_t;


int	main()
{
	vvectorbase::setCBFull(&my_full_handler);
	vvectorbase::setCBError(&my_error_handler);
	vvectorbase::margin(1024 * 100); // by K bytes
	vvectorbase::setTmpDir(".", 1); // top priority

	test_vector_t test_vector;
	for (test_vector_t::size_type iii = 0; iii < 1000; iii++) {
		double		xxx = (double)iii * 0.01;
		test_data_t	test_data;
		test_data.x = xxx;
		test_data.y = sin(xxx);
		test_data.z = cos(xxx);
		test_vector.push_back(test_data);
	}
	for (test_vector_t::size_type iii = 0; iii < 20; iii++) {
		test_data_t	test_data;
		test_data = test_vector.at(iii);
		printf("%d %g %g %g\n",
		  (int)iii, test_data.x, test_data.y, test_data.z);
	}
	return 0;
}


#endif /* TEST */
