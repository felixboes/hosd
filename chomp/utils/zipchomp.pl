#!/usr/bin/perl

chdir ("..") unless -d utils;

@which = scalar (@ARGV) ? @ARGV : ('full', 'src', 'lib', 'bin');

foreach my $which (@which)
{
	# Compress the full version of the CHomP software distribution package.
	if ($which eq 'full')
	{
		system ("zip -u -9 chomp-full lib/ bin/ obj/ obj/*/ -x */CVS/");
		system ("zip -u -9 -r chomp-full * -x _* *.zip bin/* doc/html/* obj/* obj/*/* lib/* utils/Doxyfile.warn");
	}

	# Compress the basic version of the CHomP software distribution package.
	if ($which eq 'src')
	{
		system ("zip -u -9 chomp-src lib/ bin/ obj/ obj/*/ -x */CVS/");
		system ("zip -u -9 -r chomp-src makefile licen* include/ make/ src/ programs/chomprog/* programs/chompdemo/* -x src/multiwork/* include/chomp/multiwork/* include/simplices/*");
	}

	# Compress the CHomP library only.
	if ($which eq 'lib')
	{
		system ("zip -u -9 chomp-lib lib/ bin/ obj/ obj/*/ -x */CVS/");
		system ("zip -u -9 -r chomp-lib makefile licen* include/ make/ src/");
	}

	# Compress the CHomP library only, without the CAPD part.
	if ($which eq 'cut')
	{
		system ("zip -u -9 chomp-cut lib/ bin/ obj/ obj/*/ -x */CVS/ obj/capd*/ obj/chompdemo/ obj/chomprog/");
		system ("zip -u -9 -r chomp-cut makefile licen* include/ make/ src/ -x include/capd/* include/capd/*/* src/capd*/* make/auto_dep/capd*");
	}

	# Compress the binary packages.
	if ($which eq 'bin')
	{
		chdir ("bin");
		system ("zip -u -9 ../chomp_deb64 chomp");
		chdir ("..");
		system ("zip -u -9 -r chompfull_deb64 licen* bin/* examples/* python/* -x examples/*.sh");
	}

}
