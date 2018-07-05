#!/usr/bin/perl/

system("rm COLL.NC");
system("rm coll.exe");
#system("pgf90 -r8 -o coll.exe -O2 netcdf_plugin.f90 collection_lookup_generator.f90 -lnetcdf -lnetcdff");
system("pgf90 -r8 -o coll.exe -O2 collection_lookup_generator_binary.f90");
system("./coll.exe");
#system("python agg_box_plots.py");


exit;

