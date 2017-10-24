#!/usr/bin/perl/

system("rm OUTPUT_PIRE.NC");
system("rm agg2.exe");
system("pgf90 -o agg.exe -O3 wrapper_PIRE.f90 module_mp_suliaharrington_PIRE_addagg.f90 netcdf_plugin.f90 -lnetcdf -lnetcdff");
system("./agg2.exe");
#open(my $py, "|-", "python agg_box_plots.py") or die "Cannot run Python script: $!";
system("python agg_box_plots_PIRE.py");
#close($py);


exit;

