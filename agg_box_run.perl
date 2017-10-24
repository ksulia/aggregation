#!/usr/bin/perl/

system("rm OUTPUT.NC");
system("rm agg.exe");
system("pgf90 -o agg.exe -O3 wrapper.f90 module_mp_suliaharrington.f90 netcdf_plugin.f90 -lnetcdf -lnetcdff");
system("./agg.exe");
#open(my $py, "|-", "python agg_box_plots.py") or die "Cannot run Python script: $!";
system("python agg_box_plots.py");
#close($py);


exit;

