# `~/lfs_apps`
Many of the desktop configuration files in `~/lfs_apps` generate plots of boot times and cycle through wallpapers.

Plotting files:
* [`plotbts.sh`](https://github.com/fusion809/lfs_apps/blob/master/plotbts.sh) and [`plotbts.desktop`](https://github.com/fusion809/lfs_apps/blob/master/plotbts.desktop) &mdash; boot time histogram with linear scaling on both axes; outliers excluded; not including more recent boots. 
* [`plotbtsa.sh`](https://github.com/fusion809/lfs_apps/blob/master/plotbtsa.sh) and [`plotbtsa.desktop`](https://github.com/fusion809/lfs_apps/blob/master/plotbtsa.desktop) &mdash; boot time histogram with logarithmic scaling on both axes; outliers included; including more recent boots.
* [`plotbtso.sh`](https://github.com/fusion809/lfs_apps/blob/master/plotbtso.sh) and [`plotbtso.desktop`](https://github.com/fusion809/lfs_apps/blob/master/plotbtso.desktop) &mdash; boot time histogram with linear scaling on both axes; outliers included; not including more recent boots.
These rely on [`~/lfs_gnuplot`](https://github.com/fusion809/lfs_gnuplot) Gnuplot code.

Wallpaper cycling files:
* [`cycle-wallpaper.sh`](https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper.sh) and [`cycle-wallpaper.desktop`](https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper.desktop) &mdash; moves us forward through the wallpapers in `~/wallpapers`. Keyboard shortcut: Win+W.
* [`cycle-wallpaper-previous.sh`](https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-previous.sh) and [`cycle-wallpaper-previous.desktop`](https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-previous.desktop) &mdash; moves us backward through the wallpapers in `~/wallpapers`. Keyboard shortcut: Win+Z.
* [`cycle-wallpaper-shuffle.sh`](https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-shuffle.sh) and [`cycle-wallpaper-shuffle.desktop`](https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-shuffle.desktop) &mdash; moves us randomly through the wallpapers in `~/wallpapers`. Keyboard shortcut: Win+S.
* [`specify-wallpaper.sh`](https://github.com/fusion809/lfs_apps/blob/master/specify-wallpaper.sh) and [`specify-wallpaper.desktop`](https://github.com/fusion809/lfs_apps/blob/master/specify-wallpaper.desktop) &mdash; specify the wallpaper (by number) that you want to be set as you desktop background. Keyboard shortcut: Win+N.

Some other desktop configuration files in `~/lfs_apps` open up the settings of GNOME extensions, including:
* [`arcmenu.desktop`](https://github.com/fusion809/lfs_apps/blob/master/arcmenu.desktop) &mdash; for opening the settings for ArcMenu.
* [`dash-to-dock.desktop`](https://github.com/fusion809/lfs_apps/blob/master/dash-to-dock.desktop) &mdash; for opening the settings for Dash to Dock.
* [`executor.desktop`](https://github.com/fusion809/lfs_apps/blob/master/executor.desktop) &mdash; for opening the settings for Executor.
* [`user-theme.desktop`](https://github.com/fusion809/lfs_apps/blob/master/user-theme.desktop) &mdash; for opening the user themes extension settings. 

Some desktop configuration files are instead for displaying pkgs tables:
* [`pkgs_table_alpha.desktop`](https://github.com/fusion809/lfs_apps/blob/master/pkgs_table_alpha.desktop) &mdash; for opening the package table sorted alphabetically.
* [`pkgs_table_bd.desktop`](https://github.com/fusion809/lfs_apps/blob/master/pkgs_table_bd.desktop) &mdash; for opening the package table sorted by build duration.
* [`pkgs_table_size.desktop`](https://github.com/fusion809/lfs_apps/blob/master/pkgs_table_size.desktop) &mdash; for opening the package table sorted by package size. 
Each of these call `pkgs_table_display`, a shell script. 

Other desktop configuration files in `~/lfs_apps` include:
* [`julia.desktop`](https://github.com/fusion809/lfs_apps/blob/master/julia.desktop) &mdash; for running Julia installed via `juliaup`. Without use now that the `julia` package provides Julia instead.