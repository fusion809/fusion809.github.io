# `~/lfs_apps`
The desktop configuration files in `~/lfs_apps` generate plots of boot times and cycle through wallpapers.

Plotting files:
* `plotbts.sh` and `plotbts.desktop` &mdash; boot time histogram with linear scaling on both axes; outliers excluded; not including more recent boots. 
* `plotbtsa.sh` and `plotbtsa.desktop` &mdash; boot time histogram with logarithmic scaling on both axes; outliers included; including more recent boots.
* `plotbtso.sh` and `plotbtso.desktop` &mdash; boot time histogram with linear scaling on both axes; outliers included; not including more recent boots.
These rely on [`~/lfs_gnuplot`](https://github.com/fusion809/lfs_gnuplot) Gnuplot code.

Wallpaper cycling files:
* `cycle-wallpaper.sh` and `cycle-wallpaper.desktop` &mdash; moves us forward through the wallpapers in `~/wallpapers`. Keyboard shortcut: Win+W.
* `cycle-wallpaper-previous.sh` and `cycle-wallpaper-previous.desktop` &mdash; moves us backward through the wallpapers in `~/wallpapers`. Keyboard shortcut: Win+Z.
* `cycle-wallpaper-shuffle.sh` and `cycle-wallpaper-shuffle.desktop` &mdash; moves us randomly through the wallpapers in `~/wallpapers`. Keyboard shortcut: Win+S.
* `specify-wallpaper.sh` and `specify-wallpaper.desktop` &mdash; specify the wallpaper (by number) that you want to be set as you desktop background. Keyboard shortcut: Win+N.