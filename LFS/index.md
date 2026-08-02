@def hassim=false;
@def title="Linux From Scratch (LFS) virtual machine (VM)"
@def tag=["Linux"]
@def mintoclevel=1

![LFS screenshot](https://fusion809.github.io/images/executor-raujonas.github.io/LFS_screenshot_28-07-2026.png)

**Figure 1: Screenshot of my LFS VM's GNOME session as of 28 July 2026.**

I first installed LFS 12.4 systemd edition to a virtual machine on 9 February 2026. Since then, I have upgraded the system to the development systemd branch, and kept the system up to date. It has been a challenging, yet informative journey.

\toc

# Motivations
My motivations for setting us this VM include:
* Curiosity, as I have dozens of free operating system (OS) VMs that I maintain just as a curiosity, so setting up a LFS VM keeps within this. 
* A desire to prove to myself that I can actually run and maintain LFS long term and get it to the point of being a viable daily driver. 

# GitHub repositories relating to VM and their locations on VM
* Host system [`NixOS-configs`](https://github.com/fusion809/NixOS-configs/tree/26.05/shell/user/) has shell profile for managing VM, including package management shell functions. Specifically [21-lfs.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh), [lfs-autobuild.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autobuild.sh), [lfs-updates.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh) and [lfs-vm-bootstrap.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-vm-bootstrap.sh) are the scripts for LFS management. 
* [`~/lfs_apps`](https://github.com/fusion809/lfs_apps) &mdash; desktop configuration files and shell scripts these desktop files call. 
* [`~/lfs_dotfiles`](https://github.com/fusion809/lfs_dotfiles) &mdash; Fastfetch and HyFetch configuration files for LFS VM.
* [`~/lfs_gnuplot`](https://github.com/fusion809/lfs_gnuplot) &mdash; Gnuplot files for my LFS VM.
* [`~/lfs_packaging`](https://github.com/fusion809/lfs_packaging) &mdash; which contains packaging scripts for building custom packages.
* [`~/lfs-scripts`](https://github.com/fusion809/lfs-scripts) &mdash; shell scripts (including VM shell profile and scripts called by Executor and Command Output extensions/widgets) used by LFS system. 
* [`/var/lib/book-packages`](https://github.com/fusion809/lfs_book-packages) &mdash; package inventories for LFS and BLFS packages.
* [`/var/lib/custom-packages`](https://github.com/fusion809/lfs_custom-packages) &mdash; package inventories for custom packages (`~/lfs_packaging`).

# Package management
From my NixOS host machine, I have written &mdash; with the help of artificial intelligence (AI) &mdash; several shell functions that are imported into my LFS VM and provide basic package management functionality. These functions are part of both my host's and VM's shell profile. These functions can be found in my [NixOS configuration user shell profile](https://github.com/fusion809/NixOS-configs/tree/26.05/shell/user/). 

~~~
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 1: Shell functions used for package management within the LFS VM.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center;" colspan="2">
            <b>Function name</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;" rowspan="2">
            <b>Defined in</b>
        </td>
        <td style="font-size: 20px; padding: 10px;" rowspan="2">
            <b>Syntax</b>
        </td>
        <td style="font-size: 20px; padding: 10px;" rowspan="2">
            <b>Description</b>
        </td>
    </tr>
    <tr>
        <td style="font-size: 20px; padding: 10px;">
            <b>Guest</b>
        </td>
        <td style="font-size: 20px; padding: 10px;">
            <b>Host</b>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>autobuild</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>lfs_autobuild</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autobuild.sh">lfs-autobuild.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>autobuild PACKAGE(S) [OPTION(S)]</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Default: build and install the specified package(s), if and only if the latest version of the package is not already installed. LFS/Beyond LFS (BLFS) instructions are used to build most packages. Although, some packages are built using custom build scripts defined in <a href="https://github.com/fusion809/lfs_packaging"><code>~/lfs_packaging</code></a>.<br/>
            Options:<br/>
            <code>--dry-run</code>: show what actions would be executed to build and install the package.<br/>
            <code>--strip</code>: run stripping commands after build.<br/>
            <code>--no-upstream</code>: disable upstream version searching.<br/>
            <code>--include-config</code>: include configuration commands in the LFS/BLFS book entry.<br/>
            <code>--rm-libs</code>: remove old library versions after build (disabled by default).<br/>
            <code>--lfs</code>: search only in the LFS book.<br/>
            <code>--blfs</code>: search only in the BLFS book.<br/>
            <code>--lfs-book <book></code>: specify LFS book (e.g., development, systemd, stable, or full URL).<br/>
            <code>--blfs-book <book></code>: specify BLFS book (e.g., systemd, development, stable, or full URL).<br/>
            <code>--skip-tests</code>: skip test commands (make check/test, etc.).<br/>
            <code>--ignore-test-failures</code>: ignore test failures by appending '|| true' to test commands.<br/>
            <code>-f</code>/<code>--force</code>: force rebuild and installation even when latest version is already installed.<br/>
            <code>-h</code>/<code>--help</code>: show help message.<br/>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>autoremove</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>lfs_autoremove</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>autoremove PACKAGE(S) [OPTION(S)]</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Default: remove the specified package(s), if and only if no other packages have libraries that depend on the package(s).<br/>
            Options:<br/>
            <code>--dry-run</code>: show what actions would be executed to remove the package.<br/>
            <code>-f</code>/<code>--force</code>: force removal, without regard for library dependencies.<br/>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_book_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_book_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_book_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Remove book source files.
        </td>
    </tr>  
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_lfp_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_lfp_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_lfp_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Remove custom package source tarballs.
        </td>
    </tr>    
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_doc_dirs</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_doc_dirs</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_doc_dirs</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Remove old unused documentation directories.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_kernels</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_kernels</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_kernels</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Remove old unused kernels.
        </td>
    </tr>
        <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_libraries</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_libraries</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_libraries</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Remove old unused libraries. As for used old used libraries, rebuild packages that depend on the library and then remove it.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_share_dirs</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            N/A
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-vm-bootstrap.sh">lfs-vm-bootstrap.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_old_share_dirs [OPTION(S)]</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Removes old unused <code>/usr/share</code> subdirectories.<br/>
            Options:<br/><code>--dry-run</code> shows what would be done without actually executing those actions.
        </td>
    </tr>
        <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>cleanup_src</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Remove old source archives.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>ls_old_libs</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
        N/A
        </td>
        <td style="font-size: 16px; padding: 10px;">
        <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
        <code>ls_old_libs [OPTION(S)]</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
        List old installed versions of libraries.<br/>
        Options:<br/><code>-d</code> option it list files that depend on listed files.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>update</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>lfs_update</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh">21-lfs.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>update [OPTION(s)]</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Update packages.<br/>
            Options:<br/>
            <code>--dry-run</code>: Show what would be updated without downloading/building.<br/>
            <code>--no-upstream</code>: Check only LFS/BLFS book versions (disable upstream tracking).<br/>
            <code>-h</code>/<code>--help</code>: Show help message.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>updates</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>lfs_updates</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh">lfs-updates.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>updates</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Print table that list available updates (marked with <code>[UPDATE]</code>), as well as packages with missing inventories (marked with <code>[FILES MISSING]</code>) and packages with versioning failures (marked with <code>[FAILED]</code>). 
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <code>updatec</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>lfs_updatec</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh">lfs-updates.sh</a>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>updatec</code>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Runs:

            <pre>
            function updatec {
        update
        cleanup_old_doc_dirs
        cleanup_old_kernels
        cleanup_old_libraries
        cleanup_old_share_dirs
        cleanup_src
}</pre>
In other words, it updates all packages, cleans up old libraries, share directories, documentation and source files. 
        </td>
    </tr>
</table>
~~~

# Custom desktop configuration files
The desktop configuration files in [`~/lfs_apps`](https://github.com/fusion809/lfs_apps) generate plots of boot times and cycle through wallpapers.

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

# Custom packages
Some of the packages in [`~/lfs_packaging`](https://github.com/fusion809/lfs_packaging) also have build instructions in LFS and BLFS books &mdash; such as Linux PAM, Vim, rustc and packages within the xorg-apps and xorg-libs metapackages. I provide these custom packages sometimes to overcome build errors that the book-extraction function cause and other times to more robustly ensure I have the latest version of these packages at all times. Other packages are provided in this repository because LFS and BLFS do not provide them; other books in the LFS such as SLFS do provide some of these packages, but some are unique to this repository (e.g. GNU Octave and R are).

# Fastfetch/HyFetch
I have also customized Fastfetch/HyFetch output so that it accurately prints the number of packages I have installed. The Fastfetch configuration file used is located in [`~/lfs_dotfiles/config.jsonc`](https://github.com/fusion809/lfs_dotfiles/blob/master/config.jsonc). The HyFetch configuration files are also in [`~/lfs_dotfiles/hyfetch.json`](https://github.com/fusion809/lfs_dotfiles/blob/master/hyfetch.json). 

In the screenshot above, `813 [ 552,  157,  1,  74,  29]` means that 813 packages are installed in total. Of them 552 are LFS or BLFS book packages installed via `autobuild` and its extracting build commands and source URLs from the books' webpages. A further 157 were installed via custom build scripts in [`~/lfs_packaging`](https://github.com/fusion809/lfs_packaging). 1 Julia package was installed; this package is Julia itself which was installed via `juliaup` (the compilation process of Julia is incredibly complex and even requires its own custom build of LLVM). 74 Python packages were installed via `pip`. 29 R packages were installed. 

` 383,  212` refers to number of package inventory git repository commits I have published. 383 refers to `/var/lib/book-packages` and 212 refers to `/var/lib/custom-packages`. I include it in Fastfetch output as a way of tracking the versions of custom packages.

# Shell profile
My shell profile is defined in [`~/lfs-scripts`](https://github.com/fusion809/lfs-scripts). Some scripts called for by GNOME and KDE Plasma Executor/Command Output commands are in this repository, too. 

# GNOME
GNOME was the first desktop I installed. Its packages are kept at the latest upstream version, not merely the latest version in the BLFS book. [Dash to Dock](https://github.com/micheleg/dash-to-dock) is enabled and installed, as is [WeatherPanel](https://github.com/attentivecoder/weatherpanel), [Extension List](https://github.com/tuberry/extension-list), [Kiwimenu](https://github.com/kem-a/kiwi-menu), [Show Desktop Button](https://github.com/amivaleo/Show-Desktop-Button) and [Super Into Apps](https://github.com/mikelei8291/super-into-apps). [Executor](https://github.com/raujonas/executor) is another extension I use; I've actually created my [own fork](https://github.com/fusion809/executor-raujonas.github.io) with more features.

~~~
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 2: GNOME themes.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Cursor</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Legacy applications</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Icons</b>
        </td>
        <td style="font-size: 20px; padding: 10px;">
            <b>Shell</b>
        </td>
    </tr>
    <tr>
    <td style="font-size: 16px; padding: 10px;">
    WhiteSur-cursors
    </td>
    <td style="font-size: 16px; padding: 10px;">WhiteSur-Dark
    </td>
    <td style="font-size: 16px; padding: 10px;">WhiteSur-dark
    </td>
    <td style="font-size: 16px; padding: 10px;">TST - Semi Transparent
    </td>
    </tr>
</table>
~~~

## Executor fork
*Has `~/lfs_packaging` package called [executor](https://github.com/fusion809/lfs_packaging/tree/master/executor). It can be installed via more standard ways, too.*

The base [Executor](https://github.com/raujonas/executor) extension provides up to three widgets in the GNOME panel on the left, centre and right of the panel. In these widgets is displayed the output of specified commands. The interval at which the command is re-run can also be specified. The [Executor fork](https://github.com/fusion809/executor-raujonas.github.io) I maintain provides the following additional features:
* Tooltips &mdash; which can have two components. They are, in order: (1) static text and (2) command output.
* Command execution when the panel widget is clicked, with separate commands for left-, middle-, and right-click actions. 

~~~
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 3: my Executor extension settings.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Field</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Left widget</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Centre widget</b>
        </td>
        <td style="font-size: 20px; padding: 10px;">
            <b>Right widget</b>
        </td>
    </tr>
    <tr>
    <td style="font-size: 16px; padding: 10px;">
    <b>Index in widget</b>
    </td>
    <td style="font-size: 16px; padding: 10px;">3
    </td>
    <td style="font-size: 16px; padding: 10px;">2
    </td>
    <td style="font-size: 16px; padding: 10px;">2
    </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Output command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/left_widget_command.sh" target="_blank"><code>~/lfs-scripts/left_widget_command.sh</code></a> &mdash; displays the boot time of the system. In my set up, it is used to generate output for the left widget. Runs every 60s.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/count-wallpapers.sh" target="_blank"><code>~/lfs-scripts/count-wallpapers.sh</code></a> &mdash; displays the number of the currently shown wallpaper / the total number of wallpapers in <code>~/wallpapers</code>. Runs every second.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/updates_no.sh" target="_blank"><code>~/lfs-scripts/updates_no.sh</code></a> &mdash; checks for updates using the <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh" target="_blank"><code>updates</code></a> command in the shell profile. It displays <code>$in_progress $mod_time  $no_updates 󰂕 $no_missing_total  $no_failed</code> where <code>$in_progress</code> is replaced with nothing if the <code>updates</code> command is not running, and <code>󰦕 ${percent}% </code> otherwise, where <code>$percent</code> is an approximation of how far through the running of <code>updates</code> we are. <code>$mod_time</code> is replaced with the time the <code>updates</code> command last stopped running. <code>$no_updates</code> is replaced with the number of available package updates. <code>$no_missing_total</code> is replaced with the number of packages with missing inventories. <code>$no_failed</code> is replaced with the number of package versioning failures. <code>updates</code> runs every 5 minutes - the average duration of <code>updates</code> runs. <code>updates_no.sh</code> is run every second.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Left-click command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-previous.sh" target="_blank"><code>~/lfs_apps/cycle-wallpaper-previous.sh</code></a> &mdash; show previous wallpaper.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic <a href="https://github.com/fusion809/lfs-scripts/blob/master/list-wallpapers.sh" target="_blank">~/lfs-scripts/list-wallpapers.sh</a></code> &mdash; displays the list of wallpapers in `~/wallpapers` with the currently shown wallpaper highlighted and centred.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "<a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh" target="_blank">update</a>; exec zsh"</code> &mdash; updates the system's packages, including those installed via book instructions, custom packages and pip-managed packages. 
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Middle-click command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-shuffle.sh" target="_blank"><code>~/lfs_apps/cycle-wallpaper-shuffle.sh</code></a> &mdash; show a randomly-selected wallpaper.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-extensions prefs executor@raujonas.github.io</code> &mdash; opens the settings dialog for Executor.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "source <a href="https://github.com/fusion809/lfs-scripts/blob/master/updates_no_func.sh" target="_blank">~/lfs-scripts/updates_no_func.sh</a>; silent_updates"</code> &mdash; runs <code>updates</code> to update the output shown in the widget.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Right-click command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper.sh" target="_blank"><code>~/lfs_apps/cycle-wallpaper.sh</code></a> &mdash; show next wallpaper.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/open-wallpaper.sh" target="_blank"><code>~/lfs-scripts/open-wallpaper.sh</code></a> &mdash; opens the displayed wallpaper in Eye of GNOME.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "tail -f ~/updates.log"</code> &mdash; opens a terminal and follows the output of the <code>updates</code> command being used to generate the widget content.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Tooltip text</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Left click: previous wallpaper (Win+Z).<br/>
Middle click: shuffle wallpaper (Win+S).<br/>
Right click: next wallpaper (Win+W).<br/>
Win+N: show wallpaper whose number you will be asked to specify.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Left click: list wallpapers with displayed wallpaper centred and highlighted.<br/>
Middle click: open Executor settings (Win+E).<br/>
Right click: open wallpaper in EOG.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Left click: run `update`.<br/>
Middle click: update notifications.<br/>
Right click: show log of last update check.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Tooltip command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/left_widget_tooltip_command.sh" target="_blank"><code>~/lfs-scripts/left_widget_tooltip_command.sh</code></a> &mdash; generates a line describing the version of LFS/BLFS installed, along with the number of packages installed via different means, and package inventory commit numbers in a similar format as in the Fastfetch output. 
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/centre_widget_tooltip_command_wrap.sh" target="_blank"><code>~/lfs-scripts/centre_widget_tooltip_command_wrap.sh</code></a> &mdash; lists selected wallpaper (indicated with `>`) and the 25 wallpapers before and after this one. If there are not 25 wallpapers before the current one, it will show some of the last wallpapers in the collection before the wallpaper numbered 1 to ensure that 51 wallpapers are listed (including the one set as the desktop background). If there are not 25 wallpapers after the current one, it will show some of the first wallpapers in the collection after the final one in the list to ensure that 51 wallpapers are listed in total.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs-scripts/blob/master/update-table.sh" target="_blank"><code>~/lfs-scripts/update-table.sh</code></a> &mdash; generates a more compact table of packages with updates, missing inventories and versioning failures.
        </td>
    </tr>
</table>
~~~

~~~
<br/>
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 4: Example tooltip contents.</b></caption>
    <tr>
    <td style="font-size: 20px; padding: 10px; text-align: center;">
    <b>Left</b>
    </td>
    <td style="font-size: 20px; padding: 10px; text-align: center;">
    <b>Centre</b>
    </td>
    <td style="font-size: 20px; padding: 10px; text-align: center;">
    <b>Right</b>
    </td>
    </tr>
    <tr>
    <td style="font-size: 16px; padding: 10px;">
    <img src="https://fusion809.github.io/images/executor-raujonas.github.io/Left_tooltip.png"/>
    </td>
    <td style="font-size: 16px; padding: 10px;">
    <img src="https://fusion809.github.io/images/executor-raujonas.github.io/Centre_tooltip.png"/>
    </td>
    <td style="font-size: 16px; padding: 10px;">
    <img src="https://fusion809.github.io/images/executor-raujonas.github.io/Right_tooltip.png"/>
    </td>
    </tr>
</table>
~~~

## Installing extensions via my web browser
BLFS did not provide a `gnome-browser-connector` package, which is required for installing GNOME extensions within one's browser. Manually compiling and installing it was fairly easy, however. That being said, whenever I tried to install an extension using it, I noticed that the extension was not successfully installed despite there being a folder in `~/.local/share/gnome-shell/extensions` for it. As this folder would be completely empty. Why? Well, running strace on GNOME shell revealed the problem was actually that the `gnome-browser-connector` was running `unzip` commands that assumed that Info-ZIP's unzip command was installed, not the `bsdunzip` variety provided by libarchive (which is the only one provided by BLFS or LFS). 

I tried compiling Info-ZIP's unzip, such as by following some [old BLFS instructions](https://www.linuxfromscratch.org/blfs/view/cvs/general/unzip.html) but this failed as Info-ZIP's unzip has not been updated since ~2009 and requires multiple intricate patches to get it to compile. The consolidated patch provided by BLFS was not even sufficient, even after I located the patch (the link provided in the book entry shared is actually dead, so I had to find a link to the patch elsewhere by Googling). 

Luckily, ChatGPT provided a script version of `unzip` that would run `bsdtar` in the background and could take all the arguments that `gnome-browser-connector` provided it. I have since included this script in my [custom package for libarchive](https://github.com/fusion809/lfs_packaging/tree/master/libarchive).

# KDE Plasma
KDE Plasma was the second desktop I installed. [Panel Spacer Extended](https://github.com/luisbocanegra/plasma-panel-spacer-extended) extension is installed, as is the [Command Output](https://github.com/Zren/plasma-applet-commandoutput) Plasma widget. 