@def hassim=false;
@def title="Linux From Scratch (LFS) virtual machine (VM)"
@def tag=["Linux"]
@def mintoclevel=1

![LFS screenshot](https://fusion809.github.io/images/executor-raujonas.github.io/LFS_screenshot_13-09-2026.png)

**Figure 1: Screenshot of my LFS VM's GNOME session as of 13 September 2026.**

I first installed LFS 12.4 systemd edition to a virtual machine on 9 February 2026. Since then, I have upgraded the system to the development systemd branch, and then gradually made it even more bleeding edge than this by upgrading all packages to the latest stable upstream release. Sometimes I need to keep a package back simply because its latest stable release actually depends on pre-release versions of other packages. This was the case for `gnome-control-center` on 12 September, as at this point version `51.0` of the package was out but it depended on version `51.alpha` or later of `gnome-desktop` and out of these versions only `51.alpha` was available at the time. 

\toc

# Motivations
My motivations for setting us this VM include:
* Curiosity, as I have dozens of free operating system (OS) VMs that I maintain just as a curiosity, so setting up a LFS VM keeps within this. 
* A desire to prove to myself that I can actually run and maintain LFS long term and get it to the point of being a viable daily driver. 

# Package management
From my NixOS host machine, I have written &mdash; with the help of artificial intelligence (AI) &mdash; several shell functions that are imported into my LFS VM and provide basic package management functionality. These functions are part of both my host's and VM's shell profile. These functions can be found in my [NixOS configuration user shell profile](https://github.com/fusion809/NixOS-configs/tree/26.05/shell/user/). [lfs-custom-updates.py](https://github.com/fusion809/NixOS-configs/blob/26.05/python/lfs-custom-updates.py) is used to parallelize and more efficiently check the versions of all custom packages to see if updates are available.

~~~
<table style="border-collapse: collapse; width: 100%;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 1: Shell functions used for package management within the LFS VM.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center; white-space: nowrap;">
            <b>Syntax</b> (definition file hyperlink)
        </td>
        <td style="font-size: 20px; padding: 10px; overflow-wrap: break-word;">
            <b>Description</b>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autobuild-func.sh"><code>autobuild PACKAGE(S) [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Default: build and install the specified package(s), if and only if the latest version of the package is not already installed. LFS/Beyond LFS (BLFS) instructions were used to build most packages. Although, some packages were built using custom build scripts defined in <a href="https://github.com/fusion809/lfs_packaging"><code>~/lfs_packaging</code></a>. These custom build scripts have since become the primary source of software for the system, as it tends to be more reliable and easier to tweak.<br/>
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
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autoremove.sh"><code>autoremove PACKAGE(S) [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Default: remove the specified package(s), if and only if no other packages have libraries that depend on the package(s).<br/>
            Options:<br/>
            <code>--dry-run</code>: show what actions would be executed to remove the package.<br/>
            <code>-f</code>/<code>--force</code>: force removal, without regard for library dependencies.<br/>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-libs.sh"><code>ls_old_libs [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            List old installed versions of libraries.<br/>
            Options:<br/><code>-d</code> option it list files that depend on listed files.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/Shell/02-pms.sh"><code>rm_book_src</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove book source files.
        </td>
    </tr>  
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/Shell/02-pms.sh"><code>rm_lfp_src</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove custom package source tarballs.
        </td>
    </tr>    
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-share.sh"><code>rm_old_docs</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old unused documentation directories.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-kerns.sh"><code>rm_old_kerns</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old unused kernels.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-libs.sh"><code>rm_old_libs</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old unused libraries. As for used old used libraries, rebuild packages that depend on the library and then remove it.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-share.sh"><code>rm_old_share [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Removes old unused <code>/usr/share</code> subdirectories.<br/>
            Options:<br/><code>--dry-run</code> shows what would be done without actually executing those actions.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/Shell/02-pms.sh"><code>rm_src</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old source archives and directories (not including git repos).
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-sync-to-vm.sh"><code>sync_to_vm</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Only available on host; synchronize scripts from host to virtual machine. 
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-update.sh"><code>update [OPTION(s)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Update packages.<br/>
            Options:<br/>
            <code>--dry-run</code>: Show what would be updated without downloading/building.<br/>
            <code>-h</code>/<code>--help</code>: Show help message.<br/>
            <code>--no-upstream</code>: Check only LFS/BLFS book versions (disable upstream tracking).<br/>
            <code>-v</code>/<code>--verbose</code>: Show custom package local and upstream versions during version checking.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-update.sh"><code>updatec [OPTION(s)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Runs, in order, <code>rm_old_docs</code>, <code>rm_old_kerns</code>, <code>rm_old_libs</code>, <code>rm_old_share</code> and <code>rm_src</code> if and only if <code>update</code> runs without error. Options are passed directory to <code>update</code>.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh"><code>updates</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Print table that list available updates (marked with <code>[UPDATE]</code>), as well as packages with missing inventories (marked with <code>[FILES MISSING]</code>) and packages with versioning failures (marked with <code>[FAILED]</code>). 
        </td>
    </tr>
</table>
~~~

# GitHub repositories relating to VM and their locations on VM
* Host system [`NixOS-configs`](https://github.com/fusion809/NixOS-configs/tree/26.05/shell/user/) has shell profile for managing VM, including package management shell functions. Specifically [21-lfs.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh), [lfs-autobuild.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autobuild.sh), [lfs-updates.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh) and [lfs-vm-bootstrap.sh](https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-vm-bootstrap.sh) are the scripts for LFS management. 
* [`~/build_duration`](https://github.com/fusion809/lfs_package_build_times) &mdash; contains text files that contain the duration, in seconds, of building each package.
* [`~/lfs_apps`](https://github.com/fusion809/lfs_apps) &mdash; desktop configuration files and shell scripts these desktop files call. 
* [`~/lfs_dotfiles`](https://github.com/fusion809/lfs_dotfiles) &mdash; Fastfetch, HyFetch and systemd configuration files for LFS VM.
* [`~/lfs_gnuplot`](https://github.com/fusion809/lfs_gnuplot) &mdash; Gnuplot files for my LFS VM.
* [`~/lfs_packaging`](https://github.com/fusion809/lfs_packaging) &mdash; which contains packaging scripts for building custom packages.
* [`~/lfs_scripts`](https://github.com/fusion809/lfs_scripts) &mdash; shell scripts (including VM shell profile and scripts called by Executor and Command Output extensions/widgets) used by LFS system. 
* [`/usr/share/gnome-shell/extensions/executor@raujonas.github.io`](https://github.com/fusion809/executor-raujonas.github.io) &mdash; customized verison of the [`executor@raujonas.github.io`](https://github.com/raujonas/executor) I use under my GNOME session (which is the main session I boot).
* [`/var/lib/book-packages`](https://github.com/fusion809/lfs_book_packages) &mdash; package inventories for LFS and BLFS packages. Now empty as all packages are now provided by custom build scripts. 
* [`/var/lib/custom-packages`](https://github.com/fusion809/lfs_custom_packages) &mdash; package inventories for custom packages (those in `~/lfs_packaging`).

# [`~/build_duration`](https://github.com/fusion809/build_duration)
`~/build_duration` merely contains the logs of how long each completed build has taken. Its contents are created by `~/lfs_scripts/autobuild-log.sh`, which is in turn automatically run by `~/lfs_dotfiles/systemd/user/autobuild-log.service` (which is symlinked to `~/.config/systemd/user/autobuild-log.service`). As for 13 September 2026, it is new and is far from being complete, so many packages' build times are not logged. `build_time pkg` is a shell script function defined in `~/lfs_scripts` that converts the build time into hours, minutes and seconds. `lfs_commit` commits changes made to this repository, along with changes made to `/var/lib/book-packages` and `/var/lib/custom-packages`. 

# [`~/lfs_apps`](https://github.com/fusion809/lfs_apps)
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

# [`~/lfs_dotfiles`](https://github.com/fusion809/lfs_dotfiles)
I have also customized Fastfetch/HyFetch output so that it accurately prints the number of packages I have installed. The Fastfetch configuration file used is located in [`~/lfs_dotfiles/config.jsonc`](https://github.com/fusion809/lfs_dotfiles/blob/master/config.jsonc). The HyFetch configuration files are also in [`~/lfs_dotfiles/hyfetch.json`](https://github.com/fusion809/lfs_dotfiles/blob/master/hyfetch.json). 

In the screenshot above, `838 [ 726,  1,  82,  29]` means that 837 packages are installed in total. Of them 725 were installed via custom build scripts in [`~/lfs_packaging`](https://github.com/fusion809/lfs_packaging). 1 Julia package was installed; this package is Julia itself which was installed via `juliaup` (the compilation process of Julia is incredibly complex and even requires its own custom build of LLVM). 82 Python packages were installed via `pip`. 29 R packages were installed. 

` 576,  498` refers to number of package inventory git repository commits I have published. 576 refers to `/var/lib/book-packages` and 498 refers to `/var/lib/custom-packages`. I include it in Fastfetch output as a way of tracking the versions of custom packages.

There is one systemd service file in [`~/lfs_dotfiles/systemd/user/autobuild-log.service`](https://github.com/fusion809/lfs_dotfiles/blob/master/systemd/user/autobuild-log.service) to autostart [`~/lfs_scripts/autobuild-log.sh`](https://github.com/fusion809/lfs_scripts/blob/master/autobuild-log.sh). 

# [`~/lfs_packaging`](https://github.com/fusion809/lfs_packaging)
`~/lfs_packaging` is presently used to provide all the packages of my LFS virtual machine. It contains directories whose names match the name of the package the `build.sh` script within provides. These scripts cannot be manually executed; instead packages are built using `autobuild pkg` (with the `-f` option required if the latest version of the package is already installed). That being said, `autobuild` does have the capacity to install BLFS and LFS packages from the development book instructions, too, but I prefer the custom package script approach as it allows me to easily edit the build commands and set the version of the package as the latest upstream stable release. Consequently, `autobuild pkg` always defaults to using the `~/lfs_packaging` package when one is available. 

Most of the build instructions in `build.sh` scripts are based on LFS and BLFS build instructions; some are based on SlackBuilds developed, by other packagers, to provide the package for Slackware. Some are also based on Arch Linux PKGBUILDs. `version=` lines in these scripts typically, as their first port of call, will opt to determine the latest upstream version from the source archive website of the package. Failing this, it will use the tags of its git repository. Failing this, it will use [tox-wtf's Version Aggregator and Tracker](https://github.com/tox-wtf/vat). Failing this, it will use Arch Linux's PKGBUILD for the package. Failing this, it will use the LFS or BLFS development book. If all of these methods fail, it will just print the installed version and write to `~/logs/failed_versioning.log` the time and the name of the package whose upstream versioning failed. 

# [`~/lfs_scripts`](https://github.com/fusion809/lfs_scripts)
My shell profile is defined in `~/lfs_scripts`. Some scripts called for by GNOME and KDE Plasma Executor/Command Output commands are in this repository, too. 

# GNOME
GNOME was the first desktop I installed and is the main user interface I use in the virtual machine. My NixOS system uses Hyprland instead, but I have struggled to get Hyprland to actually work in a KVM/QEMU virtual machine, so I decided to just use GNOME in the LFS VM. 

[Dash to Dock](https://github.com/micheleg/dash-to-dock) is enabled and installed, as is [WeatherPanel](https://github.com/attentivecoder/weatherpanel), [Extension List](https://github.com/tuberry/extension-list), [Kiwimenu](https://github.com/kem-a/kiwi-menu), [Show Desktop Button](https://github.com/amivaleo/Show-Desktop-Button) and [Super Into Apps](https://github.com/mikelei8291/super-into-apps). As previously mentioned, I also use my own [own fork](https://github.com/fusion809/executor-raujonas.github.io) of the Executor extension. 

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
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/left_widget_command.sh" target="_blank"><code>~/lfs_scripts/left_widget_command.sh</code></a> &mdash; displays the boot time and age of the system. The age is displayed as days/minutes/years hours:minutes:seconds. In my set up, it is used to generate output for the left widget. Runs every 60ms.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/centre_widget_command.sh" target="_blank"><code>~/lfs_scripts/centre_widget_command.sh</code></a> &mdash; displays CPU, RAM and root filesystem usage percentage and the number of the currently shown wallpaper / the total number of wallpapers in <code>~/wallpapers</code>. Runs every millisecond.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/updates_no.sh" target="_blank"><code>~/lfs_scripts/updates_no.sh</code></a> &mdash; checks for updates using the <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh" target="_blank"><code>updates</code></a> command in the shell profile. It displays <code>$in_progress󰔚 $updates_avg_duration  $mod_time  $no_updates 󰂕 $no_missing_total  $no_failed$failed_version</code> where <code>$in_progress</code> is replaced with nothing if the <code>updates</code> command is not running, and <code>󰦕 ${percent}% </code> otherwise, where <code>$percent</code> is an approximation of how far through the running of <code>updates</code> we are. <code>$updates_avg_duration</code> is the average duration, in minutes and seconds, of the run of <code>updates</code> based on <code>~/logs/updates_duration.log</code>. <code>$mod_time</code> is replaced with the time the <code>updates</code> command last stopped running. <code>$no_updates</code> is replaced with the number of available package updates. <code>$no_missing_total</code> is replaced with the number of packages with missing inventories. <code>$no_failed</code> is replaced with the number of package versioning failures. <code>$failed_version</code>, if <code>~/log/failed_versioning.log</code> is not empty, is replaced by F and the number of packages with version failures in <code>~/logs/failed_versioning.log</code>. <code>updates</code> runs every 5 minutes - the average duration of <code>updates</code> runs. <code>updates_no.sh</code> is run every millisecond.
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
            <code>gnome-terminal -- zsh -ic <a href="https://github.com/fusion809/lfs_scripts/blob/master/list-wallpapers.sh" target="_blank">~/lfs_scripts/list-wallpapers.sh</a></code> &mdash; displays the list of wallpapers in `~/wallpapers` with the currently shown wallpaper highlighted and centred.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "<a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh" target="_blank">updatec</a>; exec zsh"</code> &mdash; updates the system's packages, including those installed via book instructions, custom packages and pip-managed packages and removes unneeded files. 
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
            <code>gnome-terminal -- zsh -ic "source <a href="https://github.com/fusion809/lfs_scripts/blob/master/updates_no_func.sh" target="_blank">~/lfs_scripts/updates_no_func.sh</a>; silent_updates"</code> &mdash; runs <code>updates</code> to update the output shown in the widget.
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
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/open-wallpaper.sh" target="_blank"><code>~/lfs_scripts/open-wallpaper.sh</code></a> &mdash; opens the displayed wallpaper in Eye of GNOME.
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
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/left_widget_tooltip_command.sh" target="_blank"><code>~/lfs_scripts/left_widget_tooltip_command.sh</code></a> &mdash; generates a line describing the version of LFS/BLFS installed, along with the number of packages installed via different means, and package inventory commit numbers in a similar format as in the Fastfetch output. Also includes lines indicating how far into the current run of <code>autobuild &lt;package&gt;</code> the system is. 
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/centre_widget_tooltip_command_wrap.sh" target="_blank"><code>~/lfs_scripts/centre_widget_tooltip_command_wrap.sh</code></a> &mdash; lists selected wallpaper (indicated with <code>></code>) and the 25 wallpapers before and after this one. If there are not 25 wallpapers before the current one, it will show some of the last wallpapers in the collection before the wallpaper numbered 1 to ensure that 51 wallpapers are listed (including the one set as the desktop background). If there are not 25 wallpapers after the current one, it will show some of the first wallpapers in the collection after the final one in the list to ensure that 51 wallpapers are listed in total.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/update-table.sh" target="_blank"><code>~/lfs_scripts/update-table.sh</code></a> &mdash; generates a more compact table of packages with updates, missing inventories and versioning failures.
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
